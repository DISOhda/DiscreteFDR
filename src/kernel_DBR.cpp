#include "kernel.h"

NumericVector kernel_DBR_fast(const List &pCDFlist, const NumericVector &sorted_pv, const double &lambda, const Nullable<NumericVector> pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  // array to store F_i
  NumericVector *sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // vector to store transformed p-values
  NumericVector pval_transf(numTests, 1.0);
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numTests - 1) / size + 1;
  // last positions in step function evaluations
  int *pos = new int[numCDF]();
  // order of current p-value chunk
  IntegerVector ord;
  
  for(int i = 0; i < chunks; i++) {
    checkUserInterrupt();
    // the min( , numTests) is here for the last chunk
    NumericVector pv = sorted_pv[Range(i * size, std::min<int>((i + 1) * size, numTests) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numCDF; j++) {
      int len = sfuns[j].length();
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], len, pos[j]);
    }
    
    int j;
    for(j = 0; j < numPV; j++) {
      checkUserInterrupt();
      // index of current value in 'sorted_pv'
      // is also the number of values to be *left out* of the current sum!
      int idx_pval = i * size + j;
      // store column values in a vector
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // get order
      if(pCDFcounts.isNull()) {
        std::sort(pv.begin(), pv.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(pv, true);
      // this p-value must satisfy condition; else stop computations
      if(pv(ord[0]) <= lambda) {
        // number of remaining needed values
        int rem = numTests - idx_pval;
        // compute sum
        int k = 0;
        pval_transf[idx_pval] = 0;
        while(k < numCDF && CDFcounts[ord[k]] < rem) {
          pval_transf[idx_pval] += CDFcounts[ord[k]] * pv[ord[k]];
          rem -= CDFcounts[ord[k]];
          k++;
        }
        pval_transf[idx_pval] += rem * pv[ord[k]];
        pval_transf[idx_pval] /= (1 - lambda) * (idx_pval + 1);
      } else break;
    }
    // no more p-values will satisfy condition => stop computations
    if(j < numPV && mat(ord[0], j) > lambda) break;
  }
  
  // garbage collection
  delete[] pos;
  delete[] sfuns;
  
  return pval_transf;
}

List kernel_DBR_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const double &lambda, const double &alpha, const Nullable<NumericVector> pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  // array to store F_i
  NumericVector *sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // support size
  int numValues = support.length();
  // apply the shortcut drawn from Corollary 3, that is
  // c.1 >= the effective critical value associated to min((1 - lambda) * alpha/numTests , lambda)
  //int i = numValues - 1;
  //while(i > 0 && support[i] >= std::min<double>(lambda, (1 - lambda) * alpha / numTests)) i--;
  int i = binary_search(support, std::min<double>(lambda, (1 - lambda) * alpha / numTests), numValues);
  NumericVector pv_list = support[Range(i, numValues - 1)];
  pv_list = sort_combine(pv_list, sorted_pv);
  
  // number of p-values to be transformed
  numValues = pv_list.length();
  // critical values indices
  NumericVector crit(numTests);
  // transformed observed p-values
  NumericVector pval_transf(numTests, 1.0);
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current value in 'pv_list'
  int idx_pval = 0;
  // index of current critical value
  int idx_crit = 0;
  // index of current observed p-value to be transformed
  int idx_transf = 0;
  
  // last positions in step function evaluations
  int *pos = new int[numCDF]();
  // order of current p-value chunk
  IntegerVector ord;
  
  for(int i = 0; i < chunks; i++) {
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numCDF; j++) {
      int len = sfuns[j].length();
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], len, pos[j]);
    }
    
    // loop variable
    int j = 0;
    while(j < numPV && idx_crit < numTests) {
      checkUserInterrupt();
      // index of current value in 'pv_list'
      // is also the number of values to be *left out* of the current sum!
      idx_pval = i * size + j;
      // store column values in a vector
      // (re-use variable "pv"; previous values are no longer needed)
      NumericVector temp = NumericVector(mat(_, j));
      // get order
      if(pCDFcounts.isNull()) {
        std::sort(temp.begin(), temp.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(temp, true);
      // this p-value must satisfy condition; else stop computations
      if(temp(ord[0]) <= lambda) {
        // number of remaining needed values
        int rem = numTests - idx_crit;
        // sum for transformations
        double s = 0;
        // compute sum
        int k = 0;
        while(k < numCDF && CDFcounts[ord[k]] < rem) {
          s += CDFcounts[ord[k]] * temp[ord[k]];
          rem -= CDFcounts[ord[k]];
          k++;
        }
        s += rem * temp[ord[k]];
        s /= (1 - lambda) * (idx_crit + 1);
        
        if(idx_crit < numTests && s <= alpha) {
          while(idx_transf < numTests && pv[j] > sorted_pv[idx_transf]) idx_transf++;
          while(idx_transf < numTests && pv[j] == sorted_pv[idx_transf]) {
            pval_transf[idx_transf] = 0;
            rem = numTests - idx_transf;
            int k = 0;
            while(k < numCDF && CDFcounts[ord[k]] < rem) {
              pval_transf[idx_transf] += CDFcounts[ord[k]] * temp[ord[k]];
              rem -= CDFcounts[ord[k]];
              k++;
            }
            pval_transf[idx_transf] += rem * temp[ord[k]];
            pval_transf[idx_transf] /= (1 - lambda) * (idx_transf + 1);
            idx_transf++;
          }
          // go to next p-value in this chunk
          j++;
        } else {
          // current p-value does not satisfiy condition
          // => save previous p-value as critical value
          if(idx_pval) crit[idx_crit++] = pv_list[idx_pval - 1];
          else crit[idx_crit++] = 0;
        }
      } else break;
    }
    // no more p-values will satisfy condition => stop computations
    if(j < numPV && max(mat(_, j)) > lambda) {
      while(idx_crit < numTests) crit[idx_crit++] = pv_list[idx_pval - 1];
      break;
    }
  }
  
  // garbage collection
  delete[] pos;
  delete[] sfuns;
  
  // allow critical values to be 0
  return List::create(Named("crit.consts") = crit, Named("pval.transf") = pval_transf);
}