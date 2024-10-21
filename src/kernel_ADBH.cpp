#include "kernel.h"

NumericVector kernel_ADBH_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool stepUp, const double alpha, const NumericVector &support, const Nullable<NumericVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  
  // intermediate results
  NumericVector pv_list;
  // vector to store transformed p-values
  NumericVector pval_transf;
  // vector to store F_i(tau_m) for SU case only
  std::vector<double> f_denom;
  // lengths of the CDF vectors
  int* lens = new int[numCDF];
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  if(!stepUp) {
    // SD case
    // do not reduce p-value set
    pv_list = sorted_pv;
  } else {
    // SU case
    // get tau_m
    tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, lens, support, numTests, alpha);
    // restrict attention to these values, because tau_k needs to be <= tau_m
    //int i = numTests - 1;
    //while(i > 0 && sorted_pv[i] > tau_m.value) i--;
    pv_list = sorted_pv[Range(0, binary_search(sorted_pv, tau_m.value, numTests))];
    // get pre-computed F_i(tau_m)
    f_denom = tau_m.eval;
  }
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
   // vector to store transformed p-values
  pval_transf = NumericVector(numValues, 0.0);
  // last positions in step function evaluations
  int *pos = new int[numCDF]();
  
  for(int i = 0; i < chunks; i++) {
    checkUserInterrupt();
    // the min(... , numValues) is here for the last chunk
    NumericVector pv = sorted_pv[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(pv))
    for(int j = 0; j < numCDF; j++) {
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], lens[j], pos[j]);
      
      if(stepUp) // SU case, see (10)
        mat(j, _) = mat(j, _) / (1 - f_denom[j]);
      else // SD case, see (11)
        mat(j, _) = mat(j, _) / (1 - mat(j, _));
    }
    // compute transformed p-values
    for(int j = 0; j < numPV; j++) {
      checkUserInterrupt();
      // store column values in a vector
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // get order
      IntegerVector ord;
      if(pCDFcounts.isNull()) {
        std::sort(pv.begin(), pv.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(pv, true);
      // index of current value in pv_list(!)
      // is also the number of values to be *left out* of the current sum!
      int idx_pval = i * size + j;
      // number of remaining needed values
      int rem = numTests - idx_pval;
      // comupte sum
      int k = 0;
      while(k < numCDF && CDFcounts[ord[k]] <= rem) {
        pval_transf[idx_pval] += CDFcounts[ord[k]] * pv[ord[k]];
        rem -= CDFcounts[ord[k]];
        k++;
      }
      if(rem > 0) pval_transf[idx_pval] += rem * pv[ord[k]];
    }
  }
  
  // garbage collection
  delete[] pos;
  delete[] sfuns;
  delete[] lens;
  
  return pval_transf;
}

List kernel_ADBH_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool stepUp, const double alpha, const Nullable<NumericVector> pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  
  // support size
  int numValues = support.length();
  // intermediate results
  NumericVector pv_list;
  // vector to store F_i(tau_m) for SU case only
  std::vector<double> f_denom;
  // lengths of the CDF vectors
  int* lens = new int[numCDF];
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  if(!stepUp) {
    // SD case
    // apply the shortcut drawn from Lemma 3, that is
    // c.1 >= the effective critical value associated to alpha / (numTests + alpha)
    //int i = numValues - 1;
    //while(i > 0 && support[i] >= alpha / (numTests + alpha)) i--;
    int i = binary_search(support, alpha / (numTests + alpha), numValues);
    pv_list = support[Range(i, numValues - 1)];
    // then re-add the observed p-values (needed to compute the adjusted p-values),
    // because we may have removed some of them by the shortcut
    pv_list = rev(sort_combine(sorted_pv, pv_list));
  } else {
    // SU case
    // get tau_m
    tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, lens, support, numTests, alpha);
    // get pre-computed F_i(tau_m)
    f_denom = tau_m.eval;
    // restrict attention p-values <= tau_m, because tau_k needs to be <= tau_m
    pv_list = support[Range(0, tau_m.index)];
    // apply the shortcut drawn from Lemma 4, that is
    // c.1 >= the effective critical value associated to min((1 - tau_m) * alpha/numTests, tau_m)
    //int i = tau_m.index;
    //while(i > 0 && pv_list[i] >= std::min<double>(tau_m.value, (1 - tau_m.value) * alpha / numTests)) i--;
    int i = binary_search(pv_list[Range(0, tau_m.index)], std::min<double>(tau_m.value, (1 - tau_m.value) * alpha / numTests), tau_m.index + 1);    pv_list = rev(pv_list[Range(i, tau_m.index)]);
  }
  
  // number of p-values to be transformed
  numValues = pv_list.length();
  // critical values indices
  NumericVector crit(numTests);
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current critical value
  int idx_crit = numTests - 1;
  // index of current raw p-value to be transformed
  int idx_transf = numTests - 1;
  
  // last positions in step function evaluations
  int *pos = new int[numCDF];
  for(int i = 0; i < numCDF; i++) pos[i] = lens[i] - 1;
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks; i++) {
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(pv))
    for(int j = 0; j < numCDF; j++) {
      checkUserInterrupt();
      for(int k = 0; k < numPV; k++)
        eval_pv_rev(mat(j, k), pv[k], sfuns[j], pos[j]);
      
      if(stepUp) // SU case, see (10)
        mat(j, _) = mat(j, _) / (1 - f_denom[j]);
      else // SD case, see (11)
        mat(j, _) = mat(j, _) / (1 - mat(j, _));
    }
    // compute transformed p-value support (as in pv_list)
    int j = 0;
    while(j < numPV && ((!stepUp && (idx_transf >= 0 || idx_crit >= 0)) || (stepUp && idx_crit >= 0))) {
      checkUserInterrupt();
      // store column values in a vector
      NumericVector temp = NumericVector(mat(_, j));
      // get order
      IntegerVector ord;
      if(pCDFcounts.isNull()) {
        std::sort(temp.begin(), temp.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(temp, true);
      // number of remaining needed values
      int rem = numTests - idx_crit;
      // sum
      double s = 0;
      // comupte sum
      if(idx_crit >= 0) {
        int k = 0;
        while(k < numCDF && CDFcounts[ord[k]] <= rem && s <= alpha * (idx_crit + 1)) {
          s += CDFcounts[ord[k]] * temp[ord[k]];
          rem -= CDFcounts[ord[k]];
          k++;
        }
        if(rem > 0 && s <= alpha * (idx_crit + 1)) s += rem * temp[ord[k]];
      }
      
      // check satisfaction of condition
      if(idx_crit >= 0 && s <= alpha * (idx_crit + 1)) {
        // current p-value satisfies condition
        // => save index of current p-value as critical value
        crit[idx_crit] = pv_list[i * size + j];
        // go to next critical value index to search for
        idx_crit--;
      } else {
        // current p-value does not satisfy condition
        // compute transformed raw p-value for step-down, if there is at least
        // one equal to current support p-value
        if(!stepUp) {
          while(idx_transf >= 0 && pv[j] < sorted_pv[idx_transf]) idx_transf--;
          while(idx_transf >= 0 && pv[j] == sorted_pv[idx_transf]) {
            rem = numTests - idx_transf;
            int k = 0;
            while(k < numCDF && CDFcounts[ord[k]] < rem) {
              pval_transf[idx_transf] += CDFcounts[ord[k]] * temp[ord[k]];
              rem -= CDFcounts[ord[k]];
              k++;
            }
            pval_transf[idx_transf] += rem * temp[ord[k]];
            idx_transf--;
          }
        }
        // go to next p-value in this chunk
        j++;
      }
    }
  }
  // fill remaing critical values (if any)
  while(idx_crit >= 0) {
    crit[idx_crit] = crit[idx_crit + 1];
    idx_crit--;
  }
  
  // garbage collection
  delete[] pos;
  delete[] sfuns;
  delete[] lens;
  
  // output
  if(!stepUp)
    return List::create(Named("crit.consts") = crit, Named("pval.transf") = pval_transf);
  else
    return List::create(Named("crit.consts") = crit);
}
