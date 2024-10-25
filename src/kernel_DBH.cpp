#include "kernel.h"

NumericVector kernel_DBH_fast(const List &pCDFlist, const NumericVector &pvalues, const bool stepUp, const Nullable<NumericVector> tau_max, const double alpha, const NumericVector &support, const Nullable<NumericVector> &pCDFcounts) {
  // Number of p-values
  int numValues = pvalues.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution and number of tests
  int numTests;
  NumericVector CDFcounts;
  if(numCDF == numValues || pCDFcounts.isNull()) {
    CDFcounts = NumericVector(numCDF, 1.0);
    numTests = numCDF;
  } else {
    CDFcounts = pCDFcounts;
    numTests = sum(CDFcounts);
  }
  // vector to store transformed p-values
  NumericVector pval_transf;
  // p-values to be processed
  NumericVector pv_list;
  // lengths of the CDF vectors
  int* lens = new int[numCDF];
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  // tau_m for SU case
  tau_m_results tau_m;
  
  if(stepUp) {
    // SU case, see (10)
    // get tau_m
    if(tau_max.isNull()) {
      tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, lens, support, numTests, alpha);
    } else {
      tau_m.value = NumericVector(tau_max)[0];
      tau_m.eval = std::vector<double>(numCDF);
      for(int i = 0; i < numCDF; i++) {
        int pos = 0;
        eval_pv(tau_m.eval[i], tau_m.value, sfuns[i], lens[i], pos);
      }
    }
    
    // restrict attention to values <= tau_m, because tau_k needs to be <= tau_m
    pv_list = pvalues[Range(0, binary_search(pvalues, tau_m.value, numValues))];
    // update number of relevant p-values
    numValues = pv_list.length();
  } else
    // SD case, see (11)
    pv_list = pvalues;
  
  pval_transf = NumericVector(numValues);
  NumericVector f_eval(numValues);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    int pos = 0;
    int len = sfuns[i].length();
    for(int j = 0; j < numValues; j++) {
      eval_pv(f_eval[j], pv_list[j], sfuns[i], len, pos);
    }
    
    if(!stepUp) // SD case, see (11)
      pval_transf += CDFcounts[i] * (f_eval / (1 - f_eval)); // the output is the vector y= \sum_{i=1}^numTests F_i(pvalues)/(1 - F_i(pvalues))
    else { // SU case, see (10)
      pval_transf += CDFcounts[i] * (f_eval / (1 - tau_m.eval[i])); // the output is the vector y= \sum_{i=1}^numTests F_i(pvalues)/(1 - F_i(tau.numTests))
    }
  }
  
  // garbage collection
  delete[] sfuns;
  delete[] lens;
  
  return pval_transf;
}

List kernel_DBH_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool stepUp, const double alpha, const Nullable<NumericVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  
  // number of p-values to be transformed
  int numValues = support.length();
  // intermediate results
  NumericVector pv_list;
  // vector to store transformed p-values
  NumericVector pval_transf;
  // vector to store critical values indices
  NumericVector crit(numTests);
  // lengths of the CDF vectors
  int* lens = new int[numCDF];
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // get tau_m
  tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, lens, support, numTests, alpha);
  // set last critical value
  crit[numTests - 1] = tau_m.value;
  
  if(!stepUp) {
    // SD case
    // apply the shortcut drawn from Lemma 2, that is tau_1 >= the effective
    // critical value associated to alpha / (m + alpha)
    //int i = numValues - 1;
    //while(i > 0 && support[i] >= alpha / (numTests + alpha)) i--;
    int i = binary_search(support, alpha / (numTests + alpha), tau_m.index + 1);
    // then re-add the observed p-values (needed to compute the adjusted
    // p-values), because we may have removed some of them by the shortcut
    pv_list = sort_combine(sorted_pv, support[Range(i, tau_m.index)]);
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, false, R_NilValue, alpha, NumericVector(), pCDFcounts);
  } else { 
    // SU case
    // apply the shortcut drawn from Lemma 2, that is tau_1 >= the effective
    // critical value associated to (alpha / numTests) / (1 + alpha)
    //int i = numValues - 1;
    //while(i > 0 && support[i] > alpha / (numTests + numTests * alpha)) i--;
    int i = binary_search(support, alpha / (numTests + numTests * alpha), numValues);
    // then re-add the observed p-values (needed to compute the adjusted
    // p-values), because we may have removed some of them by the shortcut
    pv_list = sort_combine(sorted_pv, support[Range(i, tau_m.index)]);
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, true, NumericVector(1, tau_m.value), alpha, pv_list, pCDFcounts);
    //pval_transf = kernel_DBH_fast(pCDFlist, pv_list, true, R_NilValue, alpha, pv_list, pCDFcounts);
  }
  
  // get number of transformed p-values
  numValues = pval_transf.length();
  // get indices of critical values
  int idx_pval = 0;
  for(int i = 0; i < numTests - 1; i++) {
    checkUserInterrupt();
    while(idx_pval < numValues && pval_transf[idx_pval] <= (i + 1) * alpha) idx_pval++;
    crit[i] = pv_list[idx_pval - 1];
  }
  
  // garbage collection
  delete[] sfuns;
  delete[] lens;
  
  if(!stepUp) {
    // SD case
    // store transformed sorted pvalues
    NumericVector transf(numTests);
    // search for sorted p-values in 'pv_list' and get their transformations
    idx_pval = 0;
    for(int i = 0; i < numTests; i++) {
      checkUserInterrupt();
      while(idx_pval < numValues && pv_list[idx_pval] != sorted_pv[i]) idx_pval++;
      transf[i] = pval_transf[idx_pval];
    }
    // return critical values and transformed sorted p-values
    return List::create(Named("crit.consts") = crit, Named("pval.transf") = transf);
  } else {
    // SU case
    // return critical values (transformed values are not necessary)
    return List::create(Named("crit.consts") = crit);
  }
}
