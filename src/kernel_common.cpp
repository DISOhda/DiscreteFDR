#include "kernel.h"

tau_m_results DBH_tau_m(const NumericVector* sfuns, const NumericVector& CDFcounts, const int numCDF, const int* lens, const NumericVector& support, const int numTests, const double alpha) {
  // size of the p-value support
  int numValues = support.length();
  
  // index of tau_m
  int idx_tau = numValues;
  // tau_m
  double tau_m = 1;
  // step function evaluations
  std::vector<double> tau_m_eval(numCDF);
  
  // last evaluation positions of step functions
  int* pos = new int[numCDF];
  for(int i = 0; i < numCDF; i++) pos[i] = lens[i] - 1;
  
  // first valid p-value in support
  int idx_first_valid = binary_search(support, alpha / (1 + alpha), numValues) + 1;
  // positions for binary search
  int idx_left = idx_first_valid, idx_right = numValues - 1, idx_mid = numValues - 1;
  // bool variable to indicate end of search
  bool stop = false;
  // traverse step functions increasingly (start with "false")
  bool sf_incr = false;
  // int loop variable
  int i = 0;
  // sum
  double sum = numTests;
  
  // search for tau_m
  while(!stop) {
    // check if user wishes to abort computations
    checkUserInterrupt();
    
    // evaluate pCDFs
    i = 0;
    sum = 0;
    while(i < numCDF) {
      if(sf_incr)
        eval_pv(tau_m_eval[i], support[idx_mid], sfuns[i], lens[i], pos[i]);
      else
        eval_pv_rev(tau_m_eval[i], support[idx_mid], sfuns[i], pos[i]);
    
      sum += CDFcounts[i] * tau_m_eval[i] / (1 - tau_m_eval[i]);
      i++;
    }
    
    if(sum > numTests * alpha) {
      // tau_m MUST be smaller than the current p-value
      if(idx_mid == idx_first_valid) {
        // no p-value can satisfy condition
        idx_tau = idx_mid;
        tau_m = support[idx_tau];
        for(i = 0; i < numCDF; i++) 
          eval_pv_rev(tau_m_eval[i], support[idx_tau], sfuns[i], pos[i]);
        stop = true;
      } else if(idx_mid - idx_left == 1) {
        stop = true;
        // left p-value is the last one to satisfy condition
        idx_tau = idx_left;
        tau_m = support[idx_left];
        for(i = 0; i < numCDF; i++) 
          eval_pv_rev(tau_m_eval[i], support[idx_left], sfuns[i], pos[i]);
      } else {
        // tau_m MUST be smaller than the current p-value
        idx_right = idx_mid;
        idx_mid = idx_left + (idx_mid - idx_left) / 2;
        sf_incr = false;
      }
    } else {
      // tau_m COULD be larger than the current p-value
      if(idx_mid == numValues - 1 || sum == numTests * alpha || idx_right - idx_mid == 1) {
        // if difference between mid and right position equals 1 or the largest
        //   p-value satisfies the condition or sum equals threshold, then we
        //   found tau_m
        stop = true;
        idx_tau = idx_mid;
        tau_m = support[idx_mid];
      } else {
        idx_left = idx_mid;
        idx_mid = idx_left + (idx_right - idx_mid + 1) / 2;
        sf_incr = true;
      }
    }
  }
  
  // garbage collection
  delete[] pos;
  
  return {tau_m, idx_tau, tau_m_eval};
}

/*tau_m_results DBH_tau_m2(const NumericVector *sfuns, const NumericVector &CDFcounts, const int numCDF, const NumericVector &support, const int numTests, const double alpha, bool adaptive) {
  // size of the p-value support
  int len = support.length();
  
  // search for tau_m
  int* pos = new int[numCDF];
  for(int j = 0; j < numCDF; j++) pos[j] = sfuns[j].length() - 1;
  double tau_m = 1, sum = numTests * alpha + 1, eval = 0;
  int idx_tau = len;
  std::vector<double> tau_m_eval(numCDF);
  while(idx_tau > 0 && sum > numTests * alpha){
    idx_tau--;
    sum = 0;
    int j = 0;
    while(j < numCDF && sum <= numTests * alpha){
      eval_pv_rev(eval, support[idx_tau], sfuns[j], pos[j]);
      //if(adaptive)
      //  sum = std::max<double>(sum, CDFcounts[j] * (eval / (1 - eval)));
      //else
      sum += CDFcounts[j] * (eval / (1 - eval));
      j++;
    }
  }
  // if tau_m is found, compute evaluated F_i
  if(sum <= numTests * alpha){
    tau_m = support[idx_tau];
    for(int j = 0; j < numCDF; j++){
      eval_pv_rev(tau_m_eval[j], tau_m, sfuns[j], pos[j]);
    }
  }
  
  // garbage collection
  delete[] pos;
  
  return {tau_m, idx_tau, tau_m_eval};
}*/
