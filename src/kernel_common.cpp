#include "kernel.h"

tau_m_results DBH_tau_m(const NumericVector* sfuns, const NumericVector& CDFcounts, const int numCDF, const int* lens, const NumericVector& support, const int numTests, const double alpha) {
  // size of the p-value support
  int numValues = support.length();
  
  // index of tau_m
  int index = numValues;
  // tau_m
  double tau_m = 1;
  // step function evaluations
  double eval = 0;
  std::vector<double> tau_m_eval(numCDF);
  // sum
  double sum = numTests * alpha + 1;
  
  // last evaluation positions of step functions
  int* pos = new int[numCDF];
  for(int i = 0; i < numCDF; i++) pos[i] = lens[i] - 1;
  
  // positions for binary search
  int pos_left = 0, pos_right = numValues - 1, pos_mid = numValues - 1;
  // bool variable to indicate end of search
  bool stop = false;
  // traverse step functions increasingly (start with "false")
  bool sf_incr = false;
  // int loop variable
  int i = 0;
  
  // search for tau_m
  while(!stop) {
    // check if user wishes to abort computations
    checkUserInterrupt();
    
    // check if current p-value is smaller than minimum
    if(support[pos_mid] < alpha / (1 + alpha)) {
      pos_left = pos_mid;
      pos_mid = pos_left + (pos_right - pos_mid + 1) / 2;
      continue;
    } else {
      // evaluate pCDFs
      i = 0;
      sum = 0;
      while(i < numCDF && sum <= numTests * alpha) {
        if(sf_incr)
          eval_pv(tau_m_eval[i], support[pos_mid], sfuns[i], lens[i], pos[i]);
        else
          eval_pv_rev(tau_m_eval[i], support[pos_mid], sfuns[i], pos[i]);
        
        sum += CDFcounts[i] * tau_m_eval[i] / (1 - tau_m_eval[i]);
        i++;
      }
      
      if(sum > numTests * alpha) {
        if(pos_mid == 0) {
          // no p-value can satisfy condition
          index = -1;
          tau_m = 0;
          std::fill(tau_m_eval.begin(), tau_m_eval.end(), 0.0);
          stop = true;
        } else if(pos_mid == numValues - 1 && pos_mid - pos_left == 1) {
          index = pos_left;
          tau_m = support[pos_left];
          stop = true;
          for(i = 0; i < numCDF; i++) 
            eval_pv_rev(tau_m_eval[i], support[pos_mid], sfuns[i], pos[i]);
        } else {
          // tau_m MUST be smaller than the current p-value
          pos_right = pos_mid;
          pos_mid = pos_left + (pos_mid - pos_left) / 2;
          sf_incr = false;
        }
      } else {
        // tau_m COULD be larger than the current p-value
        if(pos_mid == numValues - 1 || sum == numTests * alpha || pos_right - pos_mid == 1) {
          // if difference between mid and right position equals 1 or the largest
          //   p-value satisfies the condition, we found tau_m
          index = pos_mid;
          tau_m = support[pos_mid];
          stop = true;
        } else {
          pos_left = pos_mid;
          pos_mid = pos_left + (pos_right - pos_mid + 1) / 2;
          sf_incr = true;
        }
      }
    }
  }
  
  // garbage collection
  delete[] pos;
  
  return {tau_m, index, tau_m_eval};
}

tau_m_results DBH_tau_m2(const NumericVector *sfuns, const NumericVector &CDFcounts, const int numCDF, const NumericVector &support, const int numTests, const double alpha, bool adaptive) {
  // size of the p-value support
  int len = support.length();
  
  // search for tau_m
  int* pos = new int[numCDF];
  for(int j = 0; j < numCDF; j++) pos[j] = sfuns[j].length() - 1;
  double tau_m = 1, sum = numTests * alpha + 1, eval = 0;
  int index = len;
  std::vector<double> tau_m_eval(numCDF);
  while(index > 0 && sum > numTests * alpha){
    index--;
    sum = 0;
    int j = 0;
    while(j < numCDF && sum <= numTests * alpha){
      eval_pv_rev(eval, support[index], sfuns[j], pos[j]);
      //if(adaptive)
      //  sum = std::max<double>(sum, CDFcounts[j] * (eval / (1 - eval)));
      //else
      sum += CDFcounts[j] * (eval / (1 - eval));
      j++;
    }
  }
  // if tau_m is found, compute evaluated F_i
  if(sum <= numTests * alpha){
    tau_m = support[index];
    for(int j = 0; j < numCDF; j++){
      eval_pv_rev(tau_m_eval[j], tau_m, sfuns[j], pos[j]);
    }
  }
  
  // garbage collection
  delete[] pos;
  
  return {tau_m, index, tau_m_eval};
}