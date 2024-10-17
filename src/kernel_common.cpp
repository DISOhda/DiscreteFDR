#include "kernel.h"

tau_m_results DBH_tau_m(const NumericVector *sfuns, const NumericVector &CDFcounts, const int numCDF, const NumericVector &support, const int numTests, const double alpha, bool adaptive) {
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