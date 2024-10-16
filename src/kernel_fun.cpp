#include "helper.h"

//' @name kernel
//' @rdname kernel
//' 
//' @title
//' Kernel Functions
//' 
//' @description
//'
//' Kernel functions that transform observed p-values or their support according
//' to \[HSU\], \[HSD\], \[AHSU\], \[AHSD\] and \[HBR-\eqn{\lambda}\]. The
//' output is used by [discrete.BH] or [DBR], respectively.
//' `kernel_DBH_crit`, `kernel_ADBH_crit` and `kernel_DBR_crit` additionally
//' compute and return the critical constants. 
//' The end user should not use these functions directly.
//' 
//' **Note**: As of version 2.0, these functions are purely internal functions!
//' As a consequence, they have to be called directly via `:::`, e.g. 
//' `DiscreteFDR:::kernel_DBH_fast()`. But users should **not** rely on them, as
//' parameters (including their names, order, etc.) may be changed without
//' notice!
//' 
//' @templateVar pCDFlist TRUE
//' @template param
//' 
//' @param pvalues      numeric vector, sorted in increasing order, that either
//'                     must contain the entirety of all observable values of
//'                     the p-value supports (when computing critical constants)
//'                     or only the sorted raw p-values.
//' @param stepUp       boolean specifying whether to conduct the step-up
//'                     (`TRUE`) or step-down (`FALSE`; the default)
//'                     procedure.
//' @param alpha        single real number strictly between 0 and 1 indicating
//'                     the target FDR level; for `*.fast` kernels, it is only
//'                     needed, if `stepUp = TRUE`.
//' @param support      numeric vector, sorted in increasing order, that
//'                     contains the entirety of all observable values of the
//'                     p-value supports; for `*.fast` kernels, it is ignored if
//'                     `stepUp = FALSE`.
//' @param pCDFcounts   integer vector of counts that indicates to how many
//'                     p-values each **unique** p-value distributions belongs.
//' @param sorted_pv    numeric vector containing the raw p-values, sorted in
//'                     increasing order.
//' @param lambda       real number strictly between 0 and 1 specifying the DBR
//'                     tuning parameter.
//' 
//' @details
//' When computing critical constants under step-down, that is, when using
//' `kernel_DBH_crit`, `kernel_ADBH_crit` or `kernel_DBR_crit` with
//' `stepUp = FALSE` (i.e. the step-down case), we still need to get transformed
//' p-values to compute the adjusted p-values.
//' 
//' @return
//' For `kernel.DBH.fast`, `kernel.ADBH.fast` and `kernel.DBR.fast`, a vector
//' of transformed p-values is returned. `kernel.DBH.crit`, `kernel.ADBH.crit`
//' `kernel.DBR.crit` return a list with critical constants (`$crit.consts`)
//' and transformed p-values (`$pval.transf`), but if `stepUp = FALSE`, there
//' are critical values only.
//' 
//' @seealso
//' [`discrete.BH`], [`fast.Discrete`], [`DBR`]
//' 
//' @examples \dontrun{
//' X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
//' X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
//' N1 <- rep(148, 9)
//' N2 <- rep(132, 9)
//' Y1 <- N1 - X1
//' Y2 <- N2 - X2
//' df <- data.frame(X1, Y1, X2, Y2)
//' df
//' 
//' # Compute p-values and their supports of Fisher's exact test
//' test.result <- generate.pvalues(df, "fisher")
//' raw.pvalues <- test.result$get_pvalues()
//' pCDFlist <- test.result$get_pvalue_supports()
//' 
//' alpha <- 0.05
//' 
//' # Compute the step functions from the supports
//' 
//' # If not searching for critical constants, we use only the observed p-values
//' sorted.pvals   <- sort(raw.pvalues)
//' y.DBH.sd.fast  <- DiscreteFDR:::kernel_DBH_fast(pCDFlist, sorted.pvals)
//' y.ADBH.sd.fast <- DiscreteFDR:::kernel_ADBH_fast(pCDFlist, sorted.pvals)
//' y.DBR.fast     <- DiscreteFDR:::kernel_DBR_fast(pCDFlist, sorted.pvals)
//' # transformed values
//' y.DBH.sd.fast
//' y.ADBH.sd.fast
//' y.DBR.fast
//' 
//' # compute transformed support
//' pv.list        <- sort(unique(unlist(pCDFlist)))
//' y.DBH.sd.crit  <- DiscreteFDR:::kernel_DBH_crit(pCDFlist, pv.list, sorted.pvals)
//' y.ADBH.sd.crit <- DiscreteFDR:::kernel_ADBH_crit(pCDFlist, pv.list, sorted.pvals)
//' y.DBR.crit     <- DiscreteFDR:::kernel_DBR_crit(pCDFlist, pv.list, sorted.pvals)
//' # critical constants
//' y.DBH.sd.crit$crit.consts
//' y.ADBH.sd.crit$crit.consts
//' y.DBR.crit$crit.consts
//' # The following exist only for step-down direction or DBR
//' y.DBH.sd.crit$pval.transf
//' y.ADBH.sd.crit$pval.transf
//' y.DBR.crit$pval.transf
//' }
//' 
struct tau_m_results{
  double value;
  int index;
  std::vector<double> eval;
};

tau_m_results DBH_tau_m(const NumericVector *sfuns, const NumericVector &CDFcounts, const int numCDF, const NumericVector &support, const int numTests, const double alpha, bool adaptive = false){
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

///' @export
//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DBH_fast(const List &pCDFlist, const NumericVector &pvalues, const bool stepUp = false, const double &alpha = 0.05, const NumericVector &support = NumericVector(), const Nullable<NumericVector> &pCDFcounts = R_NilValue){
  // Number of p-values
  int numValues = pvalues.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution and number of tests
  int numTests;
  NumericVector CDFcounts;
  if(numCDF == numValues || pCDFcounts.isNull()){
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
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  // tau_m for SU case
  tau_m_results tau_m;
  
  if(stepUp){
    // SU case, see (10)
    // get tau_m
    tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, support, numTests, alpha);
    // restrict attention to values <= tau_m, because tau_k needs to be <= tau_m
    int j = numValues - 1;
    while(j > 0 && pvalues[j] > tau_m.value) j--;
    pv_list = pvalues[Range(0, j)];
    numValues = j + 1;
  }else
    // SD case, see (11)
    pv_list = pvalues;
  
  pval_transf = NumericVector(numValues);
  NumericVector f_eval(numValues);
  for(int i = 0; i < numCDF; i++){
    checkUserInterrupt();
    int pos = 0;
    int len = sfuns[i].length();
    for(int j = 0; j < numValues; j++){
      eval_pv(f_eval[j], pv_list[j], sfuns[i], len, pos);
    }
    
    if(!stepUp) // SD case, see (11)
      pval_transf += CDFcounts[i] * (f_eval / (1 - f_eval)); // the output is the vector y= \sum_{i=1}^numTests F_i(pvalues)/(1 - F_i(pvalues))
    else{ // SU case, see (10)
      pval_transf += CDFcounts[i] * (f_eval / (1 - tau_m.eval[i])); // the output is the vector y= \sum_{i=1}^numTests F_i(pvalues)/(1 - F_i(tau.numTests))
    }
  }
  
  // garbage collection
  delete[] sfuns;
  
  return pval_transf;
}

///' @export
//' @rdname kernel
// [[Rcpp::export]]
List kernel_DBH_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool stepUp = false, const double &alpha = 0.05, const Nullable<NumericVector> &pCDFcounts = R_NilValue){//const List pCDFindices = R_NilValue
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
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  if(!stepUp){
    // SD case
    // get tau_m
    tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, support, numTests, alpha);
    crit[numTests - 1] = tau_m.value;
    // apply the shortcut drawn from Lemma 2, that is
    // c.1 >= the effective critical value associated to alpha / (m + alpha)
    int i = numValues - 1;
    while(i > 0 && support[i] >= alpha / (numTests + alpha)) i--;
    // then re-add the observed p-values (needed to compute the adjusted
    // p-values), because we may have removed some of them by the shortcut
    pv_list = sort_combine(sorted_pv, support[Range(i, tau_m.index)]);
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, false, alpha, NumericVector(), pCDFcounts);
  }
  else{
    // SU case
    // get tau_m
    tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, support, numTests, alpha);
    crit[numTests - 1] = tau_m.value;
    // apply the shortcut drawn from Lemma 2, that is tau_1 >= the effective
    // critical value associated to (alpha / numTests) / (1 + alpha)
    int i = numValues - 1;
    while(i > 0 && support[i] > alpha / (numTests + numTests * alpha)) i--;
    // then re-add the observed p-values (needed to compute the adjusted
    // p-values), because we may have removed some of them by the shortcut
    pv_list = sort_combine(sorted_pv, support[Range(i, tau_m.index)]);
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, true, alpha, pv_list, pCDFcounts);
  }
  
  // get number of transformed p-values
  numValues = pval_transf.length();
  // get indices of critical values
  int idx_pval = 0;
  for(int i = 0; i < numTests - 1; i++){
    checkUserInterrupt();
    while(idx_pval < numValues && pval_transf[idx_pval] <= (i + 1) * alpha) idx_pval++;
    crit[i] = pv_list[idx_pval - 1];
  }
  
  // garbage collection
  delete[] sfuns;
  
  if(!stepUp){
    // SD case
    // store transformed sorted pvalues
    NumericVector transf(numTests);
    // search for sorted p-values in 'pv_list' and get their transformations
    idx_pval = 0;
    for(int i = 0; i < numTests; i++){
      checkUserInterrupt();
      while(idx_pval < numValues && pv_list[idx_pval] != sorted_pv[i]) idx_pval++;
      transf[i] = pval_transf[idx_pval];
    }
    // return critical values and transformed sorted p-values
    return List::create(Named("crit.consts") = crit, Named("pval.transf") = transf);
  }
  else{
    // SU case
    // return critical values (transformed values are not necessary)
    return List::create(Named("crit.consts") = crit);
  }
}

///' @export
//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_ADBH_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool stepUp = false, const double alpha = 0.05, const NumericVector &support = NumericVector(), const Nullable<NumericVector> &pCDFcounts = R_NilValue){//const List pCDFindices = R_NilValue
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
  // vector to store F_i
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  if(!stepUp){
    // SD case
    // do not reduce p-value set
    pv_list = sorted_pv;
  }
  else{
    // SU case
    // get tau_m
    tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, support, numTests, alpha);
    // restrict attention to these values, because tau_k needs to be <= tau_m
    int i = numTests - 1;
    while(i > 0 && sorted_pv[i] > tau_m.value) i--;
    pv_list = sorted_pv[Range(0, i)];
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
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min(... , numValues) is here for the last chunk
    NumericVector pv = sorted_pv[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(pv))
    for(int j = 0; j < numCDF; j++){
      int len = sfuns[j].length();
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], len, pos[j]);
      
      if(stepUp) // SU case, see (10)
        mat(j, _) = mat(j, _) / (1 - f_denom[j]);
      else // SD case, see (11)
        mat(j, _) = mat(j, _) / (1 - mat(j, _));
    }
    // compute transformed p-values
    for(int j = 0; j < numPV; j++){
      checkUserInterrupt();
      // store column values in a vector
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // get order
      IntegerVector ord;
      if(pCDFcounts.isNull()){
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
      while(k < numCDF && CDFcounts[ord[k]] <= rem){
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
  
  return pval_transf;
}

/// '@export
//' @rdname kernel
// [[Rcpp::export]]
List kernel_ADBH_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool stepUp = false, const double alpha = 0.05, const Nullable<NumericVector> pCDFcounts = R_NilValue){
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
  // vector to store F_i
  NumericVector *sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  if(!stepUp){
    // SD case
    // apply the shortcut drawn from Lemma 3, that is
    // c.1 >= the effective critical value associated to alpha / (numTests + alpha)
    int i = numValues - 1;
    while(i > 0 && support[i] >= alpha / (numTests + alpha)) i--;
    pv_list = support[Range(i, numValues - 1)];
    // then re-add the observed p-values (needed to compute the adjusted p-values),
    // because we may have removed some of them by the shortcut
    pv_list = rev(sort_combine(sorted_pv, pv_list));
  }
  else{
    // SU case
    // get tau_m
    tau_m_results tau_m = DBH_tau_m(sfuns, CDFcounts, numCDF, support, numTests, alpha);
    // get pre-computed F_i(tau_m)
    f_denom = tau_m.eval;
    // restrict attention p-values <= tau_m, because tau_k needs to be <= tau_m
    pv_list = support[Range(0, tau_m.index)];
    // apply the shortcut drawn from Lemma 4, that is
    // c.1 >= the effective critical value associated to min((1 - tau_m) * alpha/numTests, tau_m)
    int i = tau_m.index;
    while(i > 0 && pv_list[i] >= std::min<double>(tau_m.value, (1 - tau_m.value) * alpha / numTests)) i--;
    pv_list = rev(pv_list[Range(i, tau_m.index)]);
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
  for(int i = 0; i < numCDF; i++) pos[i] = sfuns[i].length() - 1;
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks; i++){
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(pv))
    for(int j = 0; j < numCDF; j++){
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
    while(j < numPV && ((!stepUp && (idx_transf >= 0 || idx_crit >= 0)) || (stepUp && idx_crit >= 0))){
      checkUserInterrupt();
      // store column values in a vector
      NumericVector temp = NumericVector(mat(_, j));
      // get order
      IntegerVector ord;
      if(pCDFcounts.isNull()){
        std::sort(temp.begin(), temp.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(temp, true);
      // number of remaining needed values
      int rem = numTests - idx_crit;
      // sum
      double s = 0;
      // comupte sum
      if(idx_crit >= 0){
        int k = 0;
        while(k < numCDF && CDFcounts[ord[k]] <= rem && s <= alpha * (idx_crit + 1)){
          s += CDFcounts[ord[k]] * temp[ord[k]];
          rem -= CDFcounts[ord[k]];
          k++;
        }
        if(rem > 0 && s <= alpha * (idx_crit + 1)) s += rem * temp[ord[k]];
      }
      
      // check satisfaction of condition
      if(idx_crit >= 0 && s <= alpha * (idx_crit + 1)){
        // current p-value satisfies condition
        // => save index of current p-value as critical value
        crit[idx_crit] = pv_list[i * size + j];
        // go to next critical value index to search for
        idx_crit--;
      }else{
        // current p-value does not satisfy condition
        // compute transformed raw p-value for step-down, if there is at least
        // one equal to current support p-value
        if(!stepUp){
          while(idx_transf >= 0 && pv[j] < sorted_pv[idx_transf]) idx_transf--;
          while(idx_transf >= 0 && pv[j] == sorted_pv[idx_transf]){
            rem = numTests - idx_transf;
            int k = 0;
            while(k < numCDF && CDFcounts[ord[k]] < rem){
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
  for( ; idx_crit >= 0; idx_crit--) crit[idx_crit] = crit[idx_crit + 1];
  
  // garbage collection
  delete[] pos;
  delete[] sfuns;
  
  // output
  if(!stepUp)
    return List::create(Named("crit.consts") = crit, Named("pval.transf") = pval_transf);
  else
    return List::create(Named("crit.consts") = crit);
}

///' @export
//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DBR_fast(const List &pCDFlist, const NumericVector &sorted_pv, const double &lambda = 0.05, const Nullable<NumericVector> pCDFcounts = R_NilValue){
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
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numTests) is here for the last chunk
    NumericVector pv = sorted_pv[Range(i * size, std::min<int>((i + 1) * size, numTests) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numCDF; j++){
      int len = sfuns[j].length();
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], len, pos[j]);
    }
    
    int j;
    for(j = 0; j < numPV; j++){
      checkUserInterrupt();
      // index of current value in 'sorted_pv'
      // is also the number of values to be *left out* of the current sum!
      int idx_pval = i * size + j;
      // store column values in a vector
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // get order
      if(pCDFcounts.isNull()){
        std::sort(pv.begin(), pv.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(pv, true);
      // this p-value must satisfy condition; else stop computations
      if(pv(ord[0]) <= lambda){
        // number of remaining needed values
        int rem = numTests - idx_pval;
        // compute sum
        int k = 0;
        pval_transf[idx_pval] = 0;
        while(k < numCDF && CDFcounts[ord[k]] < rem){
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

///' @export
//' @rdname kernel
// [[Rcpp::export]]
List kernel_DBR_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const double &lambda = 0.05, const double &alpha = 0.05, const Nullable<NumericVector> pCDFcounts = R_NilValue){
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
  int i = numValues - 1;
  while(i > 0 && support[i] >= std::min<double>(lambda, (1 - lambda) * alpha / numTests)) i--;
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
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numCDF; j++){
      int len = sfuns[j].length();
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], len, pos[j]);
    }
    
    // loop variable
    int j = 0;
    while(j < numPV && idx_crit < numTests){
      checkUserInterrupt();
      // index of current value in 'pv_list'
      // is also the number of values to be *left out* of the current sum!
      idx_pval = i * size + j;
      // store column values in a vector
      // (re-use variable "pv"; previous values are no longer needed)
      NumericVector temp = NumericVector(mat(_, j));
      // get order
      if(pCDFcounts.isNull()){
        std::sort(temp.begin(), temp.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(temp, true);
      // this p-value must satisfy condition; else stop computations
      if(temp(ord[0]) <= lambda){
        // number of remaining needed values
        int rem = numTests - idx_crit;
        // sum for transformations
        double s = 0;
        // compute sum
        int k = 0;
        while(k < numCDF && CDFcounts[ord[k]] < rem){
          s += CDFcounts[ord[k]] * temp[ord[k]];
          rem -= CDFcounts[ord[k]];
          k++;
        }
        s += rem * temp[ord[k]];
        s /= (1 - lambda) * (idx_crit + 1);
        
        if(idx_crit < numTests && s <= alpha){
          while(idx_transf < numTests && pv[j] > sorted_pv[idx_transf]) idx_transf++;
          while(idx_transf < numTests && pv[j] == sorted_pv[idx_transf]){
            pval_transf[idx_transf] = 0;
            rem = numTests - idx_transf;
            int k = 0;
            while(k < numCDF && CDFcounts[ord[k]] < rem){
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
    if(j < numPV && max(mat(_, j)) > lambda){
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