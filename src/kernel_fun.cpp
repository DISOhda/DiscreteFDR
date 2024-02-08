//' @name kernel
//' 
//' @title
//' Kernel functions
//' 
//' @description
//' `r lifecycle::badge('deprecated')`
//'
//' Kernel functions that transform observed p-values or their support according
//' to \[HSU\], \[HSD\], \[AHSU\], \[AHSD\] and \[HBR-\eqn{\lambda}\]. The
//' output is used by [discrete.BH] or [DBR], respectively.
//' `kernel_DBH_crit`, `kernel_ADBH_crit` and `kernel_DBR_crit` additionally
//' compute and return the critical constants. 
//' The end user should not use these functions directly.
//' 
//' **Note**: In future versions, these functions will no longer be exported to
//' the global namespace. Instead, they will be purely internal functions and
//' will have to be called directly via `:::`, e.g.
//' `DiscreteFDR:::kernel_DBH_fast()`.
//' 
//' @details
//' When computing critical constants under step-down, that is, when using
//' `kernel_DBH_crit`, `kernel_ADBH_crit` or `kernel_DBR_crit` with
//' `stepUp = FALSE` (i.e. the step-down case), we still need to get transformed
//' p-values to compute the adjusted p-values.
//' 
//' @seealso
//' [discrete.BH], [fast.Discrete], [DBR]
//' 
//' @templateVar pCDFlist TRUE
//' @templateVar stepf FALSE
//' @templateVar pvalues TRUE
//' @templateVar stepUp TRUE
//' @templateVar alpha FALSE
//' @templateVar alpha2 TRUE
//' @templateVar lambda TRUE
//' @templateVar support TRUE
//' @templateVar sorted_pv TRUE
//' @templateVar raw.pvalues FALSE
//' @templateVar direction FALSE
//' @templateVar ret.crit.consts FALSE
//' @templateVar adaptive FALSE
//' @template param 
//' 
//' @template example
//' @examples
//' 
//' alpha <- 0.05
//' 
//' # Compute the step functions from the supports
//' 
//' # We stay in a step-down context, where pv.numer = pv.denom,
//' # for the sake of simplicity
//' 
//' # If not searching for critical constants, we use only the observed p-values
//' sorted.pvals <- sort(raw.pvalues)
//' y.DBH.fast <- kernel_DBH_fast(pCDFlist, sorted.pvals)
//' y.ADBH.fast <- kernel_ADBH_fast(pCDFlist, sorted.pvals)
//' # transformed values
//' y.DBH.fast
//' y.ADBH.fast
//' 
//' # compute transformed support
//' pv.list <- sort(unique(unlist(pCDFlist)))
//' y.DBH.crit <- kernel_DBH_crit(pCDFlist, pv.list, sorted.pvals)
//' y.ADBH.crit <- kernel_ADBH_crit(pCDFlist, pv.list, sorted.pvals)
//' y.DBR.crit <- kernel_DBR_crit(pCDFlist, pv.list, sorted.pvals)
//' # critical constants
//' y.DBH.crit$crit.consts
//' y.ADBH.crit$crit.consts
//' y.DBR.crit$crit.consts
//' # The following exist only for step-down direction or DBR
//' y.DBH.crit$pval.transf
//' y.ADBH.crit$pval.transf
//' y.DBR.crit$pval.transf
//' 
//' @return
//' For `kernel.DBH.fast`, `kernel.ADBH.fast` and `kernel.DBR.fast`, a vector
//' of transformed p-values is returned. `kernel.DBH.crit`, `kernel.ADBH.crit`
//' `kernel.DBR.crit` return a list with critical constants (`$crit.consts`)
//' and transformed p-values (`$pval.transf`), but if `stepUp = FALSE`, there
//' are critical values only.
//' 
#include <Rcpp.h>
using namespace Rcpp;
 
// fast step function evaluation
// step function is represented by a single numeric vector under the conditions
// a) f(x) = x and b) x is sorted
// this is much faster than passing and evaluating R step function objects
/*IntegerVector stepfun_new(const NumericVector &x, const NumericVector &sfun_x, const int &start = 0){
  // index variables and vector lengths
  int pos = start, size = sfun_x.length(), len = x.length();
  // output vector of the same length as 'x'
  IntegerVector out(len);
  
  // computing results
  for(int i = 0; i < len; i++){
    while(pos < size - 1 && sfun_x[pos] < x[i]) pos++;
    if(sfun_x[pos] == x[i]) out[i] = pos;
    else if(pos) out[i] = pos - 1; else out[i] = -1;
  }
  
  return out;
}*/

struct stepfun_result{
  NumericVector y;
  int last_index;
};

stepfun_result stepfun(const NumericVector &x, const NumericVector &sfun, const int &start = 0){
  // index variables and vector lengths
  int pos = start, size = sfun.length(), len = x.length();
  // output vector of the same length as 'x'
  NumericVector out(len, 1.0);
  
  // computing results
  for(int i = 0; i < len && x[i] < 1; i++){
    while(pos < size - 1 && sfun[pos] < x[i]) pos++;
    if(sfun[pos] == x[i]) out[i] = sfun[pos];
    else if(pos) out[i] = sfun[pos - 1]; else out[i] = 0;
  }
  
  stepfun_result res = {out, pos};
  
  return res;
}

// shortcut function that eliminates all values of a SORTED vector that
// are < limit, except the largest value <= limit
NumericVector short_eff(const NumericVector &x, const double limit){
  // length of the vector
  int len = x.length();
  // identify values <= limit
  NumericVector out = x[x <= limit];
  // eliminate values, but keep their maximum
  out = x[Range(which_max(out), len - 1)];
  
  return out;
}

// function that binds two vectors, sorts it and eliminates duplications 
NumericVector sort_combine(const NumericVector &x, const NumericVector &y){
  // vector lengths
  int lenA = x.length(), lenB = y.length();
  // output vector of the combined lengths of 'x' and 'y'
  NumericVector out(lenA + lenB);
  
  // bind vectors 'x' and 'y'
  for(int i = 0; i < lenA; i++) out[i] = x[i];
  for(int i = 0; i < lenB; i++) out[lenA + i] = y[i];
  // sort and eliminate duplicates
  out = sort_unique(out);
  
  return out;
}

// sort columns of a matrix in descending order
// using an intermediate numeric vector is necessary, because "in-column"
// sorting does not always work as expected, especially for large columns
void colsortdec(NumericMatrix &mat){
  // intermediate vector to store column values (necessary!)
  NumericVector vec;
  for(int i = 0; i < mat.ncol(); i++){
    // store values in a vector
    vec = NumericVector(mat(_, i));
    // sort values in descending order
    std::sort(vec.begin(), vec.end(), std::greater<double>());
    // write sorted values back to column
    mat(_, i) = vec;
  }
}

//'@rdname kernel
//'@export
// [[Rcpp::export]]
NumericVector kernel_DBH_fast(const List &pCDFlist, const NumericVector &pvalues, const bool stepUp = false, const double &alpha = 0.05, const NumericVector &support = 0, const List &pCDFlistIndices = 0){
  // number of p-value CDFs
  //int numCDF = pCDFlist.length();
  // if there are no CDF indices, all CDFs are unique
  //if(pCDFlistIndices.length()){
    
  //}
  // Number of tests
  int numTests = pCDFlist.length();
  // vector to store transformed p-values
  NumericVector pval_transf;
  // p-values to be processed
  NumericVector pv_list;
  // tau_m for SU case
  double tau_m = 1;
  
  if(stepUp){
    // SU case, see (10)
    // apply the shortcut drawn from Lemma 2, that is
    // c.numTests >= the effective critical value associated to alpha / (1 + alpha)
    pv_list = short_eff(support, alpha / (1 + alpha));
    // compute transformed support
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, false, alpha, 0);
    // search the values of the vector <= numTests * alpha
    pv_list = pv_list[pval_transf <= numTests * alpha];
    // get the greatest value (note that Ceiling is sorted in non-decreasing order)
    tau_m = pv_list[pv_list.length() - 1];
    // restrict attention to these values, because tau_k needs to be <= tau_m
    pv_list = pvalues[pvalues <= tau_m];
  }else
    // SD case, see (11)
    pv_list = pvalues;
  
  pval_transf = NumericVector(pv_list.length());
  for(int i = 0; i < numTests; i++){
    checkUserInterrupt();
    NumericVector sfun = as<NumericVector>(pCDFlist[i]);
    NumericVector f_eval = stepfun(pv_list, sfun).y;
    
    if(!stepUp) // SD case, see (11)
      pval_transf += f_eval / (1 - f_eval); // the output is the vector y= \sum_{i=1}^numTests F_i(pvalues)/(1 - F_i(pvalues))
    else{ // SU case, see (10)
      double denom = 1 - stepfun(NumericVector(1, tau_m), sfun).y[0];
      pval_transf += f_eval / denom; // the output is the vector y= \sum_{i=1}^numTests F_i(pvalues)/(1 - F_i(tau.numTests))
    }
    
  }
   
  return pval_transf;
}

//'@rdname kernel
//'@export
// [[Rcpp::export]]
List kernel_DBH_crit(const List &pCDFlist, const NumericVector &pvalues, const NumericVector &sorted_pv, const bool stepUp = false, const double &alpha = 0.05){
  // number of tests
  int numTests = pCDFlist.length();
  // number of p-values to be transformed
  int numValues = pvalues.length();
  // intermediate results
  NumericVector pv_list;
  // vector to store transformed p-values
  NumericVector pval_transf;
  // vector to store critical values indices
  NumericVector crit(numTests);
  
  if(!stepUp){
    // SD case
    // apply the shortcut drawn from Lemma 3, that is tau_1 >= the effective
    // critical value associated to (alpha / numTests) / (1 + alpha / numTests)
    pv_list = short_eff(pvalues, alpha / (numTests + alpha));
    // then re-add the observed p-values (needed to compute the adjusted
    // p-values), because we may have removed some of them by the shortcut
    pv_list = sort_combine(sorted_pv, pv_list);
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, false, alpha, 0);
  }
  else{
    // SU case
    // apply the shortcut drawn from Lemma 2, that is tau_1 >= the effective
    // critical value associated to (alpha / numTests) / (1 + alpha)
    pv_list = short_eff(pvalues, alpha / (numTests + numTests * alpha));
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, true, alpha, pv_list);
  }
  
  // get number of transformed p-values
  numValues = pval_transf.length();
  // get indices of critical values
  int idx_pval = 0;
  for(int i = 1; i <= numTests; i++){
    checkUserInterrupt();
    while(idx_pval < numValues && pval_transf[idx_pval] <= i * alpha) idx_pval++;
    crit[i - 1] = pv_list[idx_pval - 1];
  }
  
  if(!stepUp){
    // SD case
    // store transformed sorted pvalues
    NumericVector transf(numTests);
    // search for sorted p-values in 'pv_list' and get their transformations
    idx_pval = 0;
    for(int i = 0; i < numTests; i++){
      checkUserInterrupt();
      while(pv_list[idx_pval] != sorted_pv[i]) idx_pval++;
      transf[i] = pval_transf[idx_pval];
    }
    // return critical values and transformed sorted p-values
    return List::create(Named("crit.consts") = crit, Named("pval.transf") = transf);
  }
  else
    // SU case
    // return critical values (transformed values are not necessary)
    return List::create(Named("crit.consts") = crit);
}

//'@rdname kernel
//'@export
// [[Rcpp::export]]
NumericVector kernel_ADBH_fast(const List &pCDFlist, const NumericVector &pvalues, const bool stepUp = false, const double alpha = 0.05, const NumericVector &support = 0){
  // number of tests
  int numTests = pCDFlist.length();
  // intermediate results
  NumericVector pv_list;
  // vector to store transformed p-values
  NumericVector pval_transf;
  // vector to store F_i(tau_m) for SU case only
  std::vector<double> f_denom(numTests, 0.0);
  // vector to store F_i
  std::vector<NumericVector> sfun(numTests);
  for(int i = 0; i < numTests; i++)
    sfun[i] = as<NumericVector>(pCDFlist[i]);
  
  if(!stepUp){
    // SD case
    // do not reduce p-value set
    pv_list = pvalues;
  }
  else{
    // SU case
    // apply the shortcut drawn from Lemma 2, that is
    // tau_m >= the effective critical value associated to alpha / (1 + alpha)
    pv_list = short_eff(support, alpha / (numTests + numTests * alpha));
    // get the transformations of the observed p-values inside pv_list
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, false, alpha, 0);
    // search the values of the vector <= numTests * alpha
    pv_list = pv_list[pval_transf <= numTests * alpha];
    // get the greatest value (note that ceiling is sorted in non-decreasing order)
    double tau_m = pv_list[pv_list.length() - 1];
    // restrict attention to these values, because tau_k needs to be <= tau_m
    pv_list = pvalues[pvalues <= tau_m];
    // pre-compute F_i(tau_m) and selection thresholds to avoid unnecessary re-computing
    for(int i = 0; i < numTests; i++){
      f_denom[i] = 1 - stepfun(NumericVector(1, tau_m), sfun[i]).y[0];
    }
  }
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
   // vector to store transformed p-values
  pval_transf = NumericVector(numValues);
  // last positions in step function evaluations
  IntegerVector last(numTests, 0);
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min(... , numValues) is here for the last chunk
    NumericVector pv = pvalues[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(pv))
    for(int j = 0; j < numTests; j++){
      stepfun_result res = stepfun(pv, sfun[j], last[j]);
      last[j] = res.last_index;
      if(stepUp) // SU case, see (10)
        mat(j, _) = res.y / f_denom[j];
      else // SD case, see (11)
        mat(j, _) = res.y / (1 - res.y);
    }
    // sort columns in descending order
    colsortdec(mat);
    // compute transformed p-values
    for(int j = 0; j < len; j++){
      checkUserInterrupt();
      // index of current value in pv_list (!)
      int idx_pval = i * size + j;
      // evaluated F_j in descending order
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // compute sum
      pval_transf[idx_pval] = sum(pv[Range(0, numTests - idx_pval - 1)]);
    }
  }
  
  return pval_transf;
}

//'@rdname kernel
//'@export
// [[Rcpp::export]]
List kernel_ADBH_crit(const List &pCDFlist, const NumericVector &pvalues, const NumericVector &sorted_pv, const bool stepUp = false, const double alpha = 0.05){
  // index variables and lengths
  int numTests = pCDFlist.length();
  // intermediate results
  NumericVector pv_list;
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // vector to store F_i(tau_m) for SU case only
  std::vector<double> f_denom(numTests, 0.0);
  // vector to store F_i
  std::vector<NumericVector> sfun(numTests);
  for(int i = 0; i < numTests; i++)
    sfun[i] = as<NumericVector>(pCDFlist[i]);
  
  if(!stepUp){
    // SD case
    // apply the shortcut drawn from Lemma 3, that is
    // c.1 >= the effective critical value associated to (alpha / numTests) / (1 + alpha / numTests)
    double tau_1 = alpha / (numTests + alpha);
    pv_list = short_eff(pvalues, tau_1);
    // then re-add the observed p-values (needed to compute the adjusted p-values),
    // because we may have removed some of them by the shortcut
    pv_list = rev(sort_combine(sorted_pv, pv_list));
  }
  else{
    // SU case
    // apply the shortcut drawn from Lemma 2, that is
    // c.1 >= the effective critical value associated to (alpha / numTests) / (1 + alpha)
    pv_list = short_eff(pvalues, alpha / (numTests + numTests * alpha));
    // compute transformed support
    pval_transf = kernel_DBH_fast(pCDFlist, pv_list, false, alpha, 0);
    // get tau_m: search largest value of the vector which is <= numTests * alpha
    pv_list = pv_list[pval_transf <= numTests * alpha];
    int i = pval_transf.length() - 1;
    double tau_m = pv_list[i];
    while(i > 0 && pval_transf[i] > numTests * alpha) i--;
    if(pval_transf[i] <= numTests * alpha) tau_m = pv_list[i];
    // restrict attention to these values, because tau_k needs to be <= tau_m
    i = pvalues.length() - 1;
    while(i > 0 && pvalues[i] > tau_m) i--;
    pv_list = pvalues[Range(0, i)];
    // apply the shortcut drawn from Lemma 4, that is
    // c.1 >= the effective critical value associated to min((1 - tau_m) * alpha/numTests, tau_m)
    pv_list = rev(short_eff(pv_list, std::min<double>(tau_m, (1 - tau_m) * alpha / numTests)));
    // pre-compute F_i(tau_m) to avoid unnecessary re-computing
    for(i = 0; i < numTests; i++){
      f_denom[i] = 1 - stepfun(NumericVector(1, tau_m), sfun[i]).y[0];
    }
  }
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // critical values indices
  IntegerVector crit(numTests);
  // initialize critical values indices
  //crit.fill(pv_list.length() - 1);
  // index of current critical value
  int idx_crit = numTests - 1;
  // index of current raw p-value to be transformed
  int idx_transf = numTests - 1;
  // last positions in step function evaluations
  IntegerVector last(numTests, 0);
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    for(int j = 0; j < numTests; j++){
      stepfun_result res = stepfun(rev(pv), sfun[j], last[j]);
      last[j] = res.last_index;
      NumericVector f_eval = rev(res.y);
      if(stepUp) // SU case, see (12)
        mat(j, _) = f_eval / f_denom[j];   // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(tau_m))
      else // SD case, see (13)
        mat(j, _) = f_eval / (1 - f_eval); // compute columns \sum_{j=1}^numTests F_j(pv)/(1 - F_j(pv))
    }
    // sort columns in descending order
    colsortdec(mat);
    // compute transformed p-value support (as in pv_list)
    int j = 0;
    while(j < len && ((!stepUp && (idx_transf >= 0 || idx_crit >= 0)) || (stepUp && idx_crit >= 0))){
      checkUserInterrupt();
      // evaluated F_j in descending order
      NumericVector temp = NumericVector(mat(_, j));
      // sum
      double s;
      // compute sum
      if(idx_crit >= 0){
        s = sum(temp[Range(0, numTests - idx_crit - 1)]);
      }
      // check satisfaction of condition
      if(idx_crit >= 0 && s <= alpha * (idx_crit + 1)){
        // current p-value satisfies condition
        // => save index of current p-value as critical value
        crit[idx_crit] = i * size + j;
        // go to next critical value index to search for
        --idx_crit;
      }else{
        // current p-value does not satisfy condition
        // go to next p-value in this chunk
        // compute transformed raw p-value for step-down, if there is at least
        // one equal to current support p-value
        if(!stepUp){
          while(idx_transf >= 0 && pv[j] < sorted_pv[idx_transf]) idx_transf--;
          while(idx_transf >= 0 && pv[j] == sorted_pv[idx_transf]){
            pval_transf[idx_transf] = sum(temp[Range(0, numTests - idx_transf - 1)]);
            idx_transf--;
          }
        }
        j++;
      }
    }
  }
  
  for(;idx_crit >= 0; idx_crit--) crit[idx_crit] = crit[idx_crit + 1];
  
  // output
  if(!stepUp)
    return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf);
  else
    return List::create(Named("crit.consts") = pv_list[crit]);
}

//'@rdname kernel
//'@export
// [[Rcpp::export]]
NumericVector kernel_DBR_fast(const List &pCDFlist, const NumericVector &pvalues, const double &lambda = 0.05){
  // index variables and lengths
  int numTests = pCDFlist.length();
  // number of p-values to be transformed
  int numValues = pvalues.length();
  // intermediate results
  NumericVector pv_list;
  // vector to store transformed p-values
  NumericVector pval_transf(numValues, 1.0);
  // vector to store F_i
  std::vector<NumericVector> sfun(numTests);
  for(int i = 0; i < numTests; i++)
    sfun[i] = as<NumericVector>(pCDFlist[i]);
  
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // last positions in step function evaluations
  IntegerVector last(numTests, 0);
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pvalues[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numTests; j++){
      stepfun_result res = stepfun(pv, sfun[j], last[j]);
      last[j] = res.last_index;
      mat(j, _) = res.y;
    }
    
    // sort columns in descending order
    colsortdec(mat);
    
    int j = 0;
    while(j < len && mat(0, j) <= lambda){
      checkUserInterrupt();
      // index of current value in 'pv_list'
      int idx_pval = i * size + j;
      // evaluated F_j in descending order
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // compute sum
      pval_transf[idx_pval] = sum(pv[Range(0, numTests - idx_pval - 1)]) / ((1 - lambda) * (idx_pval + 1));
      // go to next p-value in this chunk
      j++;
    }
    // no more p-values will satisfy condition => stop computations
    if(j < len && mat(0, j) > lambda) break;
  }
  
  return pval_transf;
}

//'@rdname kernel
//'@export
// [[Rcpp::export]]
List kernel_DBR_crit(const List &pCDFlist, const NumericVector &pvalues, const NumericVector &sorted_pv, const double &lambda = 0.05, const double &alpha = 0.05){
  // number of tests
  int numTests = pCDFlist.length(), k_star = 0;
  // critical values indices
  IntegerVector crit(numTests);
  // transformed p-values
  NumericVector pval_transf(numTests, 1.0);
  // vector to store F_i
  std::vector<NumericVector> sfun(numTests);
  for(int i = 0; i < numTests; i++)
    sfun[i] = as<NumericVector>(pCDFlist[i]);
  
  // apply the shortcut drawn from Corollary 3, that is
  // c.1 >= the effective critical value associated to min((1 - lambda) * alpha/numTests , lambda)
  NumericVector pv_list = short_eff(pvalues, std::min<double>(lambda, (1 - lambda) * alpha / numTests));
  pv_list = sort_combine(pv_list, sorted_pv);
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current value in 'pv_list'
  int idx_pval = 0;
  // index of current critical value
  int idx_crit = 0;
  // index of current raw p-value to be transformed
  int idx_transf = 0;
  
  // last positions in step function evaluations
  IntegerVector last(numTests, 0);
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numTests; j++){
      stepfun_result res = stepfun(pv, sfun[j], last[j]);
      last[j] = res.last_index;
      mat(j, _) = res.y;
    }
    // sort columns in descending order
    colsortdec(mat);
    
    int j = 0;
    // sum for transformations
    double s = 0;
    while(j < len && mat(0, j) <= lambda && idx_crit < numTests){
      checkUserInterrupt();
      // index of current value in 'pv_list'
      idx_pval = i * size + j;
      // evaluated F_j in descending order
      // (re-use variable "pv"; previous values are no longer needed)
      NumericVector temp = NumericVector(mat(_, j));
      // compute sum
      s = sum(temp[Range(0, numTests - idx_crit - 1)]) / (1 - lambda);
      if(idx_crit < numTests && s <= alpha * (idx_crit + 1)){
        while(idx_transf < numTests && pv[j] > sorted_pv[idx_transf]) idx_transf++;
        while(idx_transf < numTests && pv[j] == sorted_pv[idx_transf]){
          pval_transf[idx_transf] = sum(temp[Range(0, numTests - idx_transf - 1)]) / ((1 - lambda) * (idx_transf + 1));
          k_star = ++idx_transf;
        }
        // go to next p-value in this chunk
        j++;
      }else{
        // current p-value does not satisfiy condition
        // => save index of previous p-value as critical value
        crit[idx_crit++] = idx_pval - 1;
      }
    }
    if(j < len && mat(0, j) > lambda){
      while(idx_crit < numTests) crit[idx_crit++] = idx_pval;
      break;
    }
  }
  
  // allow critical values equal to 0
  pv_list.push_front(0);
  return List::create(Named("crit.consts") = pv_list[crit + 1], Named("pval.transf") = pval_transf, Named("m.lambda") = k_star);
}