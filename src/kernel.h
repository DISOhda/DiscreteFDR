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
//' @param tau_max      single real number strictly between 0 and 1 indicating
//'                     the largest critical value for step-up procedures; if
//'                     `NULL` (the default), it is computed automatically,
//'                     otherwise it needs to be computed manually by the user;
//'                     ignored if `stepUp = FALSE`.
//' @param alpha        single real number strictly between 0 and 1 indicating
//'                     the target FDR level; for `*_fast` kernels, it is only
//'                     needed, if `stepUp = TRUE`.
//' @param support      numeric vector, sorted in increasing order, that
//'                     contains the entirety of all observable values of the
//'                     p-value supports; for `*_fast` kernels, it is ignored if
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
//' For `kernel_DBH_fast()`, `kernel_ADBH_fast()` and `kernel_DBR_fast()`, a
//' vector of transformed p-values is returned. `kernel_DBH_crit`,
//' `kernel_ADBH_crit` and `kernel_DBR_crit` return a list with critical
//' constants (`$crit.consts`) and transformed p-values (`$pval.transf`), but if
//' `stepUp = FALSE`, there are critical values only.
//' 
//' @seealso
//' [`discrete.BH()`], [`direct.discrete.BH()`], [`DBR()`]
//'
struct tau_m_results {
  double value;
  int index;
  std::vector<double> eval;
};

tau_m_results DBH_tau_m(
  const NumericVector* sfuns,
  const NumericVector& CDFcounts,
  const int numCDF,
  const int* lens,
  const NumericVector& support,
  const int numTests,
  const double alpha
);

//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DBH_fast(
  const List& pCDFlist,
  const NumericVector& pvalues,
  const bool stepUp = false,
  const Nullable<NumericVector> tau_max = R_NilValue,
  const double alpha = 0.05,
  const NumericVector& support = NumericVector(),
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
List kernel_DBH_crit(
  const List& pCDFlist,
  const NumericVector& support,
  const NumericVector& sorted_pv,
  const bool stepUp = false,
  const double alpha = 0.05,
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_ADBH_fast(
  const List& pCDFlist,
  const NumericVector& sorted_pv,
  const bool stepUp = false,
  const double alpha = 0.05,
  const NumericVector& support = NumericVector(),
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
List kernel_ADBH_crit(
  const List& pCDFlist,
  const NumericVector& support,
  const NumericVector& sorted_pv,
  const bool stepUp = false,
  const double alpha = 0.05,
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DBR_fast(
  const List& pCDFlist,
  const NumericVector& sorted_pv,
  const double lambda = 0.05,
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
List kernel_DBR_crit(
  const List& pCDFlist,
  const NumericVector& support,
  const NumericVector& sorted_pv,
  const double lambda = 0.05,
  const double alpha = 0.05,
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DBY_fast(
  const List& pCDFlist,
  const NumericVector& pvalues,
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);

//' @rdname kernel
// [[Rcpp::export]]
List kernel_DBY_crit(
  const List& pCDFlist,
  const NumericVector& support,
  const NumericVector& sorted_pv,
  const double alpha = 0.05,
  const Nullable<NumericVector>& pCDFcounts = R_NilValue
);
