#include "helper.h"

// function that binds two vectors, sorts it and eliminates duplications 
NumericVector sort_combine(const NumericVector& x, const NumericVector& y) {
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

// sort order
IntegerVector order(const NumericVector& x, bool descending){
  arma::vec y = as<arma::vec>(x);
  IntegerVector ord = as<IntegerVector>(wrap(arma::sort_index(y)));

  if(descending)
    return rev(ord);
  else
    return ord;
}
