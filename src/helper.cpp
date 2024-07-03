#include "helper.h"

// fast step function evaluation
// step function is represented by a single numeric vector under the conditions
// a) f(x) = x and b) x is sorted
// this is much faster than passing and evaluating R step function objects
// stepfun_result stepfun(const NumericVector &x, const NumericVector &sfun, const int &start){
//   // index variables and vector lengths
//   int pos = start, size = sfun.length(), len = x.length();
//   // if start is below 0, set it to 0
//   //if(start < 0) pos = 0;
//   // output vector of the same length as 'x'
//   NumericVector out(len, 1.0);
//   // if starting position exceeds size, return last position
//   //if(start >= size) return {out, size};
//   
//   // computing results
//   for(int i = 0; i < len && x[i] < 1; i++){
//     while(pos < size && sfun[pos] <= x[i]) pos++;
//     if(pos) out[i] = sfun[pos - 1];
//     else out[i] = 0;
//   }
//   
//   stepfun_result res = {out, pos};
//   
//   return res;
// }

/*
int stepfun_idx(const double &x, const NumericVector &sfun, const int &start = 0){
  // if x is 0, then the step function must be 0, too, i.e. return position -1
  if(x <= 0) return -1;
  else{
  // index variables and vector lengths
    int pos = start, size = sfun.length();
    // if start is below 0, set it to 0
    if(start < 0) pos = 0;
    // if starting position exceeds size, return last position
    if(start >= size) return size - 1;
   
    // compute position
    while(pos <= size - 1 & sfun[pos] <= x && sfun[pos] <= 1) pos++;
    // return position
    return pos - 1;
  }
}
*/

// shortcut function that eliminates all values of a SORTED vector that
// are < limit, except the largest value <= limit
// NumericVector short_eff(const NumericVector &x, const double limit){
//   // length of the vector
//   int len = x.length();
//   // identify values <= limit
//   NumericVector out = x[x <= limit];
//   // eliminate values, but keep their maximum
//   out = x[Range(which_max(out), len - 1)];
// 
//   return out;
// }

// double search_pv(const NumericVector &vec, const double &val, int &pos, bool reverse){
//   if(val < 1){
//     if(!reverse){
//       int len = vec.length();
//       while(pos < len && vec[pos] <= 1 && vec[pos] <= val) pos++;
//       if(pos) return vec[pos - 1]; else return 0;
//     } else {
//       while(pos > 0 && vec[pos] > val) pos--;
//       if(vec[pos] <= val) return vec[pos];
//       else return 0;
//     }
//   } else return 1;
// }

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
// void colsortdec(NumericMatrix &mat){
//   // intermediate vector to store column values (necessary!)
//   NumericVector vec;
//   for(int i = 0; i < mat.ncol(); i++){
//     // store values in a vector
//     vec = NumericVector(mat(_, i));
//     // sort values in ascending order
//     std::sort(vec.begin(), vec.end(), std::greater<double>());
//     // write sorted values back to column
//     mat(_, i) = vec;
//   }
// }

// sort order
IntegerVector order(const NumericVector &x, bool descending){
  arma::vec y = as<arma::vec>(x);
  IntegerVector ord = as<IntegerVector>(wrap(arma::sort_index(y)));

  if(descending)
    return rev(ord);
  else
    return ord;
}
