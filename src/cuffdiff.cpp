#include "cuffdiff.h"

SEXP cuffdiff_wrapper(SEXP a, SEXP b, SEXP c, SEXP d) {
	using namespace Rcpp;
  // step 0: convert input to C++ types
  double RPM_A = as<double>(a), RPMvar_A = as<double>(b), RPM_B = as<double>(c), RPMvar_B = as<double>(d);
  // step 1: call the underlying C++ function
  double res = cuffdiff( RPM_A, RPMvar_A, RPM_B, RPMvar_B );
  // step 2: return the result as a SEXP
  return wrap( res );
}
