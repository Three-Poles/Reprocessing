#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP mycov(SEXP y, SEXP x, SEXP nrow, SEXP ncol)
{
  int i, j, n = length(y);
  int k;
  k = length(x)/n;

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, k));

  double *rvalptr = REAL(rval);
  double *yptr = REAL(y);
  double *xptr = REAL(x);

  for(j = 0; j < k; j++) {
    rvalptr[j] = 0;
    for(i = 0; i < n; i++) {
      rvalptr[j] = rvalptr[j] + yptr[i]*xptr[j*n + i];
    }
    rvalptr[j] = rvalptr[j]/n;
  }
  
  UNPROTECT(1);
  return rval;
}

