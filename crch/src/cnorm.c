#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP cdcnorm(SEXP y, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP give_log)
{
  int i, n = length(y);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);
  int *give_logptr = INTEGER(give_log);

  for(i = 0; i < n; i++) {
    if(yptr[i] <= leftptr[i]) {
      rvalptr[i] = pnorm5((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, *give_logptr);
    } else {
      if(yptr[i] >= rightptr[i]){
        rvalptr[i] = pnorm5((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, *give_logptr);
      } else {
        if(*give_logptr == 0) {
          rvalptr[i] = dnorm((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i]; 
        } else {
          rvalptr[i] = dnorm((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1) - log(sigmaptr[i]);
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cpcnorm(SEXP q, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP lower_tail, SEXP log_p)
{
  int i, n = length(q);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *qptr = REAL(q);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);
  int *lower_tailptr = INTEGER(lower_tail);
  int *log_pptr = INTEGER(log_p);
  
  if(*lower_tailptr == 1) {
    for(i = 0; i < n; i++) {
     if(qptr[i] < leftptr[i]) {
        if(*log_pptr == 1) {
          rvalptr[i] = log(0);
        } else {
          rvalptr[i] = 0;
        }
      } else {
        if(qptr[i] >= rightptr[i]){
          rvalptr[i] = 1 * (1 - *log_pptr);
        } else {
          rvalptr[i] = pnorm5((qptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, *log_pptr); 
        }
      }
    }
  } else {
    for(i = 0; i < n; i++) {
      if(qptr[i] <= leftptr[i]) {
        rvalptr[i] = 1 * (1 - *log_pptr);
      } else {
        if(qptr[i] > rightptr[i]){
          if(*log_pptr == 1) {
            rvalptr[i] = log(0);
          } else {
            rvalptr[i] = 0;
          }
        } else {
          rvalptr[i] = pnorm5((qptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, *log_pptr); 
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP scnorm_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double ddist, pdist, mills;

  for(i = 0; i < n; i++) {
    if(xptr[i] <= leftptr[i]) {
      ddist = dnorm((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = ddist / pdist;
      rvalptr[i] = -1 * mills;
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dnorm((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = pnorm5((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        mills = ddist / pdist;
        rvalptr[i] = mills;
      } else {
        rvalptr[i] = (xptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP scnorm_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double ddist, pdist, mills;

  for(i = 0; i < n; i++) {
    if(xptr[i] <= leftptr[i]) {
      ddist = dnorm((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = ddist / pdist;
      rvalptr[i] = -1 * mills * (leftptr[i] - muptr[i]) / sigmaptr[i];
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dnorm((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = pnorm5((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        mills = ddist / pdist;
        rvalptr[i] = mills * (rightptr[i] - muptr[i]) / sigmaptr[i];
      } else {
        rvalptr[i] = (pow((xptr[i] - muptr[i]), 2.0) - pow(sigmaptr[i], 2.0)) / pow(sigmaptr[i], 3.0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP hcnorm_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double ddist, pdist, sdist, mills;

  for(i = 0; i < n; i++) {
    if(xptr[i] <= leftptr[i]) {
      ddist = dnorm((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      sdist = (leftptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0);
      mills = ddist / pdist;
      rvalptr[i] = - sdist * mills - pow(mills, 2.0);
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dnorm((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = pnorm5((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        sdist = (rightptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0);
        mills = ddist / pdist;
        rvalptr[i] = sdist * mills - pow(mills, 2.0);
      } else {
        rvalptr[i] = -1 / pow(sigmaptr[i], 2.0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP hcnorm_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double ddist, pdist, sdist, mills, sd2, dcm, dcm2;

  for(i = 0; i < n; i++) {
    if(xptr[i] <= leftptr[i]) {
      ddist = dnorm((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      sd2 = pow(sigmaptr[i], 2.0);
      dcm = leftptr[i] - muptr[i];
      dcm2 = pow(dcm, 2.0);
      sdist = dcm / sd2;
      mills = ddist / pdist;
      rvalptr[i] = (2 * dcm/sd2 - dcm2/sd2*sdist)*mills- pow(mills, 2.0)*dcm2/sd2;
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dnorm((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = pnorm5((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        sd2 = pow(sigmaptr[i], 2.0);
        dcm = rightptr[i] - muptr[i];
        dcm2 = pow(dcm, 2.0);
        sdist = dcm / sd2;
        mills = ddist / pdist;
        rvalptr[i] = (- 2 * dcm/sd2 + dcm2/sd2*sdist)*mills - pow(mills, 2.0)*dcm2/sd2;
      } else {
        rvalptr[i] = (pow(sigmaptr[i], 2.0) - 3 * pow((xptr[i] - muptr[i]), 2.0))/pow(sigmaptr[i], 4.0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP hcnorm_musigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double ddist, pdist, sdist, mills, dcm;

  for(i = 0; i < n; i++) {
    if(xptr[i] <= leftptr[i]) {
      ddist = dnorm((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      dcm = leftptr[i] - muptr[i];
      sdist = dcm / pow(sigmaptr[i], 2.0);
      mills = ddist / pdist;
      rvalptr[i] = (1/sigmaptr[i] - dcm/sigmaptr[i]*sdist) * mills - dcm/sigmaptr[i] * pow(mills, 2.0);
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dnorm((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = pnorm5((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        dcm = rightptr[i] - muptr[i];
        sdist = dcm / pow(sigmaptr[i], 2.0);
        mills = ddist / pdist;
        rvalptr[i] = (- 1/sigmaptr[i] + dcm/sigmaptr[i]*sdist) * mills - dcm/sigmaptr[i] * pow(mills, 2.0);
      } else {
        rvalptr[i] = -2 * (xptr[i] - muptr[i]) / pow(sigmaptr[i], 3.0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}




