#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP cdclogis(SEXP y, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP give_log)
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
      rvalptr[i] = plogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, *give_logptr);
    } else {
      if(yptr[i] >= rightptr[i]){
        rvalptr[i] = plogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, *give_logptr);
      } else {
        if(*give_logptr == 0) {
          rvalptr[i] = dlogis((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i]; 
        } else {
          rvalptr[i] = dlogis((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1) - log(sigmaptr[i]);
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cpclogis(SEXP q, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP lower_tail, SEXP log_p)
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
          rvalptr[i] = plogis((qptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, *log_pptr); 
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
          rvalptr[i] = plogis((qptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, *log_pptr); 
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP sclogis_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
      ddist = dlogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = plogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = ddist / pdist;
      rvalptr[i] = -1 * mills;
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dlogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = plogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        mills = ddist / pdist;
        rvalptr[i] = mills;
      } else {
        rvalptr[i] = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP sclogis_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
      ddist = dlogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = plogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = ddist / pdist;
      rvalptr[i] = -1 * mills * (leftptr[i] - muptr[i]) / sigmaptr[i];
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dlogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = plogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        mills = ddist / pdist;
        rvalptr[i] = mills * (rightptr[i] - muptr[i]) / sigmaptr[i];
      } else {
        rvalptr[i] = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))*(xptr[i]-muptr[i])/pow(sigmaptr[i], 2) - 1/sigmaptr[i];
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP hclogis_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
      ddist = dlogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = plogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      sdist = (1 - 2 * plogis(-(leftptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      mills = ddist / pdist;
      rvalptr[i] = - sdist * mills - pow(mills, 2.0);
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dlogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = plogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        sdist = (1 - 2 * plogis(-(rightptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
        mills = ddist / pdist;
        rvalptr[i] = sdist * mills - pow(mills, 2.0);
      } else {
        rvalptr[i] = - 2/pow(sigmaptr[i], 2) * dlogis((xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP hclogis_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double ddist, pdist, sdist, sdist2, mills, sd2, dcm, dcm2;

  for(i = 0; i < n; i++) {
    if(xptr[i] <= leftptr[i]) {
      ddist = dlogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = plogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      sd2 = pow(sigmaptr[i], 2.0);
      dcm = leftptr[i] - muptr[i];
      dcm2 = pow(dcm, 2.0);
      sdist = (1 - 2 * plogis(-(leftptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      mills = ddist / pdist;
      rvalptr[i] = (2 * dcm/sd2 - dcm2/sd2*sdist)*mills- pow(mills, 2.0)*dcm2/sd2;
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dlogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = plogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        sd2 = pow(sigmaptr[i], 2.0);
        dcm = rightptr[i] - muptr[i];
        dcm2 = pow(dcm, 2.0);
        sdist = (1 - 2 * plogis(-(rightptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
        mills = ddist / pdist;
        rvalptr[i] = (- 2 * dcm/sd2 + dcm2/sd2*sdist)*mills - pow(mills, 2.0)*dcm2/sd2;
      } else {
        sdist = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
        sdist2 = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))*(xptr[i]-muptr[i])/pow(sigmaptr[i], 2) - 1/sigmaptr[i];
        dcm = xptr[i] - muptr[i];
        sd2 = pow(sigmaptr[i], 2.0);
        rvalptr[i] = - sdist*dcm/sd2 - 2 * pow(dcm/sd2, 2) * dlogis(dcm / sigmaptr[i], 0.0, 1.0, 0) - sdist2/sigmaptr[i];
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP hclogis_musigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
      ddist = dlogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = plogis((leftptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 0);
      dcm = leftptr[i] - muptr[i];
      sdist = (1 - 2 * plogis(-dcm/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      mills = ddist / pdist;
      rvalptr[i] = (1/sigmaptr[i] - dcm/sigmaptr[i]*sdist) * mills - dcm/sigmaptr[i] * pow(mills, 2.0);
    } else {
      if(xptr[i] >= rightptr[i]){
        ddist = dlogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
        pdist = plogis((rightptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0, 0);
        dcm = rightptr[i] - muptr[i];
        sdist = (1 - 2 * plogis(-dcm/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
        mills = ddist / pdist;
        rvalptr[i] = (- 1/sigmaptr[i] + dcm/sigmaptr[i]*sdist) * mills - dcm/sigmaptr[i] * pow(mills, 2.0);
      } else {
        dcm = xptr[i] - muptr[i];
        sdist = (1 - 2 * plogis(-dcm/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
        rvalptr[i] = -sdist/sigmaptr[i] - 2*dcm/pow(sigmaptr[i], 3)*dlogis(dcm/sigmaptr[i], 0.0, 1.0, 0);
      }
    }
  }

  UNPROTECT(1);
  return rval;
}




