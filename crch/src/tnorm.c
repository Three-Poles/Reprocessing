#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP cdtnorm(SEXP y, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP give_log)
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

  double denom;

  for(i = 0; i < n; i++) {
    if((yptr[i] < leftptr[i]) | (yptr[i] > rightptr[i])) {
      if(*give_logptr == 0) {
        rvalptr[i] = 0.0;
      } else {
        rvalptr[i] = log(0.0);
      }
    } else {
      denom = pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
      if(*give_logptr == 0) {
        rvalptr[i] = dnorm((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i]/denom; 
      } else {
        rvalptr[i] = dnorm((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1) - log(sigmaptr[i]*denom);
      }
    }
  }


  UNPROTECT(1);
  return rval;
}


SEXP cptnorm(SEXP q, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP lower_tail, SEXP log_p)
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

  double qtmp, denom;
  
  if(*lower_tailptr == 1) {
    if(*log_pptr == 1) {
      for(i = 0; i < n; i++) {
        if(qptr[i] < leftptr[i]) {
          rvalptr[i] = log(0);
        } else {
          if(qptr[i] >= rightptr[i]){
            rvalptr[i] = 0;
          } else {
            qtmp = -(pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
            rvalptr[i] = log(qtmp) - log(denom); 
          }
        } 
      }
    } else {      
      for(i = 0; i < n; i++) { 
        if(qptr[i] < leftptr[i]) {
          rvalptr[i] = 0;
        } else {
          if(qptr[i] >= rightptr[i]){
            rvalptr[i] = 1;
          } else {
            qtmp = -(pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
            rvalptr[i] = qtmp/denom;
          }
        }
      }
    }
  } else {
    if(*log_pptr == 1) {
      for(i = 0; i < n; i++) {
        if(qptr[i] < leftptr[i]) {
          rvalptr[i] = 0;
        } else {
          if(qptr[i] >= rightptr[i]){
            rvalptr[i] = log(0);
          } else {
            qtmp = (pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
            rvalptr[i] = log(qtmp) - log(denom); 
          }
        } 
      }
    } else {      
      for(i = 0; i < n; i++) { 
        if(qptr[i] < leftptr[i]) {
          rvalptr[i] = 1;
        } else {
          if(qptr[i] >= rightptr[i]){
            rvalptr[i] = 0;
          } else {
            qtmp = (pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = pnorm5((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              pnorm5((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
            rvalptr[i] = qtmp/denom;
          }
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP stnorm_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdist, enum1, denom, drm, dlm;

  for(i = 0; i < n; i++) {
    if(xptr[i] < leftptr[i]) {
      rvalptr[i] = 0;
    } else {
      if(xptr[i] > rightptr[i]){
        rvalptr[i] = 0;
      } else {
        drm = rightptr[i] - muptr[i];
        dlm = leftptr[i] - muptr[i];
        denom = (pnorm(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - pnorm(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
        sdist = (xptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0);
        enum1 = (dnorm(drm/sigmaptr[i], 0.0, 1.0, 0) - dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
        rvalptr[i] = sdist + enum1/denom;
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP stnorm_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdist, sd2, drm, dlm, enum1, enum2, denom;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sdist = (pow((xptr[i] - muptr[i]), 2.0) - pow(sigmaptr[i], 2.0)) / pow(sigmaptr[i], 3.0);
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      denom = (pnorm(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - pnorm(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      if(finite(rightptr[i])) {
        enum1 = drm*dnorm(drm/sigmaptr[i], 0.0, 1.0, 0);
      } else {
        enum1 = 0;
      }
      if(finite(leftptr[i])) {
        enum2 = dlm*dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0);
      } else {
        enum2 = 0;
      }
      rvalptr[i] = sdist + (enum1-enum2)/sd2/denom;

    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htnorm_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdistl, sdistr, denom, enum1, enum3, hdist, dlm, drm, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = drm / sd2;
      } else {
        sdistr = 0;
      }
      if(finite(dlm)) {
        sdistl = dlm / sd2;
      } else {
        sdistl = 0;
      }
      hdist = -1/sd2;
      denom = (pnorm(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - pnorm(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      enum1 = (dnorm(drm/sigmaptr[i], 0.0, 1.0, 0) - dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
      enum3 = sdistr*dnorm(drm/sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i] - sdistl*dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i];
      rvalptr[i] = hdist + pow(enum1/denom, 2.0) + enum3/denom;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htnorm_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdistl, sdistr, denom, enum2, enum4, hdist, dlm, drm, dlm2, drm2, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = (pow(drm, 2.0) - sd2) / pow(sigmaptr[i], 3.0);
        drm2 = drm;
      } else {
        sdistr = 0;
        drm2 = 0;
      }
      if(finite(dlm)) {
        sdistl = (pow(dlm, 2.0) - sd2) / pow(sigmaptr[i], 3.0);
        dlm2 = dlm;
      } else {
        sdistl = 0;
        dlm2 = 0;
      }
      hdist = (sd2 - 3 * pow((xptr[i] - muptr[i]), 2.0))/pow(sd2, 2.0);
      denom = (pnorm(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - pnorm(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      enum2 = (drm2*dnorm(drm/sigmaptr[i], 0.0, 1.0, 0) - dlm2*dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0))/sd2;
      enum4 = drm2/sd2*dnorm(drm/sigmaptr[i], 0.0, 1.0, 0)*(sdistr - 1/sigmaptr[i]) - 
             dlm2/sd2*dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0)*(sdistl - 1/sigmaptr[i]);
      rvalptr[i] = hdist + pow(enum2/denom, 2.0) + enum4/denom;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htnorm_musigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdistl, sdistr, denom, enum1, enum2, enum5, hdist, dlm, drm, dlm2, drm2, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = (pow(drm, 2.0) - sd2) / pow(sigmaptr[i], 3.0);
        drm2 = drm;
      } else {
        sdistr = 0;
        drm2 = 0;
      }
      if(finite(dlm)) {
        sdistl = (pow(dlm, 2.0) - sd2) / pow(sigmaptr[i], 3.0);
        dlm2 = dlm;
      } else {
        sdistl = 0;
        dlm2 = 0;
      }
      hdist = -2 * (xptr[i] - muptr[i]) / pow(sigmaptr[i], 3.0);
      denom = (pnorm(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - pnorm(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      enum1 = (dnorm(drm/sigmaptr[i], 0.0, 1.0, 0) - dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
      enum2 = (drm2*dnorm(drm/sigmaptr[i], 0.0, 1.0, 0) - dlm2*dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0))/sd2;
      enum5 = (sdistr*dnorm(drm/sigmaptr[i], 0.0, 1.0, 0) - sdistl*dnorm(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
      rvalptr[i] = hdist + enum5/denom + enum1*enum2/pow(denom, 2.0);
    }
  }


  UNPROTECT(1);
  return rval;
}




