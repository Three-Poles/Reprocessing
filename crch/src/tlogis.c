#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP cdtlogis(SEXP y, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP give_log)
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
      denom = plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
      if(*give_logptr == 0) {
        rvalptr[i] = dlogis((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i]/denom; 
      } else {
        rvalptr[i] = dlogis((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1) - log(sigmaptr[i]*denom);
      }
    }
  }


  UNPROTECT(1);
  return rval;
}


SEXP cptlogis(SEXP q, SEXP mu, SEXP sigma, SEXP left, SEXP right, SEXP lower_tail, SEXP log_p)
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
            qtmp = -(plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
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
            qtmp = -(plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
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
            qtmp = (plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
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
            qtmp = (plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((qptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0));
            denom = plogis((rightptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0) - 
              plogis((leftptr[i] - muptr[i])/sigmaptr[i], 0.0, 1.0, 1,0);
            rvalptr[i] = qtmp/denom;
          }
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP stlogis_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
        denom = (plogis(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - plogis(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
        sdist = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
        enum1 = (dlogis(drm/sigmaptr[i], 0.0, 1.0, 0) - dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
        rvalptr[i] = sdist + enum1/denom;
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP stlogis_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      sdist = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))*(xptr[i]-muptr[i])/sd2 - 1/sigmaptr[i];
      denom = (plogis(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - plogis(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      if(finite(rightptr[i])) {
        enum1 = drm*dlogis(drm/sigmaptr[i], 0.0, 1.0, 0);
      } else {
        enum1 = 0;
      }
      if(finite(leftptr[i])) {
        enum2 = dlm*dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0);
      } else {
        enum2 = 0;
      }
      rvalptr[i] = sdist + (enum1-enum2)/sd2/denom;

    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htlogis_mu(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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
        sdistr = (1 - 2 * plogis(-drm/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      } else {
        sdistr = 0;
      }
      if(finite(dlm)) {
        sdistl = (1 - 2 * plogis(-dlm/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      } else {
        sdistl = 0;
      }
      hdist = - 2/sd2 * dlogis((xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 0);
      denom = (plogis(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - plogis(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      enum1 = (dlogis(drm/sigmaptr[i], 0.0, 1.0, 0) - dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
      enum3 = sdistr*dlogis(drm/sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i] - sdistl*dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0)/sigmaptr[i];
      rvalptr[i] = hdist + pow(enum1/denom, 2.0) + enum3/denom;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htlogis_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdist, sdist2, dcm, sdistl, sdistr, denom, enum2, enum4, hdist, dlm, drm, dlm2, drm2, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = (1 - 2 * plogis(-drm/sigmaptr[i], 0.0, 1.0, 1, 0))*drm/sd2 - 1/sigmaptr[i];
        drm2 = drm;
      } else {
        sdistr = 0;
        drm2 = 0;
      }
      if(finite(dlm)) {
        sdistl = (1 - 2 * plogis(-dlm/sigmaptr[i], 0.0, 1.0, 1, 0))*dlm/sd2 - 1/sigmaptr[i];
        dlm2 = dlm;
      } else {
        sdistl = 0;
        dlm2 = 0;
      }
      sdist = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      sdist2 = (1 - 2 * plogis(-(xptr[i]-muptr[i])/sigmaptr[i], 0.0, 1.0, 1, 0))*(xptr[i]-muptr[i])/pow(sigmaptr[i], 2) - 1/sigmaptr[i];
      dcm = xptr[i] - muptr[i];
      hdist = - sdist*dcm/sd2 - 2 * pow(dcm/sd2, 2) * dlogis(dcm / sigmaptr[i], 0.0, 1.0, 0) - sdist2/sigmaptr[i];
      denom = (plogis(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - plogis(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      enum2 = (drm2*dlogis(drm/sigmaptr[i], 0.0, 1.0, 0) - dlm2*dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0))/sd2;
      enum4 = drm2/sd2*dlogis(drm/sigmaptr[i], 0.0, 1.0, 0)*(sdistr - 1/sigmaptr[i]) - 
             dlm2/sd2*dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0)*(sdistl - 1/sigmaptr[i]);
      rvalptr[i] = hdist + pow(enum2/denom, 2.0) + enum4/denom;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htlogis_musigma(SEXP x, SEXP mu, SEXP sigma, SEXP left, SEXP right)
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

  double sdistl, sdistr, sdist, dcm, denom, enum1, enum2, enum5, hdist, dlm, drm, dlm2, drm2, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = (1 - 2 * plogis(-drm/sigmaptr[i], 0.0, 1.0, 1, 0))*drm/sd2 - 1/sigmaptr[i];
        drm2 = drm;
      } else {
        sdistr = 0;
        drm2 = 0;
      }
      if(finite(dlm)) {
        sdistl = (1 - 2 * plogis(-dlm/sigmaptr[i], 0.0, 1.0, 1, 0))*dlm/sd2 - 1/sigmaptr[i];
        dlm2 = dlm;
      } else {
        sdistl = 0;
        dlm2 = 0;
      }
      dcm = xptr[i] - muptr[i];
      sdist = (1 - 2 * plogis(-dcm/sigmaptr[i], 0.0, 1.0, 1, 0))/sigmaptr[i];
      hdist = -sdist/sigmaptr[i] - 2*dcm/pow(sigmaptr[i], 3)*dlogis(dcm/sigmaptr[i], 0.0, 1.0, 0);
      denom = (plogis(drm/sigmaptr[i], 0.0, 1.0, 1, 0) - plogis(dlm/sigmaptr[i], 0.0, 1.0, 1, 0));
      enum1 = (dlogis(drm/sigmaptr[i], 0.0, 1.0, 0) - dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
      enum2 = (drm2*dlogis(drm/sigmaptr[i], 0.0, 1.0, 0) - dlm2*dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0))/sd2;
      enum5 = (sdistr*dlogis(drm/sigmaptr[i], 0.0, 1.0, 0) - sdistl*dlogis(dlm/sigmaptr[i], 0.0, 1.0, 0))/sigmaptr[i];
      rvalptr[i] = hdist + enum5/denom + enum1*enum2/pow(denom, 2.0);
    }
  }


  UNPROTECT(1);
  return rval;
}




