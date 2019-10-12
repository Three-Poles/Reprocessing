#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP cdtt(SEXP y, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right, SEXP give_log)
{
  int i, n = length(y);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
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
      denom = pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1,0);
      if(*give_logptr == 0) {
        rvalptr[i] = dt((yptr[i] - muptr[i]) / sigmaptr[i], *dfptr, 0)/sigmaptr[i]/denom; 
      } else {
        rvalptr[i] = dt((yptr[i] - muptr[i]) / sigmaptr[i], *dfptr, 1) - log(sigmaptr[i]*denom);
      }
    }
  }


  UNPROTECT(1);
  return rval;
}


SEXP cptt(SEXP q, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right, SEXP lower_tail, SEXP log_p)
{
  int i, n = length(q);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *qptr = REAL(q);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
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
            qtmp = -(pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((qptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0));
            denom = pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1,0);
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
            qtmp = -(pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((qptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0));
            denom = pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1,0);
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
            qtmp = (pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((qptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0));
            denom = pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1,0);
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
            qtmp = (pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((qptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0));
            denom = pt((rightptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1, 0) - 
              pt((leftptr[i] - muptr[i])/sigmaptr[i], *dfptr, 1,0);
            rvalptr[i] = qtmp/denom;
          }
        }
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP stt_mu(SEXP x, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
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
        denom = (pt(drm/sigmaptr[i], *dfptr, 1, 0) - pt(dlm/sigmaptr[i], *dfptr, 1, 0));
        sdist = (xptr[i]- muptr[i])/pow(sigmaptr[i],2) * (*dfptr + 1) / (*dfptr + pow((xptr[i]-muptr[i])/sigmaptr[i], 2));
        enum1 = (dt(drm/sigmaptr[i], *dfptr, 0) - dt(dlm/sigmaptr[i], *dfptr, 0))/sigmaptr[i];
        rvalptr[i] = sdist + enum1/denom;
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP stt_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
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
      sdist = pow((xptr[i] - muptr[i]), 2.0)/ pow(sigmaptr[i], 3.0)* (*dfptr + 1)/ (*dfptr +  pow((xptr[i] - muptr[i]), 2.0)/ sd2) - 1/sigmaptr[i];
      denom = (pt(drm/sigmaptr[i], *dfptr, 1, 0) - pt(dlm/sigmaptr[i], *dfptr, 1, 0));
      if(finite(rightptr[i])) {
        enum1 = drm*dt(drm/sigmaptr[i], *dfptr, 0);
      } else {
        enum1 = 0;
      }
      if(finite(leftptr[i])) {
        enum2 = dlm*dt(dlm/sigmaptr[i], *dfptr, 0);
      } else {
        enum2 = 0;
      }
      rvalptr[i] = sdist + (enum1-enum2)/sd2/denom;

    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htt_mu(SEXP x, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
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
        sdistr = drm/sd2 * (*dfptr + 1) / (*dfptr + pow(drm, 2.0)/sd2);
      } else {
        sdistr = 0;
      }
      if(finite(dlm)) {
        sdistl = dlm/sd2 * (*dfptr + 1) / (*dfptr + pow(dlm, 2.0)/sd2);
      } else {
        sdistl = 0;
      }
      hdist = (*dfptr + 1)* (pow((xptr[i] - muptr[i]), 2.0) - *dfptr*sd2)/pow((*dfptr*sd2 + pow((xptr[i] - muptr[i]), 2.0)), 2.0);
      denom = (pt(drm/sigmaptr[i], *dfptr, 1, 0) - pt(dlm/sigmaptr[i], *dfptr, 1, 0));
      enum1 = (dt(drm/sigmaptr[i], *dfptr, 0) - dt(dlm/sigmaptr[i], *dfptr, 0))/sigmaptr[i];
      enum3 = sdistr*dt(drm/sigmaptr[i], *dfptr, 0)/sigmaptr[i] - sdistl*dt(dlm/sigmaptr[i], *dfptr, 0)/sigmaptr[i];
      rvalptr[i] = hdist + pow(enum1/denom, 2.0) + enum3/denom;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htt_sigma(SEXP x, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double dcm, dcm2, sdistl, sdistr, denom, enum2, enum4, hdist, dlm, drm, dlm2, drm2, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = pow(drm, 2.0)/ pow(sigmaptr[i], 3.0)* (*dfptr + 1)/ (*dfptr +  pow(drm, 2.0)/ sd2) - 1/sigmaptr[i];
        drm2 = drm;
      } else {
        sdistr = 0;
        drm2 = 0;
      }
      if(finite(dlm)) {
        sdistl = pow(dlm, 2.0)/ pow(sigmaptr[i], 3.0)* (*dfptr + 1)/ (*dfptr +  pow(dlm, 2.0)/ sd2) - 1/sigmaptr[i];
        dlm2 = dlm;
      } else {
        sdistl = 0;
        dlm2 = 0;
      }
      dcm = xptr[i] - muptr[i];
      dcm2 = pow(dcm, 2.0);
      hdist = dcm2 * (*dfptr + 1) * (-3 * sd2 * *dfptr - dcm2) / (sd2 * pow((*dfptr*sd2 + dcm2), 2.0)) + 1/sd2;
      denom = (pt(drm/sigmaptr[i], *dfptr, 1, 0) - pt(dlm/sigmaptr[i], *dfptr, 1, 0));
      enum2 = (drm2*dt(drm/sigmaptr[i], *dfptr, 0) - dlm2*dt(dlm/sigmaptr[i], *dfptr, 0))/sd2;
      enum4 = drm2/sd2*dt(drm/sigmaptr[i], *dfptr, 0)*(sdistr - 1/sigmaptr[i]) - 
             dlm2/sd2*dt(dlm/sigmaptr[i], *dfptr, 0)*(sdistl - 1/sigmaptr[i]);
      rvalptr[i] = hdist + pow(enum2/denom, 2.0) + enum4/denom;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP htt_musigma(SEXP x, SEXP mu, SEXP sigma, SEXP df, SEXP left, SEXP right)
{
  int i, n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *rvalptr = REAL(rval);
  double *xptr = REAL(x);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *dfptr = REAL(df);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);

  double sdistl, sdistr, dcm, denom, enum1, enum2, enum5, hdist, dlm, drm, dlm2, drm2, sd2;

  for(i = 0; i < n; i++) {
    if((xptr[i] < leftptr[i]) | (xptr[i] > rightptr[i])) {
      rvalptr[i] = 0;
    } else {
      sd2 = pow(sigmaptr[i], 2.0);
      drm = rightptr[i] - muptr[i];
      dlm = leftptr[i] - muptr[i];
      if(finite(drm)) {
        sdistr = pow(drm, 2.0)/ pow(sigmaptr[i], 3.0)* (*dfptr + 1)/ (*dfptr +  pow(drm, 2.0)/ sd2) - 1/sigmaptr[i];
        drm2 = drm;
      } else {
        sdistr = 0;
        drm2 = 0;
      }
      if(finite(dlm)) {
        sdistl = pow(dlm, 2.0)/ pow(sigmaptr[i], 3.0)* (*dfptr + 1)/ (*dfptr +  pow(dlm, 2.0)/ sd2) - 1/sigmaptr[i];
        dlm2 = dlm;
      } else {
        sdistl = 0;
        dlm2 = 0;
      }
      dcm = xptr[i] - muptr[i];
      hdist = - 2* dcm * (*dfptr + 1) *sigmaptr[i] * *dfptr / pow((*dfptr*sd2 + pow(dcm, 2.0)), 2.0);
      denom = (pt(drm/sigmaptr[i], *dfptr, 1, 0) - pt(dlm/sigmaptr[i], *dfptr, 1, 0));
      enum1 = (dt(drm/sigmaptr[i], *dfptr, 0) - dt(dlm/sigmaptr[i], *dfptr, 0))/sigmaptr[i];
      enum2 = (drm2*dt(drm/sigmaptr[i], *dfptr, 0) - dlm2*dt(dlm/sigmaptr[i], *dfptr, 0))/sd2;
      enum5 = (sdistr*dt(drm/sigmaptr[i], *dfptr, 0) - sdistl*dt(dlm/sigmaptr[i], *dfptr, 0))/sigmaptr[i];
      rvalptr[i] = hdist + enum5/denom + enum1*enum2/pow(denom, 2.0);
    }
  }


  UNPROTECT(1);
  return rval;
}




