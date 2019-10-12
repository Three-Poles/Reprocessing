#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cdclogis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdcnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdct(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdtlogis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdtnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdtt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpclogis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpcnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpct(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cptlogis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cptnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cptt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hclogis_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hclogis_musigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hclogis_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hcnorm_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hcnorm_musigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hcnorm_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hct_mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hct_musigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hct_sigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htlogis_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htlogis_musigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htlogis_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htnorm_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htnorm_musigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htnorm_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htt_mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htt_musigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP htt_sigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mycov(SEXP, SEXP);
extern SEXP sclogis_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sclogis_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scnorm_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scnorm_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sct_mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sct_sigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stlogis_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stlogis_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stnorm_mu(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stnorm_sigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stt_mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stt_sigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cdclogis",        (DL_FUNC) &cdclogis,        6},
    {"cdcnorm",         (DL_FUNC) &cdcnorm,         6},
    {"cdct",            (DL_FUNC) &cdct,            7},
    {"cdtlogis",        (DL_FUNC) &cdtlogis,        6},
    {"cdtnorm",         (DL_FUNC) &cdtnorm,         6},
    {"cdtt",            (DL_FUNC) &cdtt,            7},
    {"cpclogis",        (DL_FUNC) &cpclogis,        7},
    {"cpcnorm",         (DL_FUNC) &cpcnorm,         7},
    {"cpct",            (DL_FUNC) &cpct,            8},
    {"cptlogis",        (DL_FUNC) &cptlogis,        7},
    {"cptnorm",         (DL_FUNC) &cptnorm,         7},
    {"cptt",            (DL_FUNC) &cptt,            8},
    {"hclogis_mu",      (DL_FUNC) &hclogis_mu,      5},
    {"hclogis_musigma", (DL_FUNC) &hclogis_musigma, 5},
    {"hclogis_sigma",   (DL_FUNC) &hclogis_sigma,   5},
    {"hcnorm_mu",       (DL_FUNC) &hcnorm_mu,       5},
    {"hcnorm_musigma",  (DL_FUNC) &hcnorm_musigma,  5},
    {"hcnorm_sigma",    (DL_FUNC) &hcnorm_sigma,    5},
    {"hct_mu",          (DL_FUNC) &hct_mu,          6},
    {"hct_musigma",     (DL_FUNC) &hct_musigma,     6},
    {"hct_sigma",       (DL_FUNC) &hct_sigma,       6},
    {"htlogis_mu",      (DL_FUNC) &htlogis_mu,      5},
    {"htlogis_musigma", (DL_FUNC) &htlogis_musigma, 5},
    {"htlogis_sigma",   (DL_FUNC) &htlogis_sigma,   5},
    {"htnorm_mu",       (DL_FUNC) &htnorm_mu,       5},
    {"htnorm_musigma",  (DL_FUNC) &htnorm_musigma,  5},
    {"htnorm_sigma",    (DL_FUNC) &htnorm_sigma,    5},
    {"htt_mu",          (DL_FUNC) &htt_mu,          6},
    {"htt_musigma",     (DL_FUNC) &htt_musigma,     6},
    {"htt_sigma",       (DL_FUNC) &htt_sigma,       6},
    {"mycov",           (DL_FUNC) &mycov,           2},
    {"sclogis_mu",      (DL_FUNC) &sclogis_mu,      5},
    {"sclogis_sigma",   (DL_FUNC) &sclogis_sigma,   5},
    {"scnorm_mu",       (DL_FUNC) &scnorm_mu,       5},
    {"scnorm_sigma",    (DL_FUNC) &scnorm_sigma,    5},
    {"sct_mu",          (DL_FUNC) &sct_mu,          6},
    {"sct_sigma",       (DL_FUNC) &sct_sigma,       6},
    {"stlogis_mu",      (DL_FUNC) &stlogis_mu,      5},
    {"stlogis_sigma",   (DL_FUNC) &stlogis_sigma,   5},
    {"stnorm_mu",       (DL_FUNC) &stnorm_mu,       5},
    {"stnorm_sigma",    (DL_FUNC) &stnorm_sigma,    5},
    {"stt_mu",          (DL_FUNC) &stt_mu,          6},
    {"stt_sigma",       (DL_FUNC) &stt_sigma,       6},
    {NULL, NULL, 0}
};

void R_init_crch(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

