#include "APLCON.h"

typedef double doublereal;
typedef int integer;
typedef float real;

#include "../../APLCON/aploop.f90.h"

// main routines

void c_aplcon_aplcon(int NVAR, int MCST) {
    aplcon::aplcon_(NVAR, MCST);
}

int c_aplcon_aploop(std::vector<double>& X, std::vector<double>& VX, double F[]) {
    vecdr x(X);
    vecdr vx(VX);
    return aplcon::aploop_(x, vx, F);
}

// routines to obtain results

void c_aplcon_chndpv(double* CHI2, int* ND, double* PVAL) {
    aplcon::chndpv_(CHI2, ND, PVAL);
}

void c_aplcon_apstat(double* FOPT, int* NFUN, int* NITER) {
    aplcon::apstat_(FOPT, NFUN, NITER);
}

void c_aplcon_appull(double* PULLS) {
    aplcon::appull_(PULLS);
}

// setup/config routines

void c_aplcon_apdeps(double EPSF) {
   aplcon::apdeps_(&EPSF);
}

void c_aplcon_apstep(int I, double STEP) {
    aplcon::apstep_(&I, &STEP);
}

void c_aplcon_apoiss(int I) {
    aplcon::apoiss_(&I);
}