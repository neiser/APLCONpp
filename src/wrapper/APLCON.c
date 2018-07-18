typedef double doublereal;
typedef int integer;
typedef float real;

// main routines

void c_aplcon_aplcon(int NVAR, int MCST) {
    extern int aplcon_(integer *, integer *);
    aplcon_(&NVAR, &MCST);
}

void c_aplcon_aploop(double X[], double VX[], double F[], int* IRET) {
    extern int aploop_(doublereal*, doublereal*, doublereal*, integer*);
    aploop_(X, VX, F, IRET);
}

// routines to obtain results

void c_aplcon_chndpv(float* CHI2, int* ND, float* PVAL) {
    extern int chndpv_(real*, integer *, real *);
    chndpv_(CHI2, ND, PVAL);
}

void c_aplcon_apstat(double* FOPT, int* NFUN, int* NITER) {
    extern int apstat_(doublereal*, integer*, integer*);
    apstat_(FOPT, NFUN, NITER);
}

void c_aplcon_appull(double* PULLS) {
    extern int appull_(doublereal*);
    appull_(PULLS);
}

// setup/config routines

void c_aplcon_apdeps(double EPSF) {
    extern int apdeps_(doublereal*);
    apdeps_(&EPSF);
}

void c_aplcon_apstep(int I, double STEP) {
    extern int apstep_(integer*, doublereal*);
    apstep_(&I, &STEP);
}

void c_aplcon_apoiss(int I) {
    extern int apoiss_(integer*);
    apoiss_(&I);
}