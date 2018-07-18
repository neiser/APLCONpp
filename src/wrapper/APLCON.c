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

void c_aplcon_apepschi(double EPSCHI) {
    extern int apepschi_(doublereal*);
    apepschi_(&EPSCHI);
}

void c_aplcon_apderf(double DERFAC) {
    extern int apderf_(doublereal*);
    apderf_(&DERFAC);
}

void c_aplcon_apderu(double DERUFC) {
    extern int apderu_(doublereal*);
    apderu_(&DERUFC);
}

void c_aplcon_apdlow(double DERLOW) {
    extern int apdlow_(doublereal*);
    apdlow_(&DERLOW);
}

void c_aplcon_apiter(int ITERMX) {
    extern int apiter_(integer*);
    apiter_(&ITERMX);
}

void c_aplcon_apstep(int I, double STEP) {
    extern int apstep_(integer*, doublereal*);
    apstep_(&I, &STEP);
}

void c_aplcon_apfix(int I) {
    extern int apfix_(integer*);
    apfix_(&I);
}

void c_aplcon_aplimt(int I, double XLOW, double XHIG) {
    extern int aplimt_(integer*, doublereal*, doublereal*);
    aplimt_(&I, &XLOW, &XHIG);
}

void c_aplcon_apoiss(int I) {
    extern int apoiss_(integer*);
    apoiss_(&I);
}

void c_aplcon_apsqrt(int I) {
    // nothing
}

void c_aplcon_aplogn(int I) {
    // nothing
}