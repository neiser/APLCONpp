// main routines

void c_aplcon_aplcon(const int NVAR, const int MCST) {
    aplcon_(&NVAR, &MCST);
}

void c_aplcon_aploop(double X[], double VX[], double F[], int* IRET) {
    aploop_(X, VX, F, IRET);
}

// routines to obtain results

void c_aplcon_chndpv(float* CHI2, int* ND, float* PVAL) {
    chndpv_(CHI2, ND, PVAL);
}

void c_aplcon_apstat(double* FOPT, int* NFUN, int* NITER) {
   apstat_(FOPT, NFUN, NITER);
}

void c_aplcon_appull(double* PULLS) {
    appull_(PULLS);
}

// setup/config routines

void c_aplcon_aprint(const int LUNP, const int IPR) {
    // empty
}

void c_aplcon_apdeps(const double EPSF) {
    apdeps_(&EPSF);
}

void c_aplcon_apepschi(const double EPSCHI) {
    apepschi_(&EPSCHI);
}

void c_aplcon_apderf(const double DERFAC) {
    apderf_(&DERFAC);
}

void c_aplcon_apderu(const double DERUFC) {
    apderu_(&DERUFC);
}

void c_aplcon_apdlow(const double DERLOW) {
    apdlow_(&DERLOW);
}

void c_aplcon_apiter(const int ITERMX) {
    apiter_(&ITERMX);
}

void c_aplcon_apstep(const int I, const double STEP) {
    apstep_(&I, &STEP);
}

void c_aplcon_apfix(const int I) {
    apfix_(&I);
}

void c_aplcon_aplimt(const int I, const double XLOW, const double XHIG) {
    aplimt_(&I, &XLOW, &XHIG);
}

void c_aplcon_apoiss(const int I) {
    apoiss_(&I);
}

void c_aplcon_apsqrt(const int I) {
    // nothing
}

void c_aplcon_aplogn(const int I) {
    // nothing
}