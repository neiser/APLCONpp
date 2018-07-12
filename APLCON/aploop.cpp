#include <fem.hpp> // Fortran EMulation library of fable module

namespace aplcon {

using namespace fem::major_types;

struct common_simcom
{
  double epsf;
  double epschi;
  double chisq;
  double ftest;
  double ftestp;
  double chsqp;
  double frms;
  double frmsp;
  double derfac;
  double decxp;
  double derufc;
  double derlow;
  double weight;
  int nx;
  int nf;
  int num;
  int iflg;
  int init;
  int ipak;
  int indcf;
  int istat;
  int indst;
  int indlm;
  int ndenda;
  int ndende;
  int ndactl;
  int indas;
  int indvs;
  int indtr;
  int indfc;
  int indhh;
  int indxs;
  int inddx;
  int indxp;
  int indrh;
  int indwm;
  int india;
  int ndtot;
  int icnt;
  int nxf;
  int mxf;
  int ndf;
  int iunph;
  int ncst;
  int iter;
  int ncalls;
  int ndpdim;
  int indqn;
  int itermx;
  int nauxc;
  int indpu;
  int nfprim;
  int ndtotl;
  arr<float, 2> tab;

  common_simcom() :
    epsf(fem::double0),
    epschi(fem::double0),
    chisq(fem::double0),
    ftest(fem::double0),
    ftestp(fem::double0),
    chsqp(fem::double0),
    frms(fem::double0),
    frmsp(fem::double0),
    derfac(fem::double0),
    decxp(fem::double0),
    derufc(fem::double0),
    derlow(fem::double0),
    weight(fem::double0),
    nx(fem::int0),
    nf(fem::int0),
    num(fem::int0),
    iflg(fem::int0),
    init(fem::int0),
    ipak(fem::int0),
    indcf(fem::int0),
    istat(fem::int0),
    indst(fem::int0),
    indlm(fem::int0),
    ndenda(fem::int0),
    ndende(fem::int0),
    ndactl(fem::int0),
    indas(fem::int0),
    indvs(fem::int0),
    indtr(fem::int0),
    indfc(fem::int0),
    indhh(fem::int0),
    indxs(fem::int0),
    inddx(fem::int0),
    indxp(fem::int0),
    indrh(fem::int0),
    indwm(fem::int0),
    india(fem::int0),
    ndtot(fem::int0),
    icnt(fem::int0),
    nxf(fem::int0),
    mxf(fem::int0),
    ndf(fem::int0),
    iunph(fem::int0),
    ncst(fem::int0),
    iter(fem::int0),
    ncalls(fem::int0),
    ndpdim(fem::int0),
    indqn(fem::int0),
    itermx(fem::int0),
    nauxc(fem::int0),
    indpu(fem::int0),
    nfprim(fem::int0),
    ndtotl(fem::int0),
    tab(dimension(1000, 10), fem::fill0)
  {}
};

struct common_nauxcm
{
  static const int naux = 125000;

  arr<double> aux;

  common_nauxcm() :
    aux(dimension(naux), fem::fill0)
  {}
};

const int common_nauxcm::naux;

struct common :
  fem::common,
  common_simcom,
  common_nauxcm
{
  fem::cmn_sve anumde_sve;
  fem::cmn_sve duminv_sve;
  fem::cmn_sve dgamml_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

//C 1 "aploop.F"
//C 1 "<built-in>"
//C 1 "<command-line>"
//C 1 "aploop.F"
//C
//C dimension parameters
void
aplcon(
  common& cmn,
  int const& nvar,
  int const& mcst)
{
  // COMMON simcom
  int& nx = cmn.nx;
  int& nf = cmn.nf;
  int& indst = cmn.indst;
  int& indlm = cmn.indlm;
  int& indtr = cmn.indtr;
  int& ndtot = cmn.ndtot;
  int& nxf = cmn.nxf;
  // COMMON nauxcm
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  //
  //C     ==================================================================
  //C     initialize and define dimension  NVAR/MCST
  //C     set default parameters
  //C     define pointer to arrays within aux array
  //C     clear aux array
  //C     ==================================================================
  //C ,I1,I2,NFF,I
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 12 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 13 "aploop.F" 2
  //C
  //C number of variables
  nx = nvar;
  //C number of constraint equations
  nf = mcst;
  //C primary value of NF
  cmn.nfprim = nf;
  //C dimension of AUX array
  cmn.ndpdim = naux;
  //C derivative factor
  cmn.derfac = 1.0e-3;
  //C factor or unmeasured variable
  cmn.derufc = 1.0e-5;
  //C factor for lower limit
  cmn.derlow = 1.0e-2;
  //C accuracy limit
  cmn.epsf = 1.0e-6;
  //C chi2 accuracy limit
  cmn.epschi = 0.1e-4;
  //C max number of iterations
  cmn.itermx = 10;
  //C copy AUX dimension
  cmn.nauxc = naux;
  //C
  cmn.init = 0;
  //C init phase
  cmn.istat = 0;
  //C
  //C total number of fit equations
  nxf = nx + nf;
  //C number elements symmetric matrix
  cmn.mxf = (nxf * nxf + nxf) / 2;
  //C
  //C ______________________________________________________________________
  //C     Indices for steps, flags and limits are defined here
  //C reserve NX * (NF + 2) for Jacobian
  //C steps
  indst = nx * (nf + 2);
  //C transformation flags
  indtr = indst + nx;
  //C 2*limits for variables
  indlm = indtr + nx;
  //C space used so far
  ndtot = indlm + 2 * nx;
  if (ndtot > naux) {
    FEM_STOP(0);
  }
  //C     __________________________________________________________________
  //C     storage of initial sub-arrays
  //C           1 ...    NX * (NF+2)   Jacobian, derivative matrix   A(.)
  //C     INDST+1 ...    NX            steps                         ST(.)
  //C     INDTR+1 ...    NX            properties of variables
  //C     INDLM+1 ...    2 * NX        limits for variables          XL(2,.)
  //C                    -----------
  //C     NDTOT =        NX * (NF+6)   initial memory area
  //C NX*(NF+6)
  int ij = fem::int0;
  FEM_DO_SAFE(ij, 1, ndtot) {
    //C clear  A(.), ST(.),...,XL(2,.)
    aux(ij) = 0.0e0;
  }
  //C reset n d f
  cmn.ndf = nf;
  //C reset number of calls
  cmn.ncalls = 0;
  //C
}

struct anumde_save
{
  bool tinue;

  anumde_save() :
    tinue(fem::bool0)
  {}
};

//C
//C numerical derivatives
void
anumde(
  common& cmn,
  arr_ref<double> x,
  arr_cref<double> f,
  arr_ref<double> a,
  arr_ref<double> st,
  arr_cref<double, 2> xl,
  arr_cref<double> fc,
  arr_ref<double> hh,
  int& jret)
{
  FEM_CMN_SVE(anumde);
  x(dimension(star));
  f(dimension(star));
  a(dimension(star));
  st(dimension(star));
  xl(dimension(2, star));
  fc(dimension(star));
  hh(dimension(star));
  int& nx = cmn.nx;
  int& nf = cmn.nf;
  int& ipak = cmn.ipak;
  int& indtr = cmn.indtr;
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  //
  bool& tinue = sve.tinue;
  if (is_called_first_time) {
    tinue = false;
  }
  int i = fem::int0;
  int ntrfl = fem::int0;
  int ntvar = fem::int0;
  int ntmes = fem::int0;
  int ntder = fem::int0;
  int ntine = fem::int0;
  int ntlim = fem::int0;
  int ntprf = fem::int0;
  double xsave = fem::double0;
  int ilr = fem::int0;
  bool limdef = fem::bool0;
  double stm = fem::double0;
  arr_1d<2, double> xd(fem::fill0);
  arr_1d<2, double> xt(fem::fill0);
  int j = fem::int0;
  int ij = fem::int0;
  int nzer = fem::int0;
  int nonz = fem::int0;
  int nalz = fem::int0;
  double ratmax = fem::double0;
  double derzer = fem::double0;
  double der = fem::double0;
  double ratdif = fem::double0;
  //C     ==================================================================
  //C     calculation of numerical derivatives, one variable at a time
  //C        check limits for variable
  //C        calculate displaced value of variable
  //C        calculate derivative in Jacobian matrix A
  //C        classify derivative properties of variable
  //C        return with JRET=-1, unless finished
  //C
  //C     X(.),F(.) = variables and constraints
  //C     A(.)      = Jacobian
  //C     ST(.)     = steps
  //C     XL(2,.)   = limits
  //C     FC(.)     = Constraints at central variable values
  //C     HH(.)     = Constraints at first step (variable + step)
  //C     __________________________________________________________________
  //C     Logic:
  //C     Entry  1: start with I=0 ...
  //C            2: continue with indices
  //C     Return 1: continue with constraint evaluation
  //C            2: Jacobian ready
  //C     ==================================================================
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 317 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 318 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 319 "aploop.F" 2
  //C entry flag
  //C     ...
  //C ...means continue at return
  jret = -1;
  //C continue
  if (tinue) {
    goto statement_30;
  }
  tinue = true;
  //C initialize derivative loop
  i = 0;
  //C finished
  statement_10:
  if (i >= nx) {
    //C ... means differentation finished
    jret = 0;
    //C Jacobian ready
    tinue = false;
    return;
  }
  //C next variable
  i++;
  //C skip fixed variable
  if (st(i) == 0.0f) {
    goto statement_10;
  }
  ipak = i;
  //C
  //C 1 "unpackfl.inc" 1
  //C
  //C     unpackfl.inc = code for flag unpacking
  //C
  //C 60 "unpackfl.inc"
  //C get packed flags
  ntrfl = aux(indtr + ipak);
  //C transformation flag
  ntvar = fem::mod(ntrfl, 10);
  //C M-estimate flag
  ntmes = fem::mod(ntrfl / 10, 10);
  //C derivative type flag
  ntder = fem::mod(ntrfl / 100, 10);
  //C inequality flag
  ntine = fem::mod(ntrfl / 1000, 10);
  //C limit flag
  ntlim = fem::mod(ntrfl / 10000, 10);
  //C profile flag
  ntprf = fem::mod(ntrfl / 100000, 10);
  //C 338 "aploop.F" 2
  //C skip repeated derivative calculation
  if (ntder >= 4) {
    goto statement_10;
  }
  //C save current value of variable
  xsave = x(i);
  //C define steps
  ilr = 0;
  //C     __________________________________________________________________
  //C     check limits for variable
  //C true if limits defined
  limdef = xl(1, i) != xl(2, i);
  if (limdef) {
    if (xsave + st(i) > xl(2, i) || xsave - st(i) > xl(1, i)) {
      //C minimal step s
      stm = 0.9999f * fem::min(xl(2, i) - xsave, xsave - xl(1, i));
      if (3.0f * stm > st(i)) {
        //C use smaller symmetric step
        st(i) = stm;
      }
      else {
        //C minimal ste
        stm = 0.4999f * fem::max(xl(2, i) - xsave, xsave - xl(1, i));
        if (2.0f * stm < st(i)) {
          st(i) = stm;
        }
        if (st(i) < xl(2, i) - xsave) {
          xd(1) = xsave + st(i);
          //C + one-sided steps
          xd(2) = xsave + st(i) * 2.0f;
          ilr = 1;
        }
        else {
          xd(1) = xsave - st(i);
          //C - one-sided steps
          xd(2) = xsave - st(i) * 2.0f;
          ilr = 2;
        }
      }
    }
  }
  //C     __________________________________________________________________
  //C     define displaced values for derivative calculation
  if (ilr == 0) {
    statement_20:
    if (ntvar == 0 || ntvar == 2 || ntvar == 3) {
      //C symmetric (two-sided) steps
      xt(1) = xsave + st(i);
      xt(2) = xsave - st(i);
      xd(1) = xt(1);
      xd(2) = xt(2);
      //C 1/x
    }
    else if (ntvar == 1) {
      //C internal
      xt(1) = 1.0e0 / xsave + st(i);
      xt(2) = 1.0e0 / xsave - st(i);
      //C external
      xd(1) = 1.0e0 / xt(1);
      xd(2) = 1.0e0 / xt(2);
      //C log-normal
    }
    else if (ntvar == 4) {
      //C internal
      xt(1) = fem::log(xsave) + st(i);
      xt(2) = fem::log(xsave) - st(i);
      //C external
      xd(1) = fem::exp(xt(1));
      xd(2) = fem::exp(xt(2));
      //C sqrt
    }
    else if (ntvar == 5) {
      if (fem::pow2(st(i)) >= xsave) {
        if (xsave <= 0.0e0) {
          ntvar = 0;
          //C
          //C 1 "packfl.inc" 1
          //C
          //C     packfl.inc   = code for flag packing
          //C
          //C     explanation see:
          //C     unpackfl.inc = code for flag unpacking
          //C
          aux(indtr + ipak) = ((((ntprf * 10 + ntlim) * 10 + ntine) *
            10 + ntder) * 10 + ntmes) * 10 + ntvar;
          //C 387 "aploop.F" 2
          goto statement_20;
        }
        else {
          st(i) = 0.9e0 * fem::sqrt(xsave);
        }
      }
      //C internal
      xt(1) = fem::sqrt(xsave) + st(i);
      xt(2) = fem::sqrt(xsave) - st(i);
      //C external
      xd(1) = fem::pow2(xt(1));
      xd(2) = fem::pow2(xt(2));
      //C x**power
    }
    else if (ntvar == 6) {
      //C internal
      xt(1) = fem::pow(xsave, xl(2, i)) + st(i);
      xt(2) = fem::pow(xsave, xl(2, i)) - st(i);
      //C external
      xd(1) = fem::pow(xt(1), (1.0e0 / xl(2, i)));
      xd(2) = fem::pow(xt(2), (1.0e0 / xl(2, i)));
    }
  }
  //C     __________________________________________________________________
  //C     set variable to displaced value and return for calculation
  //C first step
  x(i) = xd(1);
  i = -i;
  return;
  //C     __________________________________________________________________
  //C     continue
  //C calculation of first step done ...
  statement_30:
  if (i < 0) {
    FEM_DO_SAFE(j, 1, nf) {
      //C save constraint values
      hh(j) = f(j);
    }
    //C reverse flag
    i = -i;
    //C set next step ...
    x(i) = xd(2);
    //C ... and return for second step
    return;
  }
  //C     __________________________________________________________________
  //C     INIT ne 0: second step done - calculate derivative
  //C restore variable I
  x(i) = xsave;
  //C derivative calculation
  ij = i;
  //C reset: number of zero derivatives
  nzer = 0;
  //C reset: number of non-zero derivatives
  nonz = 0;
  //C flag all derivatives are zero
  nalz = 0;
  //C max of diff-ratio
  ratmax = 0.0e0;
  //C abs of nonzero-derivative
  derzer = 0.0e0;
  //C loop on all constraint functions
  FEM_DO_SAFE(j, 1, nf) {
    //C symmetric formula
    if (ilr == 0) {
      //C!! internal variable
      der = (hh(j) - f(j)) / (xt(1) - xt(2));
      //C asymmetric formula
    }
    else {
      der = 0.5e0 * (3.0e0 * fc(j) + f(j) - 4.0e0 * hh(j)) / st(i);
      //C sign
      if (ilr == 2) {
        der = -der;
      }
    }
    //C      _________________________________________________________________
    //C      classify derivative properties of variable I
    //C non all zero
    if (a(ij) != 0.0e0 || der != 0.0e0) {
      nalz = 1;
    }
    if (der == 0.0e0) {
      //C derivative zero - count
      nzer++;
      //C derivative non-zero
    }
    else {
      ratdif = fem::abs(a(ij) - der) / (fem::abs(a(ij)) + fem::abs(der));
      ratmax = fem::max(ratmax, ratdif);
      //C abs value of this derivative
      derzer = fem::abs(der);
      //C count non-zero derivative
      nonz++;
    }
    //C insert into Jacobian matrix A
    a(ij) = der;
    ij += nx;
  }
  if (ntder == 0) {
    ntder = 1;
  }
  else {
    if (nonz == 1 && fem::abs(derzer - 1.0e0) < 1.0e-12) {
      ntder = fem::min(ntder + 1, 7);
    }
    else if (ratmax < 1.0e-12) {
      ntder = fem::min(ntder + 1, 7);
    }
    else {
      //C reset to 2
      ntder = 2;
    }
  }
  ntder = 0;
  //C
  //C 1 "packfl.inc" 1
  //C
  //C     packfl.inc   = code for flag packing
  //C
  //C     explanation see:
  //C     unpackfl.inc = code for flag unpacking
  //C
  aux(indtr + ipak) = ((((ntprf * 10 + ntlim) * 10 + ntine) * 10 +
    ntder) * 10 + ntmes) * 10 + ntvar;
  //C 461 "aploop.F" 2
  goto statement_10;
}

struct duminv_save
{
  double eps;

  duminv_save() :
    eps(fem::double0)
  {}
};

//C
//C matrix inve
void
duminv(
  common& cmn,
  arr_cref<double> a,
  arr_ref<double> w,
  arr_ref<double> b,
  int const& nx,
  int const& nf,
  int const& mb,
  int& nrank,
  arr_ref<double> aux,
  arr_ref<double> qnext)
{
  FEM_CMN_SVE(duminv);
  a(dimension(star));
  w(dimension(star));
  b(dimension(star));
  aux(dimension(star));
  qnext(dimension(star));
  double& eps = sve.eps;
  if (is_called_first_time) {
    eps = 1.0e-6;
  }
  int n = fem::int0;
  int i = fem::int0;
  int ij = fem::int0;
  int ia = fem::int0;
  int j = fem::int0;
  int jfirst = fem::int0;
  int nmeas = fem::int0;
  int jlast = fem::int0;
  int m = fem::int0;
  double sum = fem::double0;
  int jk = fem::int0;
  int k = fem::int0;
  bool solve = fem::bool0;
  double vkk = fem::double0;
  int last = fem::int0;
  int jj = fem::int0;
  int l = fem::int0;
  int kk = fem::int0;
  int jl = fem::int0;
  double vjk = fem::double0;
  int lk = fem::int0;
  //C     Obtain solution of a system of linear equations V *  X  =  B  with
  //C     symmetric matrix V and inverse (for M =  1)  or  matrix  inversion
  //C     only (for M = 0).
  //C
  //C                   - - - -
  //C        CALL SMINV(W,B,N,M,NRANK,AUX,QNEXT)
  //C                   - -     -----
  //C
  //C           W = symmetric N-by-N matrix in symmetric storage mode
  //C               W(1) = W11, W(2) = W12, W(3) = W22, W(4) = W13, . . .
  //C               replaced by inverse matrix
  //C           B = N-vector   (for M = 0 use a dummy argument)
  //C               replaced by solution vector
  //C           M = see above
  //C
  //C     Method of solution is by elimination selecting the  pivot  on  the
  //C     diagonal each stage. The rank of the matrix is returned in  NRANK.
  //C     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
  //C     matrix V and the corresponding elements of  B  are  set  to  zero.
  //C     SMINV can be used for a dimension up to 100 (see INVCDR).
  //C
  //C     special entry for partial inversion ******************************
  //C
  //C     ...
  nrank = 0;
  //C dimension parameter
  n = nx + nf;
  //C
  //C     make sure AUX is zero, prevents uninit access
  //C     in continued execution of DBMINV
  FEM_DO_SAFE(i, 1, n) {
    aux(i) = 0.0e0;
  }
  //C
  //C     -VX(NX-sym) is already inserted in W(NX+NF-sym) ------------------
  //C
  //C number of elements of V(.)
  ij = (nx * nx + nx) / 2;
  //C
  FEM_DO_SAFE(i, 1, n) {
    //C reset pointer
    qnext(i) = 0.0e0;
  }
  //C
  ia = 0;
  FEM_DO_SAFE(j, 1, nf) {
    FEM_DO_SAFE(i, 1, nx) {
      //C copy A(.) into W_12
      w(ij + i) = a(ia + i);
    }
    FEM_DO_SAFE(i, 1, j) {
      //C reset last submatrix W_22 of W(.)
      w(ij + nx + i) = 0.0f;
    }
    ij += nx + j;
    ia += nx;
  }
  //C
  //C     distinguish between measured and unmeasured variables ------------
  //C
  //C first index of measured variable
  jfirst = 0;
  //C number of measured variables
  nmeas = 0;
  FEM_DO_SAFE(i, 1, nx) {
    //C measured variable
    if (w((i * i + i) / 2) < 0.0f) {
      if (jfirst == 0) {
        //C first index of measured variable
        jfirst = i;
      }
      else {
        //C insert index at previous index
        qnext(jlast) = i;
      }
      //C save index
      jlast = i;
      nmeas++;
    }
  }
  //C nothing to do
  if (jlast == 0) {
    goto statement_10;
  }
  //C stop index for last measured variable
  qnext(jlast) = -1;
  //C
  //C     apply exchange algorithm to sub-matrices -------------------------
  //C
  //C loop I over constraint equations
  FEM_DO_SAFE(i, nx + 1, n) {
    //C
    //C first index of unmeasured variable
    j = jfirst;
    //C already inverted element index J
    FEM_DO_SAFE(m, 1, nmeas) {
      sum = 0.0e0;
      //C index of diagonal element before
      jk = (j * j - j) / 2;
      FEM_DO_SAFE(k, 1, nx) {
        //C index in j column
        if (k <= j) {
          jk++;
        }
        if (qnext(k) != 0.0e0) {
          sum += w(jk) * w((i * i - i) / 2 + k);
        }
        //C index in j row
        if (k >= j) {
          jk += k;
        }
      }
      //C = A-row * VX-row/col
      aux(j) = sum;
      //C next index of unmeasured variable
      j = qnext(j);
    }
    //C
    FEM_DO_SAFE(k, i, n) {
      sum = 0.0e0;
      //C first index of unmeasured variable
      j = jfirst;
      //C already inverted element index J
      FEM_DO_SAFE(m, 1, nmeas) {
        //C = A-row * H
        sum += w((k * k - k) / 2 + j) * aux(j);
        //C next index of unmeasured variable
        j = qnext(j);
      }
      //C add to diagonal W_22
      w((k * k - k) / 2 + i) += sum;
    }
    //C
    //C first index of unmeasured variable
    j = jfirst;
    FEM_DO_SAFE(m, 1, nmeas) {
      //C add to off-diagonal W_22
      w((i * i - i) / 2 + j) = -aux(j);
      //C next index of unmeasured variable
      j = qnext(j);
    }
    //C
  }
  //C
  //C     set pointer for unmeasured variables -----------------------------
  //C
  jfirst = 0;
  jlast = 0;
  FEM_DO_SAFE(i, 1, n) {
    //C unmeasured variable
    if (qnext(i) == 0.0e0) {
      if (jfirst == 0) {
        //C first index of unmeasured variable
        jfirst = i;
      }
      else {
        //C next index of unmeasured variable
        qnext(jlast) = i;
      }
      jlast = i;
    }
    else {
      //C reset index for measured variable
      qnext(i) = 0.0e0;
    }
  }
  //C no unmeasured variable
  if (jlast == 0) {
    goto statement_10;
  }
  //C end flag
  qnext(jlast) = -1;
  //C
  //C     common code for inversion and (M=1) solution of matrix equation
  //C
  statement_10:
  solve = true;
  //C solution flag
  if (mb == 0) {
    solve = false;
  }
  //C
  //C     loop begin (loop on all remaining rows/cols)
  //C
  //C loop on all remaining elements
  FEM_DO_SAFE(i, 1, n) {
    //C search for pivot element
    vkk = 0.0e0;
    //C pivot index
    k = 0;
    //C first candidate index
    j = jfirst;
    last = 0;
    //C test for linearity and zero matrix
    statement_20:
    if (j > 0) {
      //C diagonal index
      jj = (j * j + j) / 2;
      if (fem::abs(w(jj)) > fem::max(fem::abs(vkk), eps * aux(j))) {
        //C largest pivot candidate so far
        vkk = w(jj);
        //C index of largest
        k = j;
        l = last;
      }
      last = j;
      //C index of next candidate
      j = qnext(j);
      goto statement_20;
    }
    //C pivot element found - proceed
    if (k != 0) {
      //C increase rank counter
      nrank++;
      kk = (k * k + k) / 2;
      if (l == 0) {
        //C new first index
        jfirst = qnext(k);
      }
      else {
        //C bridge used index
        qnext(l) = qnext(k);
      }
      //C reset used index
      qnext(k) = 0.0e0;
      //C increase rank
      nrank++;
      //C
      //C invert pivot
      vkk = 1.0f / vkk;
      w(kk) = -vkk;
      if (solve) {
        b(k) = b(k) * vkk;
      }
      jk = kk - k;
      jl = 0;
      //C elimination
      FEM_DO_SAFE(j, 1, n) {
        if (j == k) {
          jk = kk;
          jl += j;
        }
        else {
          if (j < k) {
            jk++;
          }
          else {
            jk += j - 1;
          }
          vjk = w(jk);
          w(jk) = vkk * vjk;
          if (solve) {
            b(j) = b(j) - b(k) * vjk;
          }
          lk = kk - k;
          FEM_DO_SAFE(l, 1, j) {
            jl++;
            if (l == k) {
              lk = kk;
            }
            else {
              if (l < k) {
                lk++;
              }
              else {
                lk += l - 1;
              }
              w(jl) = w(jl) - w(lk) * vjk;
            }
          }
        }
      }
      //C no pivot candadate found - reset
    }
    else {
      FEM_DO_SAFE(k, 1, n) {
        //C undefined variable
        if (qnext(k) != 0.0e0) {
          //C clear undefined vector element
          if (solve) {
            b(k) = 0.0e0;
          }
          FEM_DO_SAFE(j, 1, k) {
            //C clear matri
            if (qnext(j) != 0.0e0) {
              w((k * k - k) / 2 + j) = 0.0e0;
            }
          }
        }
      }
      goto statement_30;
    }
    //C
    //C end of inversion loop
  }
  //C
  statement_30:
  FEM_DO_SAFE(i, 1, (n * n + n) / 2) {
    //C finally reverse sign
    w(i) = -w(i);
  }
  return;
  //C
  //C go to invert remaining parts
  goto statement_10;
}

//C
//C scalar vector product
double
scalxy(
  arr_cref<double> x,
  arr_cref<double> y,
  int const& n)
{
  double return_value = fem::double0;
  x(dimension(star));
  y(dimension(star));
  //C     Scalar product of two vectors
  //C                - - -                 T
  //C        S = VXY(X,Y,N)           S = X  * Y (scalar product)
  //C
  //C I,M
  //C     ...
  double sum = 0.0e0;
  int j = fem::int0;
  FEM_DO_SAFE(j, 1, n) {
    sum += x(j) * y(j);
  }
  return_value = sum;
  return return_value;
}

//C
//C next iteration step
void
aniter(
  common& cmn,
  arr_cref<double> x,
  arr_cref<double> vx,
  arr_cref<double> f,
  arr_cref<double> a,
  arr_ref<double> xp,
  arr_ref<double> rh,
  arr_ref<double> wm,
  arr_ref<double> dx)
{
  x(dimension(star));
  vx(dimension(star));
  f(dimension(star));
  a(dimension(star));
  xp(dimension(star));
  rh(dimension(star));
  wm(dimension(star));
  dx(dimension(star));
  // COMMON simcom
  double& chisq = cmn.chisq;
  double& chsqp = cmn.chsqp;
  double& weight = cmn.weight;
  int& nx = cmn.nx;
  int& nf = cmn.nf;
  int& ipak = cmn.ipak;
  int& indhh = cmn.indhh;
  int& iter = cmn.iter;
  // COMMON nauxcm
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  //
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 468 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 469 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 470 "aploop.F" 2
  //C     ...
  //C start next iteration
  iter++;
  //C save current chi^2
  chsqp = chisq;
  //C
  //C     __________________________________________________________________
  //C     right-hand side of equation
  cmn.ncst = 0;
  //C first NX components
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, nx) {
    //C define right hand side of equation
    rh(i) = 0.0e0;
  }
  //C next NF components
  int j = fem::int0;
  FEM_DO_SAFE(j, 1, nf) {
    aux(cmn.indfc + j) = f(j);
    rh(nx + j) = -f(j);
  }
  int ia = 0;
  //C "subtract" actual step
  FEM_DO_SAFE(j, 1, nf) {
    rh(nx + j) += scalxy(a(ia + 1), dx, nx);
    ia += nx;
    //C right hand side for chi**2
    aux(indhh + j) = rh(nx + j);
  }
  //C     __________________________________________________________________
  //C     form matrix and solve
  FEM_DO_SAFE(i, 1, (nx * nx + nx) / 2) {
    //C copy -VX(.) into W_11
    wm(i) = -vx(i);
  }
  //C modify V for Poisson variables
  int ii = 0;
  int ntrfl = fem::int0;
  int ntvar = fem::int0;
  int ntmes = fem::int0;
  int ntder = fem::int0;
  int ntine = fem::int0;
  int ntlim = fem::int0;
  int ntprf = fem::int0;
  FEM_DO_SAFE(i, 1, nx) {
    ii += i;
    ipak = i;
    //C
    //C 1 "unpackfl.inc" 1
    //C
    //C     unpackfl.inc = code for flag unpacking
    //C
    //C 60 "unpackfl.inc"
    //C get packed flags
    ntrfl = aux(cmn.indtr + ipak);
    //C transformation flag
    ntvar = fem::mod(ntrfl, 10);
    //C M-estimate flag
    ntmes = fem::mod(ntrfl / 10, 10);
    //C derivative type flag
    ntder = fem::mod(ntrfl / 100, 10);
    //C inequality flag
    ntine = fem::mod(ntrfl / 1000, 10);
    //C limit flag
    ntlim = fem::mod(ntrfl / 10000, 10);
    //C profile flag
    ntprf = fem::mod(ntrfl / 100000, 10);
    //C 502 "aploop.F" 2
    //C Poisson
    if (ntvar == 2) {
      //C -MAX(ABS(X(I)),1.0D0)
      wm(ii) = -fem::sqrt(1.0f + fem::pow2(x(i)));
    }
  }
  int nrank = fem::int0;
  arr<double> diag(dimension(1000), fem::fill0);
  arr<double> qnext(dimension(1000), fem::fill0);
  duminv(cmn, a, wm, rh, nx, nf, 1, nrank, diag, qnext);
  //C next chi^2
  chisq = -scalxy(aux(indhh + 1), rh(nx + 1), nf);
  if (chisq < 0.0e0) {
    chisq = 0.0e0;
  }
  //C     __________________________________________________________________
  //C     handle corrections and cutstep
  //C default weight
  weight = 1.0e0;
  if (iter > 1 && chisq >= 2.00e0 * chsqp) {
    weight = 0.1e0;
  }
  if (iter > 1 && chisq >= 3.00e0 * chsqp) {
    weight = 0.05e0;
  }
  //C
  FEM_DO_SAFE(i, 1, nx) {
    //C save previous corrections
    xp(i) = dx(i);
    //C store new corrections
    dx(i) = rh(i);
  }
}

void
addtox(
  common& cmn,
  arr_ref<double> x,
  arr_cref<double> xs,
  arr_ref<double> dx,
  arr_cref<double> xp)
{
  x(dimension(star));
  xs(dimension(star));
  dx(dimension(star));
  xp(dimension(star));
  // COMMON simcom
  double& weight = cmn.weight;
  int& ipak = cmn.ipak;
  // COMMON nauxcm
  const int naux = 125000;
  arr_cref<double> aux(cmn.aux, dimension(naux));
  //
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 525 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 526 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 527 "aploop.F" 2
  //C     ...
  int i = fem::int0;
  int ntrfl = fem::int0;
  int ntvar = fem::int0;
  int ntmes = fem::int0;
  int ntder = fem::int0;
  int ntine = fem::int0;
  int ntlim = fem::int0;
  int ntprf = fem::int0;
  FEM_DO_SAFE(i, 1, cmn.nx) {
    //C reduce step evtl.
    dx(i) = weight * dx(i) + (1.0e0 - weight) * xp(i);
    ipak = i;
    //C
    //C 1 "unpackfl.inc" 1
    //C
    //C     unpackfl.inc = code for flag unpacking
    //C
    //C 60 "unpackfl.inc"
    //C get packed flags
    ntrfl = aux(cmn.indtr + ipak);
    //C transformation flag
    ntvar = fem::mod(ntrfl, 10);
    //C M-estimate flag
    ntmes = fem::mod(ntrfl / 10, 10);
    //C derivative type flag
    ntder = fem::mod(ntrfl / 100, 10);
    //C inequality flag
    ntine = fem::mod(ntrfl / 1000, 10);
    //C limit flag
    ntlim = fem::mod(ntrfl / 10000, 10);
    //C profile flag
    ntprf = fem::mod(ntrfl / 100000, 10);
    //C 532 "aploop.F" 2
    if (ntvar == 0 || ntvar == 2 || ntvar == 3) {
      //C correct x and return to test constraints
      x(i) = xs(i) + dx(i);
      //C log-normal
    }
    else if (ntvar == 4) {
      x(i) = fem::exp(fem::log(xs(i)) + dx(i));
      //C sqrt
    }
    else if (ntvar == 5) {
      x(i) = fem::pow2((fem::sqrt(xs(i)) + dx(i)));
    }
  }
}

void
lesfcm(
  arr_ref<double> c)
{
  c(dimension(13));
  //C     ...
  //C determinant
  c(8) = c(1) * c(3) - c(2) * c(2);
  //C V_11
  c(10) = c(3) / c(8);
  //C V_12
  c(11) = -c(2) / c(8);
  //C V_22
  c(12) = c(1) / c(8);
  //C 1. parameter
  c(6) = c(4) * c(10) + c(5) * c(11);
  //C 2. parameter
  c(7) = c(4) * c(11) + c(5) * c(12);
}

//C
//C test convergence
void
antest(
  common& cmn,
  int& iret)
{
  // COMMON simcom
  double& epsf = cmn.epsf;
  double& chisq = cmn.chisq;
  double& ftest = cmn.ftest;
  double& ftestp = cmn.ftestp;
  double& chsqp = cmn.chsqp;
  double& frms = cmn.frms;
  double& frmsp = cmn.frmsp;
  double& weight = cmn.weight;
  int& iunph = cmn.iunph;
  int& ncst = cmn.ncst;
  int& iter = cmn.iter;
  //
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 545 "aploop.F" 2
  //C     ....
  //C calculate new Jacobian
  iret = -1;
  iunph = 0;
  //C     __________________________________________________________________
  //C     combined measure?
  float factor = fem::float0;
  int i = fem::int0;
  arr_1d<14, double> cm(fem::fill0);
  if (iter == 1 && ncst == 0) {
    factor = chisq / (fem::abs(ftestp - ftest) + epsf);
    FEM_DO_SAFE(i, 1, 14) {
      cm(i) = 0.0e0;
    }
    //C damping factor
    cm(9) = 0.5e0;
    cm(1) = 1.0e0 + 1.0e0;
    cm(2) = fem::pow2(frms) + fem::pow2(frmsp);
    cm(3) = fem::pow4(frms) + fem::pow4(frmsp);
    cm(4) = chisq + chsqp;
    cm(5) = chisq * fem::pow2(frms) + chsqp * fem::pow2(frmsp);
    //C fit
    lesfcm(cm);
    //C combined penalty
    cm(13) = cm(6) - cm(7) * fem::pow2(frms);
  }
  else if (ncst == 0) {
    cm(1) = cm(9) * cm(1) + 1.0e0;
    cm(2) = cm(9) * cm(2) + fem::pow2(frms);
    cm(3) = cm(9) * cm(3) + fem::pow4(frms);
    cm(4) = cm(9) * cm(4) + chisq;
    cm(5) = cm(9) * cm(5) + chisq * fem::pow2(frms);
    //C fit
    lesfcm(cm);
    cm(14) = cm(13);
    //C combined penalty
    cm(13) = cm(6) - cm(7) * fem::pow2(frms);
  }
  //C     __________________________________________________________________
  //C     cutstep
  if (ncst < 2 && (iunph != 0 || (iter > 1 && ftest > 2.0e0 * ftestp + epsf))) {
    ncst++;
    weight = 0.25e0;
    weight = 0.50e0;
    //C cutstep - add corrections
    iret = -2;
    return;
  }
  //C
  //C     __________________________________________________________________
  //C     convergent
  float dchisq = fem::float0;
  if (iter >= 2 && ncst == 0) {
    dchisq = chisq - chsqp;
    if (fem::abs(dchisq) <= cmn.epschi && ftest < epsf) {
      //C convergence
      iret = 0;
      return;
    }
  }
  //C     __________________________________________________________________
  //C     failure
  //C non-convergence
  if (iter > cmn.itermx) {
    iret = 2;
  }
}

void
acopxv(
  common& cmn,
  arr_cref<double> x,
  arr_ref<double> vx,
  arr_cref<double> dx,
  arr_ref<double> as,
  arr_ref<double> wm,
  arr_ref<double> pu)
{
  x(dimension(star));
  vx(dimension(star));
  dx(dimension(star));
  as(dimension(star));
  wm(dimension(star));
  pu(dimension(star));
  // COMMON simcom
  int& nx = cmn.nx;
  //
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 602 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 603 "aploop.F" 2
  //C     ...
  //C     __________________________________________________________________
  //C     convergence: pull calculation
  int ii = 0;
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, nx) {
    as(i) = x(i);
    ii += i;
    pu(i) = 0.0f;
    if (vx(ii) > 0.0f) {
      if (vx(ii) - wm(ii) > 0.0f) {
        pu(i) = dx(i) / fem::sqrt(vx(ii) - wm(ii));
      }
    }
  }
  //C     __________________________________________________________________
  //C     copy/exchange result/input covariance matrix
  double scopy = fem::double0;
  FEM_DO_SAFE(i, 1, (nx * nx + nx) / 2) {
    scopy = vx(i);
    //C copy fitted covariance matrix
    vx(i) = wm(i);
    //C ... and save input matrix
    wm(i) = scopy;
  }
}

//C
//C steering routine for loop
void
iploop(
  common& cmn,
  arr_ref<double> x,
  arr_ref<double> vx,
  arr_cref<double> f,
  arr_ref<double> xs,
  arr_ref<double> dx,
  arr_ref<double> fcopy,
  arr_ref<double> xp,
  arr_ref<double> rh,
  int& iret)
{
  x(dimension(star));
  vx(dimension(star));
  f(dimension(star));
  xs(dimension(star));
  dx(dimension(star));
  fcopy(dimension(star));
  xp(dimension(star));
  rh(dimension(star));
  double& ftest = cmn.ftest;
  double& frms = cmn.frms;
  int& nf = cmn.nf;
  int& indwm = cmn.indwm;
  int& iter = cmn.iter;
  int& ncalls = cmn.ncalls;
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  //
  int istatu = fem::int0;
  int nfit = fem::int0;
  int j = fem::int0;
  double fj = fem::double0;
  arr_1d<100, double> fex(fem::fill0);
  int jret = fem::int0;
  //C     ==================================================================
  //C     call IPLDER
  //C     call IPLCON
  //C     ==================================================================
  //C ,FOPT,FAC
  //C     local variables
  //C constraint values
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 111 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 112 "aploop.F" 2
  //C     __________________________________________________________________
  //C                    = -1  numerical derivatives
  //C                    =  0  constraint function evaluation
  //C                    =  1  constraint function with test afterwards
  //C                    =  2  end-of-fit
  //C     __________________________________________________________________
  //C     ...
  //C count calls
  ncalls++;
  //C
  //C     __________________________________________________________________
  //C     initialization
  if (ncalls != 1) {
    goto statement_20;
  }
  //C!!
  istatu = 0;
  //C reset fit count
  nfit = 0;
  iter = 0;
  FEM_DO_SAFE(j, 1, cmn.nx) {
    //C save initial X values
    xs(j) = x(j);
    //C reset correction DX
    dx(j) = 0.0e0;
  }
  //C     __________________________________________________________________
  //C     start/restart
  //C count fits
  nfit++;
  iter = 0;
  cmn.ncst = 0;
  cmn.chisq = 0.0e0;
  statement_20:
  FEM_DO_SAFE(j, 1, nf) {
    fj = f(j);
    //C extended F
    fex(j) = fj;
  }
  //C     __________________________________________________________________
  //C     constraint test summary
  if (istatu < 0) {
    goto statement_40;
  }
  //C save previous value
  cmn.ftestp = ftest;
  //C reset constraint tests
  ftest = 0.0e0;
  cmn.frmsp = frms;
  frms = 0.0e0;
  FEM_DO_SAFE(j, 1, nf) {
    fj = f(j);
    //C copy constraint vector
    fcopy(j) = fj;
    //C sum absolute values
    ftest += fem::abs(fj);
    //C sum squares
    frms += fem::pow2(fj);
  }
  //C average |F|
  ftest = fem::max(1.0e-16, ftest / fem::ffloat(nf));
  //C LS mean
  frms = fem::sqrt(frms / fem::ffloat(nf) + 1.0e-32);
  if (istatu == 1) {
    goto statement_60;
  }
  //C     __________________________________________________________________
  //C     start numerical derivatives
  //C!!
  statement_30:
  istatu = -1;
  //C     __________________________________________________________________
  //C     derivative calculation
  //C ISTATU=-1
  statement_40:
  if (istatu + 1 != 0) {
    goto statement_50;
  }
  //C derivative matrix A
  //C steps  ST(.)
  //C limits XL(2,.)
  //C copy FC(.) central F(.)
  //C copy HH(.) shifted F(.)
  anumde(cmn, x, fex, aux, aux(cmn.indst + 1), aux(cmn.indlm + 1),
    aux(1 + cmn.indfc), aux(1 + cmn.indhh), jret);
  iret = -1;
  //C...for constraint calculation
  if (jret < 0) {
    return;
  }
  //C
  //C     __________________________________________________________________
  //C     next iteration
  statement_50:
  aniter(cmn, x, vx, fcopy, aux, xp, rh, aux(indwm + 1), dx);
  goto statement_70;
  //C     __________________________________________________________________
  //C     test cutsteps
  statement_60:
  antest(cmn, iret);
  //C numerical derivative:   ISTATU=-1
  if (iret + 1 == 0) {
    goto statement_30;
  }
  //C convergence or failure: ISTATU= 2
  if (iret >= 0) {
    goto statement_80;
  }
  //C     __________________________________________________________________
  //C     apply corrections DX(.) to X(.) with transformations
  statement_70:
  addtox(cmn, x, xs, dx, xp);
  //C test at next entry
  istatu = 1;
  return;
  //C     __________________________________________________________________
  //C     end-of-primary-fit (NFIT=1)
  statement_80:
  istatu = 2;
  acopxv(cmn, x, vx, aux(cmn.inddx + 1), aux(cmn.indas + 1), aux(indwm + 1),
    aux(cmn.indpu + 1));
  iret = 0;
}

//C 1 "condutil.F"
//C 1 "<built-in>"
//C 1 "<command-line>"
//C 1 "condutil.F"
//C
//C index (I,J)=(J,I) in symmetric matri
int
ijsym(
  int const& i,
  int const& j)
{
  int return_value = fem::int0;
  //C,NCOUNT
  if (i <= j) {
    return_value = (j * j - j) / 2 + i;
  }
  else {
    return_value = (i * i - i) / 2 + j;
  }
  return return_value;
}

//C
//C define initial steps
void
asteps(
  common& cmn,
  arr_cref<double> x,
  arr_ref<double> vx,
  arr_ref<double> st)
{
  x(dimension(star));
  vx(dimension(star));
  st(dimension(star));
  common_write write(cmn);
  // COMMON simcom
  double& derfac = cmn.derfac;
  double& derlow = cmn.derlow;
  int& nx = cmn.nx;
  int& ipak = cmn.ipak;
  int& indtr = cmn.indtr;
  int& ndf = cmn.ndf;
  // COMMON nauxcm
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  //
  //C     ==================================================================
  //C     check transformed variables (e.g. > 0)
  //C     define initial step size for
  //C        o  measured variables
  //C        o  unmeasured variables
  //C     determine number of degrees of freedom
  //C     transform steps for transformed variables
  //C     transform covariance matrix for transformed variables
  //C     ==================================================================
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 212 "aploop.F" 2
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 213 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 214 "aploop.F" 2
  //C     ...
  int ii = 0;
  //C loop on all variables
  int i = fem::int0;
  int ntrfl = fem::int0;
  int ntvar = fem::int0;
  int ntmes = fem::int0;
  int ntder = fem::int0;
  int ntine = fem::int0;
  int ntlim = fem::int0;
  int ntprf = fem::int0;
  double vii = fem::double0;
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, nx) {
    ipak = i;
    //C
    //C 1 "unpackfl.inc" 1
    //C
    //C     unpackfl.inc = code for flag unpacking
    //C
    //C 60 "unpackfl.inc"
    //C get packed flags
    ntrfl = aux(indtr + ipak);
    //C transformation flag
    ntvar = fem::mod(ntrfl, 10);
    //C M-estimate flag
    ntmes = fem::mod(ntrfl / 10, 10);
    //C derivative type flag
    ntder = fem::mod(ntrfl / 100, 10);
    //C inequality flag
    ntine = fem::mod(ntrfl / 1000, 10);
    //C limit flag
    ntlim = fem::mod(ntrfl / 10000, 10);
    //C profile flag
    ntprf = fem::mod(ntrfl / 100000, 10);
    //C 220 "aploop.F" 2
    ii += i;
    //C original diagonal element
    vii = fem::abs(vx(ii));
    if (vii == 0.0e0) {
      //C set unmeasured flag
      ntmes = 1;
    }
    //C      _________________________________________________________________
    //C      check transformed variables
    //C logarithmic transformation - check
    if (ntvar == 4) {
      if (x(i) <= 0.0f) {
        //C reset: X(i) has to be positive
        ntvar = 0;
      }
      //C sqrt transformation - check
    }
    else if (ntvar == 5) {
      if (x(i) <= 0.0f) {
        //C reset: X(i) has to be positive
        ntvar = 0;
      }
    }
    //C      _________________________________________________________________
    //C      define step size for derivative calculation
    //C measured variable
    if (vii != 0.0e0) {
      if (st(i) > 0.0e0) {
        //C user step, if smaller
        st(i) = fem::min(st(i), derfac * fem::sqrt(vii));
        //C step is undefined
      }
      else if (st(i) <= 0.0e0) {
        //C step from cov matrix
        st(i) = derfac * fem::sqrt(vii);
        if (ntvar != 2) {
          st(i) = fem::min(st(i), derlow * fem::max(1.0e-6, fem::abs(x(i))));
          //C Poisson
        }
        else {
          st(i) = fem::min(st(i), derlow * fem::max(1.0e-6, fem::abs(
            x(i) + 1.0e0)));
        }
      }
      //C unmeasured variable
    }
    else if (vii == 0.0e0) {
      //C reduce degrees of freedom
      ndf = ndf - 1;
      FEM_DO_SAFE(j, 1, nx) {
        //C clear matrix elements
        vx(ijsym(i, j)) = 0.0e0;
      }
      st(i) = cmn.derufc * fem::max(1.0e0, fem::abs(x(i)));
    }
    //C fix
    if (ntine == 1) {
      st(i) = 0.0e0;
    }
    //C      _________________________________________________________________
    //C      transform steps for transformed variables
    //C
    if (st(i) == 0.0e0) {
      //C fixed by user
      ntine = 1;
      //C
    }
    else {
      //C lognormal variable
      if (ntvar == 3) {
        write(6, star), "Step from to ", st(i), st(i) / x(i);
        //C change step to log step
        st(i) = st(i) / x(i);
        //C sqrt variable
      }
      else if (ntvar == 4) {
        //C change step to sqrt step
        st(i) = 0.5e0 * st(i) / fem::sqrt(x(i));
      }
    }
    //C
    //C     __________________________________________________________________
    //C     transform covariance matrix for transformed variables
    //C transform covariance matrix for logn
    if (ntvar == 4) {
      write(6, star), "lognormal variable", i, x(i);
      FEM_DO_SAFE(j, 1, nx) {
        vx(ijsym(i, j)) = vx(ijsym(i, j)) / x(i);
        if (i == j) {
          vx(ijsym(i, j)) = vx(ijsym(i, j)) / x(i);
        }
      }
      //C ... and for sqrt
    }
    else if (ntvar == 5) {
      write(6, star), "sqrt variable", i, x(i);
      FEM_DO_SAFE(j, 1, nx) {
        vx(ijsym(i, j)) = vx(ijsym(i, j)) * 0.5e0 / fem::sqrt(x(i));
        if (i == j) {
          vx(ijsym(i, j)) = vx(ijsym(i, j)) * 0.5e0 / fem::sqrt(x(i));
        }
      }
    }
    //C
    //C 1 "packfl.inc" 1
    //C
    //C     packfl.inc   = code for flag packing
    //C
    //C     explanation see:
    //C     unpackfl.inc = code for flag unpacking
    //C
    aux(indtr + ipak) = ((((ntprf * 10 + ntlim) * 10 + ntine) * 10 +
      ntder) * 10 + ntmes) * 10 + ntvar;
    //C 286 "aploop.F" 2
  }
}

//C
//C steering routine for loop
void
aploop(
  common& cmn,
  arr_ref<double> x,
  arr_ref<double> vx,
  arr_cref<double> f,
  int& iret)
{
  x(dimension(star));
  vx(dimension(star));
  f(dimension(star));
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  int& nx = cmn.nx;
  int& indcf = cmn.indcf;
  int& indlm = cmn.indlm;
  int& indas = cmn.indas;
  int& indfc = cmn.indfc;
  int& indhh = cmn.indhh;
  int& indxs = cmn.indxs;
  int& inddx = cmn.inddx;
  int& indxp = cmn.indxp;
  int& indrh = cmn.indrh;
  int& indwm = cmn.indwm;
  int& india = cmn.india;
  int& ndtot = cmn.ndtot;
  int& nxf = cmn.nxf;
  int& mxf = cmn.mxf;
  int& indqn = cmn.indqn;
  int& indpu = cmn.indpu;
  int& ndtotl = cmn.ndtotl;
  //
  int nff = fem::int0;
  int i = fem::int0;
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 62 "aploop.F" 2
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 63 "aploop.F" 2
  //C
  if (cmn.ncalls != 0) {
    goto statement_10;
  }
  //C     __________________________________________________________________
  //C     indices/pointer etc at first APLOOP entry
  //C max. number of constraints
  nff = cmn.nf;
  //C total number of fit equations
  nxf = nx + nff;
  //C number elements symmetric matrix
  mxf = (nxf * nxf + nxf) / 2;
  //C pointer to FC(NF) = copy of F(NF)
  indfc = indlm + 2 * nx;
  //C pointer to FCOPY
  indcf = indfc + nff;
  //C pointer to HH(NF) = copy of F(.)
  indhh = indcf + nff;
  //C save X(.)        pointer
  indxs = indhh + nff;
  //C step
  inddx = indxs + nx;
  //C previos step
  indxp = inddx + nx;
  //C right-hand side
  indrh = indxp + nx;
  //C weight matrix
  indwm = indrh + nxf;
  //C matrix diagonal "DIAG"
  india = indwm + mxf;
  //C next pointer    "QNEXT"
  indqn = india + nxf;
  //C X result
  indas = indqn + nxf;
  //C pulls, solution X and Vx
  indpu = indas + nx;
  //C total number of words (so far)
  ndtot = indpu + nx;
  ndtotl = ndtot;
  if (ndtotl > naux) {
    FEM_STOP(0);
  }
  FEM_DO_SAFE(i, indlm + 2 * nx + 1, ndtot) {
    //C reset part of aux, unused so far
    aux(i) = 0.0e0;
  }
  //C initial steps ST(.)
  asteps(cmn, x, vx, aux(1 + cmn.indst));
  //C     __________________________________________________________________
  //C     internal APLOOP
  //C default status is -1 = continue
  statement_10:
  iret = -1;
  iploop(cmn, x, vx, f, aux(indxs + 1), aux(inddx + 1), aux(indcf + 1),
    aux(indxp + 1), aux(indrh + 1), iret);
}

struct dgamml_save
{
  arr<double> cof;
  double stp;

  dgamml_save() :
    cof(dimension(6), fem::fill0),
    stp(fem::double0)
  {}
};

//C
//C ln[Gamma(x)]
double
dgamml(
  common& cmn,
  double const& x)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(dgamml);
  // SAVE
  arr_ref<double> cof(sve.cof, dimension(6));
  double& stp = sve.stp;
  //
  if (is_called_first_time) {
    static const double values[] = {
      76.18009172947146e0, -86.50532032941677e0, 24.01409824083091e0,
        -1.231739572450155e0, 0.1208650973866179e-2,
        -.05395239384953e-5, 2.5066282746310005e0
    };
    fem::data_of_type<double>(FEM_VALUES_AND_SIZE),
      cof, stp;
  }
  //C     ...
  double xx = x;
  double yy = xx;
  double tmp = (xx + 0.5e0) * fem::log(xx + 5.5e0) - xx - 5.5e0;
  double ser = 1.000000000190015e0;
  int j = fem::int0;
  FEM_DO_SAFE(j, 1, 6) {
    yy += 1.0e0;
    ser += cof(j) / yy;
  }
  return_value = tmp + fem::log(stp * ser / xx);
  return return_value;
}

double
dgamin(
  common& cmn,
  double const& a,
  double const& x)
{
  double return_value = fem::double0;
  double gln = fem::double0;
  double ap = fem::double0;
  double sum = fem::double0;
  double del = fem::double0;
  int n = fem::int0;
  const int itmax = 300;
  const double eps = 3.0e-7f;
  double b = fem::double0;
  const double fpmin = 1.0e-30f;
  double c = fem::double0;
  double d = fem::double0;
  double h = fem::double0;
  int i = fem::int0;
  double an = fem::double0;
  //C     incomplete gamma function P(a,x)
  //C     returns -1.0 for x < 0 or a =< 0
  //C     ...
  if (x < 0.0e0 || a <= 0.0e0) {
    //C error
    return_value = -1.0e0;
  }
  else if (x == 0.0f) {
    return_value = 0.0e0;
    //C series representation
  }
  else if (x < a + 1.0e0) {
    //C ln[Gamma(a)]
    gln = dgamml(cmn, a);
    ap = a;
    sum = 1.0e0 / a;
    del = sum;
    FEM_DO_SAFE(n, 1, itmax) {
      ap += 1.f;
      del = del * x / ap;
      sum += del;
      if (fem::abs(del) < fem::abs(sum) * eps) {
        goto statement_10;
      }
    }
    FEM_STOP("DGAMIN:  error 1");
    statement_10:
    return_value = sum * fem::exp(-x + a * fem::log(x) - gln);
    //C continued fraction representation
  }
  else {
    //C ln[Gamma(a)]
    gln = dgamml(cmn, a);
    b = x + 1.0e0 - a;
    c = 1.0e0 / fpmin;
    d = 1.0e0 / b;
    h = d;
    FEM_DO_SAFE(i, 1, itmax) {
      an = -i * (i - a);
      b += 2.0e0;
      d = an * d + b;
      if (fem::abs(d) < fpmin) {
        d = fpmin;
      }
      c = b + an / c;
      if (fem::abs(c) < fpmin) {
        c = fpmin;
      }
      d = 1.0e0 / d;
      del = d * c;
      h = h * del;
      if (fem::abs(del - 1.0e0) < eps) {
        goto statement_20;
      }
    }
    return return_value;
    FEM_STOP("DGAMIN:  error 2");
    statement_20:
    return_value = 1.0e0 - fem::exp(-x + a * fem::log(x) - gln) * h;
  }
  return return_value;
}

//C 1 "chprob.F"
//C 1 "<built-in>"
//C 1 "<command-line>"
//C 1 "chprob.F"
//C
//C prob from chisquare
double
chprob(
  common& cmn,
  double const& chisq,
  int const& n)
{
  double return_value = fem::double0;
  //C     chi square probability for N degrees of freedom at CHISQ
  //C     =integral from 0 ... CHISQ
  //C     ...
  if (chisq <= 0.0f) {
    return_value = 1.0e0;
  }
  else {
    return_value = 1.0e0 - dgamin(cmn, 0.5e0 * n, 0.5e0 * chisq);
  }
  return return_value;
}

void
chndpv(
  common& cmn,
  float& chi2,
  int& nd,
  float& pval)
{
  // COMMON simcom
  double& chisq = cmn.chisq;
  int& ndf = cmn.ndf;
  //
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NX       = number of parameters
  //C     NF       = number of constraint equations
  //C     NUM      = flag for numerical differentiation
  //C                = 0  analytical derivatives
  //C                = 1  numerical derivatives
  //C                = 2  numerical derivatives, compare with analytical
  //C     IFLG     = flag for first case (check derivative)
  //C     INIT     = flag
  //C                = 1  deriving
  //C                = 0  else
  //C
  //C     EPSF     = required |F| accuracy
  //C     ISTAT    = status flag
  //C                =0
  //C                =1  derivative matrix finished
  //C                =2  X(.) corrections applied
  //C
  //C     XL(2,.)  = lower and upper values of parameters
  //C     ST(.)    = step sizes for numerical differentiation
  //C     FC(.)    = central values of parameters
  //C     H(.)     = copy
  //C
  //C     A        = derivative matrix a/flags during matrix inversion/pulls
  //C     A(NX,NF)
  //C
  //C     ******************************************************************
  //C
  //C     ND       = number of degrees of freedom
  //C              = number of constraint equations minus number of
  //C                unmeasured parameters
  //C     CHSQ     = chi square
  //C
  //C     DERFAC   = factor for standard deviation in numerical derivatives
  //C
  //C     NDENDE   = index of last used word incl. single-precision array
  //C     NDACTL   = index of actual single-precision array
  //C
  //C     =========================end=of=macro=============================
  //C 637 "aploop.F" 2
  //C     ...
  //C chi^square
  chi2 = chisq;
  //C number of degrees of freedom
  nd = ndf;
  //C p-value
  pval = chprob(cmn, chisq, ndf);
}

//C
//C inverse normal distribution
double
dingau(
  double const& p)
{
  double return_value = fem::double0;
  double a1 = fem::double0;
  double a2 = fem::double0;
  double a3 = fem::double0;
  double a4 = fem::double0;
  double a5 = fem::double0;
  double a6 = fem::double0;
  double b1 = fem::double0;
  double b2 = fem::double0;
  double b3 = fem::double0;
  double b4 = fem::double0;
  double b5 = fem::double0;
  double c1 = fem::double0;
  double c2 = fem::double0;
  double c3 = fem::double0;
  double c4 = fem::double0;
  double c5 = fem::double0;
  double c6 = fem::double0;
  double d1 = fem::double0;
  double d2 = fem::double0;
  double d3 = fem::double0;
  double d4 = fem::double0;
  double p_low = fem::double0;
  double p_high = fem::double0;
  double q = fem::double0;
  double z = fem::double0;
  double r = fem::double0;
  //C     routine written by john herrero
  //C     ...
  a1 = -39.6968302866538f;
  a2 = 220.946098424521f;
  a3 = -275.928510446969f;
  a4 = 138.357751867269f;
  a5 = -30.6647980661472f;
  a6 = 2.50662827745924f;
  b1 = -54.4760987982241f;
  b2 = 161.585836858041f;
  b3 = -155.698979859887f;
  b4 = 66.8013118877197f;
  b5 = -13.2806815528857f;
  c1 = -0.00778489400243029f;
  c2 = -0.322396458041136f;
  c3 = -2.40075827716184f;
  c4 = -2.54973253934373f;
  c5 = 4.37466414146497f;
  c6 = 2.93816398269878f;
  d1 = 0.00778469570904146f;
  d2 = 0.32246712907004f;
  d3 = 2.445134137143f;
  d4 = 3.75440866190742f;
  p_low = 0.02425f;
  p_high = 1 - p_low;
  if (p < p_low) {
    goto statement_201;
  }
  if (p >= p_low) {
    goto statement_301;
  }
  statement_201:
  q = fem::dsqrt(-2 * fem::dlog(p));
  z = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((
    d1 * q + d2) * q + d3) * q + d4) * q + 1);
  goto statement_204;
  statement_301:
  if ((p >= p_low) && (p <= p_high)) {
    goto statement_202;
  }
  if (p > p_high) {
    goto statement_302;
  }
  statement_202:
  q = p - 0.5f;
  r = q * q;
  z = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
    (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
  goto statement_204;
  statement_302:
  if ((p > p_high) && (p < 1)) {
    goto statement_203;
  }
  statement_203:
  q = fem::dsqrt(-2 * fem::dlog(1 - p));
  z = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / (((
    (d1 * q + d2) * q + d3) * q + d4) * q + 1);
  statement_204:
  return_value = z;
  return return_value;
}

} // namespace aplcon
