#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

void
a1prof(...)
{
  throw std::runtime_error(
    "Missing function implementation: a1prof");
}

void
a2prof(...)
{
  throw std::runtime_error(
    "Missing function implementation: a2prof");
}

void
aprini(...)
{
  throw std::runtime_error(
    "Missing function implementation: aprini");
}

void
atitoe(...)
{
  throw std::runtime_error(
    "Missing function implementation: atitoe");
}

void
b1prof(...)
{
  throw std::runtime_error(
    "Missing function implementation: b1prof");
}

void
b2prof(...)
{
  throw std::runtime_error(
    "Missing function implementation: b2prof");
}

void
cfcorr(...)
{
  throw std::runtime_error(
    "Missing function implementation: cfcorr");
}

void
cfgmpr(...)
{
  throw std::runtime_error(
    "Missing function implementation: cfgmpr");
}

void
cfprv(...)
{
  throw std::runtime_error(
    "Missing function implementation: cfprv");
}

void
cfprvp(...)
{
  throw std::runtime_error(
    "Missing function implementation: cfprvp");
}

double
chprob(...)
{
  throw std::runtime_error(
    "Missing function implementation: chprob");
}

void
ciprv(...)
{
  throw std::runtime_error(
    "Missing function implementation: ciprv");
}

double
dingau(...)
{
  throw std::runtime_error(
    "Missing function implementation: dingau");
}

void
duminv(...)
{
  throw std::runtime_error(
    "Missing function implementation: duminv");
}

int
ijsym(...)
{
  throw std::runtime_error(
    "Missing function implementation: ijsym");
}

double
scalxy(...)
{
  throw std::runtime_error(
    "Missing function implementation: scalxy");
}

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
  double penalt;
  int nadfs;
  int nx;
  int nf;
  int num;
  int iflg;
  int init;
  int lunsim;
  int ipr;
  int ncase;
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
    penalt(fem::double0),
    nadfs(fem::int0),
    nx(fem::int0),
    nf(fem::int0),
    num(fem::int0),
    iflg(fem::int0),
    init(fem::int0),
    lunsim(fem::int0),
    ipr(fem::int0),
    ncase(fem::int0),
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

struct common_cprofl
{
  static const int mlr = 36;
  static const int mmlr = mlr + mlr + 1;
  static const int mseca = 100;

  double xl;
  double xr;
  arr<double> clr;
  arr<double> xlr;
  arr<double> ylr;
  arr<double> flr;
  double constr;
  double constx;
  double consty;
  double chiref;
  arr<double, 2> csp;
  double center;
  double sigmax;
  double sigmay;
  double centex;
  double centey;
  double red;
  double redl;
  double redr;
  double cphi;
  double sphi;
  double crt;
  double srt;
  double cth;
  double sth;
  double cht;
  arr<float, 2> xcont;
  arr<float, 2> ycont;
  arr<float> xc;
  arr<float> yc;
  arr<int, 2> npsec;
  int ilr;
  int nlr;
  int iseca;
  int nseca;
  int nfadd;
  int ipf;
  int ndext;
  int ipfx;
  int ipfy;
  int nstar;
  int istar;

  common_cprofl() :
    xl(fem::double0),
    xr(fem::double0),
    clr(dim1(-mlr,  + mlr), fem::fill0),
    xlr(dim1(-mlr,  + mlr), fem::fill0),
    ylr(dim1(-mlr,  + mlr), fem::fill0),
    flr(dim1(-mlr,  + mlr), fem::fill0),
    constr(fem::double0),
    constx(fem::double0),
    consty(fem::double0),
    chiref(fem::double0),
    csp(dimension(5, mmlr), fem::fill0),
    center(fem::double0),
    sigmax(fem::double0),
    sigmay(fem::double0),
    centex(fem::double0),
    centey(fem::double0),
    red(fem::double0),
    redl(fem::double0),
    redr(fem::double0),
    cphi(fem::double0),
    sphi(fem::double0),
    crt(fem::double0),
    srt(fem::double0),
    cth(fem::double0),
    sth(fem::double0),
    cht(fem::double0),
    xcont(dimension(24, 3), fem::fill0),
    ycont(dimension(24, 3), fem::fill0),
    xc(dimension(101), fem::fill0),
    yc(dimension(101), fem::fill0),
    npsec(dimension(2, mseca), fem::fill0),
    ilr(fem::int0),
    nlr(fem::int0),
    iseca(fem::int0),
    nseca(fem::int0),
    nfadd(fem::int0),
    ipf(fem::int0),
    ndext(fem::int0),
    ipfx(fem::int0),
    ipfy(fem::int0),
    nstar(fem::int0),
    istar(fem::int0)
  {}
};

const int common_cprofl::mlr;
const int common_cprofl::mmlr;
const int common_cprofl::mseca;

struct common :
  fem::common,
  common_simcom,
  common_nauxcm,
  common_cprofl
{
  fem::cmn_sve aplcon_sve;
  fem::cmn_sve anumde_sve;
  fem::cmn_sve aiprin_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct aplcon_save
{
  bool start;

  aplcon_save() :
    start(fem::bool0)
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
  FEM_CMN_SVE(aplcon);
  // COMMON simcom
  int& nx = cmn.nx;
  int& nf = cmn.nf;
  int& ncase = cmn.ncase;
  int& indst = cmn.indst;
  int& indlm = cmn.indlm;
  int& indtr = cmn.indtr;
  int& ndtot = cmn.ndtot;
  int& nxf = cmn.nxf;
  int& itermx = cmn.itermx;
  // COMMON nauxcm
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  //
  // SAVE
  bool& start = sve.start;
  //
  if (is_called_first_time) {
    start = true;
  }
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 1 "cprofil.inc" 1
  //C
  //C     Common for profile likelihood analysis with APLCON
  //C
  //C     NDEXT +1     is start of save area for primary X(.) and VX(.)
  //C     ILR    =     index of current profile point
  //C     NLR    =     number of profile points
  //C     ISECA  =     index of current profile analysis
  //C     NSECA  =     number of profile analyses
  //C     NFADD  =     max number of additional constraints
  //C     IPF    =     index of profile variable
  //C     NPSEC(2,.)   indices for profile analyses
  //C     XL,XR  =     limits for profile analysis
  //C     XLR(.) =     X values for profile
  //C     YLR(.) =     Y values for profile
  //C     FLR(.) =     chi^2 - reference value for profile
  //C     CONSTR =     fixed value
  //C     CONSTX,Y     fixed values
  //C     CENTER =     fitted value (minimum)
  //C     SIGMAX =     parabolic X error
  //C     SIGMAY =     parabolic Y error
  //C     CENTEX =     X center
  //C     CENTEY =     Y center
  //C     RED    =     reduction factor, e.g. 1/2
  //C     REDL   =     reduction factor left
  //C     REDR   =     reduction factor right
  //C     CPHI   =     cos phi of eigenvector
  //C     SPHI   =     sin phi
  //C     NSTAR  =     number of directions
  //C     ISTAR  =     direction index
  //C     CRT    =     cos of rotation transformation
  //C     SRT    =     sin of rotation transformation
  //C     CTH    =     cos rotation
  //C     STH    =     sin rotation
  //C     CHT    =     copy of cos rotation
  //C     XCONT(24,3)  x contur data, max 12 points, 3 contours
  //C     YCONT(24,3)  y contur data, max 12 points, 3 contours
  //C
  //C 14 "aploop.F" 2
  //C     ...
  if (start) {
    start = false;
    //C reset counter
    ncase = 0;
    cmn.ipr = 5;
  }
  //C count cases
  ncase++;
  //C
  //C reset number of profile searches
  cmn.nseca = 0;
  //C reset number of additional constraints
  cmn.nfadd = 0;
  //C primary fit initialization
  cmn.nadfs = 0;
  //C using space for NF+2 instead of NF
  //C
  //C number of variables
  nx = nvar;
  //C number of constraint equations
  nf = mcst;
  //C primary value of NF
  cmn.nfprim = nf;
  //C dimension of AUX array
  cmn.ndpdim = naux;
  //C printout unit (default)
  cmn.lunsim = 6;
  //C     IPR=5               ! default print flag is 5
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
  itermx = 100;
  itermx = 50;
  itermx = 20;
  itermx = 10;
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
    //C error print
    aprini(1);
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
  //C initial print without X, VX
  aprini(0);
  //C      CALL APNAME(0,' ')  ! reset parameter names
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
  common_write write(cmn);
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 494 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 495 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 496 "aploop.F" 2
  //C entry flag
  //C     ...
  //C      WRITE(*,*) 'NUMDE entered',TINUE
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
  //C      WRITE(*,*) 'I ST(I)',I,ST(I)
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
  //C 517 "aploop.F" 2
  //C      WRITE(*,*) 'I NTDER',I,NTDER
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
          if (cmn.ipr > 1) {
            write(cmn.lunsim, star), "Variable", i, " is", x(i),
              " reset to normal";
          }
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
          //C 569 "aploop.F" 2
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
      //C         DER=0.5D0*(HH(J)-F(J))/ST(I) ! numerical 1. derivative
      //C!! internal variable
      der = (hh(j) - f(j)) / (xt(1) - xt(2));
      //C          WRITE(*,246) I,J,HH(J),F(J),XT(1),XT(2),DER
      //C 246      FORMAT('der I J',2I3,5F10.4)
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
  //C 646 "aploop.F" 2
  goto statement_10;
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 653 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 654 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 655 "aploop.F" 2
  //C     ...
  //C start next iteration
  iter++;
  //C save current chi^2
  chsqp = chisq;
  //C
  //C      A(2)=1.0D0
  //C      A(3)=-1.0D0
  //C      WRITE(*,321) 'Jacobian',(A(I),I=1,NX*NF)
  //C321  FORMAT(A/(4F15.6))
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
  //C      WRITE(*,*) 'RHS',(RH(I),I=1,NX+NF)
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
    //C 692 "aploop.F" 2
    //C       WRITE(*,*) 'NTVAR=2 II',II,X(I),NTVAR
    //C Poisson
    if (ntvar == 2) {
      //C -MAX(ABS(X(I)),1.0D0)
      wm(ii) = -fem::sqrt(1.0f + fem::pow2(x(i)));
    }
  }
  int nrank = fem::int0;
  arr<double> diag(dimension(1000), fem::fill0);
  arr<double> qnext(dimension(1000), fem::fill0);
  duminv(a, wm, rh, nx, nf, 1, nrank, diag, qnext);
  //C      WRITE(*,*) NRANK,' is rank. NX NF ',NX,NF
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
  //C      WRITE(*,*) 'Prev corrections',(XP(I),I=1,NX)
  //C      WRITE(*,*) 'Curr corrections',(DX(I),I=1,NX)
  //C      WRITE(*,*) 'Weight=',WEIGHT
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 720 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 721 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 722 "aploop.F" 2
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
    //C 727 "aploop.F" 2
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

struct aiprin_save
{
  arr<fem::str<19> > text;

  aiprin_save() :
    text(dimension(5), fem::fill0)
  {}
};

void
aiprin(
  common& cmn,
  arr_cref<double> x,
  arr_cref<double> vx,
  int const& iarg,
  int const& iret)
{
  FEM_CMN_SVE(aiprin);
  x(dimension(star));
  vx(dimension(star));
  common_write write(cmn);
  double& chisq = cmn.chisq;
  double& ftest = cmn.ftest;
  double& frms = cmn.frms;
  int& nx = cmn.nx;
  int& nf = cmn.nf;
  int& lunsim = cmn.lunsim;
  int& ipr = cmn.ipr;
  int& ncase = cmn.ncase;
  int& ndf = cmn.ndf;
  int& ncst = cmn.ncst;
  int& iter = cmn.iter;
  int& ncalls = cmn.ncalls;
  int& ndtotl = cmn.ndtotl;
  const int naux = 125000;
  arr_cref<double> aux(cmn.aux, dimension(naux));
  //
  str_arr_ref<1> text(sve.text, dimension(5));
  if (is_called_first_time) {
    static const char* values[] = {
      "Chisquare too high ", "Too many iterations",
        "Unphysical region  ", "NDF less or equal 0",
        "AUX dimension small"
    };
    fem::data_of_type_str(FEM_VALUES_AND_SIZE),
      text;
  }
  float f = fem::float0;
  double chp = fem::double0;
  int i = fem::int0;
  int j = fem::int0;
  static const char* format_101 = "(3x,a30,i7)";
  static const char* format_102 = "(3x,a30,5x,1p,3e8.1)";
  //C      INTEGER NDTOTL
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 908 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 909 "aploop.F" 2
  //C
  //C 1 "cprofil.inc" 1
  //C
  //C     Common for profile likelihood analysis with APLCON
  //C
  //C     NDEXT +1     is start of save area for primary X(.) and VX(.)
  //C     ILR    =     index of current profile point
  //C     NLR    =     number of profile points
  //C     ISECA  =     index of current profile analysis
  //C     NSECA  =     number of profile analyses
  //C     NFADD  =     max number of additional constraints
  //C     IPF    =     index of profile variable
  //C     NPSEC(2,.)   indices for profile analyses
  //C     XL,XR  =     limits for profile analysis
  //C     XLR(.) =     X values for profile
  //C     YLR(.) =     Y values for profile
  //C     FLR(.) =     chi^2 - reference value for profile
  //C     CONSTR =     fixed value
  //C     CONSTX,Y     fixed values
  //C     CENTER =     fitted value (minimum)
  //C     SIGMAX =     parabolic X error
  //C     SIGMAY =     parabolic Y error
  //C     CENTEX =     X center
  //C     CENTEY =     Y center
  //C     RED    =     reduction factor, e.g. 1/2
  //C     REDL   =     reduction factor left
  //C     REDR   =     reduction factor right
  //C     CPHI   =     cos phi of eigenvector
  //C     SPHI   =     sin phi
  //C     NSTAR  =     number of directions
  //C     ISTAR  =     direction index
  //C     CRT    =     cos of rotation transformation
  //C     SRT    =     sin of rotation transformation
  //C     CTH    =     cos rotation
  //C     STH    =     sin rotation
  //C     CHT    =     copy of cos rotation
  //C     XCONT(24,3)  x contur data, max 12 points, 3 contours
  //C     YCONT(24,3)  y contur data, max 12 points, 3 contours
  //C
  //C 910 "aploop.F" 2
  //C     ...
  //C      NDTOTL=120
  if (iarg == 1) {
    goto statement_10;
  }
  if (iarg == 2) {
    goto statement_20;
  }
  if (iarg == 3) {
    goto statement_30;
  }
  if (iarg == 4) {
    goto statement_40;
  }
  if (iarg == 5) {
    goto statement_50;
  }
  //C     __________________________________________________________________
  //C     APLCON Titel print
  //C print
  if (ipr >= 3) {
    write(lunsim, star), " ";
    write(lunsim, star), "APLCON - constrained least squares",
      "              Version  01/10/2011";
    write(lunsim, star), " ";
    if (ipr >= 5 && nx <= 128) {
      cfprv(lunsim, x, vx, nx);
    }
  }
  return;
  //C     __________________________________________________________________
  //C     APLCON memory space error
  statement_10:
  if (ndtotl > naux) {
    write(6, star), " ";
    write(6, star), "APLCON  - constrained least squares";
    write(6, star), "          Case ", ncase;
    write(6, star), "Insufficient space in internal array AUX(", naux, ")";
    write(6, star), "Required:", ndtotl, " elements (at least)";
    write(6, star), "-> stop";
    write(6, star), " ";
    FEM_STOP(0);
  }
  return;
  //C     __________________________________________________________________
  //C     initial printout normal start
  statement_20:
  if (ipr >= 3) {
    write(lunsim, star), " ";
    write(lunsim, format_101), " Constrained least squares fit:";
    write(lunsim, format_101), "                         case", ncase;
    write(lunsim, format_101), "              nr of variables", nx;
    write(lunsim, format_101), "            nr of constraints", nf;
    write(lunsim, format_101), "           degrees of freedom", ndf;
    write(lunsim, format_102), "                      epsilon", cmn.epsf;
    write(lunsim, format_102), "     factors for numer. diff.",
      cmn.derfac, cmn.derufc, cmn.derlow;
    write(lunsim, format_101), "          used array elements", cmn.ndtot;
    write(lunsim, format_101), "         total array elements", cmn.ndpdim;
  }
  //C
  //C print
  if (ipr >= 6 && nx <= 32) {
    //C initial parameter values
    cfprv(lunsim, x, vx, nx);
  }
  return;
  //C     __________________________________________________________________
  //C     iteration printout
  statement_30:
  if (ipr >= 4 && iter == 0) {
    write(lunsim,
      "(/,'Iteration    calls        chi^2  Frms       |F| ',6x,'dchi^2',2x,"
      "'chi^2&Frms')");
    write(lunsim, "(1x,i5,2x,i10,13x,e10.2,e10.2)"), iter, ncalls, frms, ftest;
  }
  if (ipr >= 4 && iter >= 1) {
    if (ncst == 0) {
      write(lunsim, "(1x,i5,2x,i10,f13.3,e10.2,e10.2,e11.2,f10.3)"),
        iter, ncalls, chisq, frms, ftest, chisq - cmn.chsqp,
        cmn.penalt;
    }
    else {
      write(lunsim, "(1x,i5,'.',i1,i10,13x,e10.2,e10.2)"), iter,
        ncst, ncalls, frms, ftest;
    }
  }
  if (ipr >= 7 && ncst == 0) {
    cfgmpr(lunsim, x, nx, 1, "of X-values");
    cfgmpr(lunsim, f, nf, 1, " of constraint function values");
  }
  return;
  //C     __________________________________________________________________
  //C     final return with IRET = 0,1,2,3,4,5
  //C print
  statement_40:
  if (ipr >= 1) {
    //C         WRITE(*,*) 'Printing result IPR, IRET=',IPR,IRET
    //C         IF(IRET.GT.5) IRET=0 !!!!!!!!!!!!!!!
    //C convergence (IRET=0)
    if (iret == 0) {
      //C            WRITE(*,*) 'Printing result IPR, IRET=',IPR,IRET
      if (ipr >= 3) {
        //C               WRITE(*,*) 'before print ITER ... IPR=',IPR
        write(lunsim, star), " ";
        write(lunsim,
          "('Convergence after',i7,' iterations with chi^2=',g15.8,'  ndf=',"
          "i4)"),
          iter, chisq, ndf;
        if (ndf > 0) {
          //C p-value
          chp = chprob(chisq, ndf);
          write(lunsim, "(22x,'corresponding to P-value=',f11.6)"), chp;
          if (chp < 0.05e0) {
            write(lunsim,
              "('P-value converted to number of standard deviations:')");
            write(lunsim,
              "(16x,'1-,2-sided standard deviations=',f11.6,',',f13.6)"),
              dingau(1.0e0 - 0.5e0 * chp), dingau(1.0e0 - chp);
          }
        }
      }
      if (ipr >= 4) {
        if (nx <= 128) {
          cfprvp(lunsim, x, vx, aux(cmn.indpu + 1), nx);
          if (ipr >= 5) {
            write(lunsim, star), " ";
            cfcorr(lunsim, vx, nx);
          }
        }
        else {
          write(lunsim, star), "  x-vector (fitted):";
          {
            write_loop wloop(cmn, lunsim, "(3x,5g12.5)");
            FEM_DO_SAFE(i, 1, nx) {
              wloop, x(i);
            }
          }
        }
        write(lunsim, star), " ";
      }
      //C non-convergence (IRET > 0)
    }
    else if (ipr >= 3) {
      write(lunsim,
        "(' No convergence (',i3,' ITER, IRET =',i2,')  ',a)"), iter,
        iret, text(iret);
      write(lunsim, star), " ";
    }
  }
  return;
  //C     __________________________________________________________________
  //C
  statement_50:
  return;
  //C     __________________________________________________________________
  //C
  // UNHANDLED: ENTRY aprini(iarg)
  //C print
  if (ipr >= 3) {
    {
      write_loop wloop(cmn, lunsim, star);
      FEM_DO_SAFE(j, 1, 71) {
        wloop, "_";
      }
    }
    write(lunsim, star), " ";
    write(lunsim, star), "APLCON - constrained least squares",
      "                  Version  01/10/2011";
    {
      write_loop wloop(cmn, lunsim, star);
      FEM_DO_SAFE(j, 1, 71) {
        wloop, "_";
      }
    }
    write(lunsim, star), " ";
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
  //C      WRITE(*,*) 'Solution ',C(6),C(7)
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 740 "aploop.F" 2
  //C     ....
  //C calculate new Jacobian
  iret = -1;
  iunph = 0;
  //C     __________________________________________________________________
  //C     combined measure?
  //C      WRITE(*,*) 'ITER,NCST',ITER,NCST,CHISQ,FRMS,CHSQP,FRMSP
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
    //C         WRITE(*,*) 'Penalty',CM(6),CM(7),CM(13)
  }
  //C      WRITE(*,*) 'ANTEST:',CHISQ,FTEST,CHISQ+FACTOR*FTEST
  cmn.penalt = cm(13);
  float x = fem::float0;
  float vx = fem::float0;
  aiprin(cmn, x, vx, 3, iret);
  //C     __________________________________________________________________
  //C     cutstep
  //C      WRITE(*,*) 'Test cutstep:',FTEST,FTESTP,EPSF,2.0D0*FTESTP+EPSF
  if (ncst < 2 && (iunph != 0 || (iter > 1 && ftest > 2.0e0 * ftestp + epsf))) {
    ncst++;
    weight = 0.25e0;
    weight = 0.50e0;
    //C         IF(FTEST/FTESTP.GT. 5.0D0) WEIGHT=0.10D0
    //C         IF(FTEST/FTESTP.GT.10.0D0) WEIGHT=0.05D0
    //C cutstep - add corrections
    iret = -2;
    return;
  }
  //C
  //C      IF(NCST.EQ.0.AND.PENALT.LT.0.0D0) THEN
  //C         WEIGHT=0.1D0
  //C         NCST=NCST+1
  //C         IRET=-2                  ! cutstep - add corrections
  //C         WRITE(*,*) 'Penalty!'
  //C         RETURN
  //C      END IF
  //C
  //C      IF(NCST.EQ.0.AND.CM(13).GT.CM(14)) THEN
  //C         WEIGHT=0.1D0
  //C         NCST=NCST+1
  //C         IRET=-2                  ! cutstep - add corrections
  //C         WRITE(*,*) 'Penalty increasing!'
  //C         RETURN
  //C      END IF
  //C
  //C      IF(ITER.EQ.1.AND.NCST.EQ.0) THEN
  //C         WEIGHT=0.1
  //C         NCST=NCST+1
  //C         IRET=-2                  ! cutstep - add corrections
  //C         RETURN
  //C      END IF
  //C
  //C     __________________________________________________________________
  //C     convergent
  float dchisq = fem::float0;
  if (iter >= 2 && ncst == 0) {
    dchisq = chisq - chsqp;
    if (fem::abs(dchisq) <= cmn.epschi && ftest < epsf) {
      //C convergence
      iret = 0;
      //C
      //C            WRITE(*,*) '>>>> convergence'
      //C
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 831 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 832 "aploop.F" 2
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
  //C     __________________________________________________________________
  //C     apply transformation to covariance matrix
  atitoe(x, vx);
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
  common_write write(cmn);
  double& chisq = cmn.chisq;
  double& ftest = cmn.ftest;
  double& frms = cmn.frms;
  int& nadfs = cmn.nadfs;
  int& nx = cmn.nx;
  int& nf = cmn.nf;
  int& ipr = cmn.ipr;
  int& ndende = cmn.ndende;
  int& indas = cmn.indas;
  int& indvs = cmn.indvs;
  int& indwm = cmn.indwm;
  int& ndtot = cmn.ndtot;
  int& iter = cmn.iter;
  int& ncalls = cmn.ncalls;
  int& indpu = cmn.indpu;
  int& nfprim = cmn.nfprim;
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  double& chiref = cmn.chiref;
  const int mseca = 100;
  arr_cref<int, 2> npsec(cmn.npsec, dimension(2, mseca));
  int& iseca = cmn.iseca;
  int& nseca = cmn.nseca;
  int& ipf = cmn.ipf;
  int& ipfx = cmn.ipfx;
  int& ipfy = cmn.ipfy;
  //
  int istatu = fem::int0;
  int nfit = fem::int0;
  int j = fem::int0;
  arr_1d<2, double> fadd(fem::fill0);
  double fj = fem::double0;
  arr_1d<100, double> fex(fem::fill0);
  int jret = fem::int0;
  int iprsav = fem::int0;
  int kret = fem::int0;
  //C     ==================================================================
  //C     initial print
  //C        add 1/2 constraints for profile analysis
  //C     call IPLDER
  //C     call IPLCON
  //C        profile analysis
  //C     ==================================================================
  //C      INTEGER NITER,NFIT
  //C      INTEGER J,IRET,JRET,NSECAS,IJSYM,ILRP,ILR1,ILR2,NFUN ,NN,NTLIMP
  //C ,FOPT,FAC
  //C     local variables
  //C constraint values
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 149 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 150 "aploop.F" 2
  //C
  //C 1 "cprofil.inc" 1
  //C
  //C     Common for profile likelihood analysis with APLCON
  //C
  //C     NDEXT +1     is start of save area for primary X(.) and VX(.)
  //C     ILR    =     index of current profile point
  //C     NLR    =     number of profile points
  //C     ISECA  =     index of current profile analysis
  //C     NSECA  =     number of profile analyses
  //C     NFADD  =     max number of additional constraints
  //C     IPF    =     index of profile variable
  //C     NPSEC(2,.)   indices for profile analyses
  //C     XL,XR  =     limits for profile analysis
  //C     XLR(.) =     X values for profile
  //C     YLR(.) =     Y values for profile
  //C     FLR(.) =     chi^2 - reference value for profile
  //C     CONSTR =     fixed value
  //C     CONSTX,Y     fixed values
  //C     CENTER =     fitted value (minimum)
  //C     SIGMAX =     parabolic X error
  //C     SIGMAY =     parabolic Y error
  //C     CENTEX =     X center
  //C     CENTEY =     Y center
  //C     RED    =     reduction factor, e.g. 1/2
  //C     REDL   =     reduction factor left
  //C     REDR   =     reduction factor right
  //C     CPHI   =     cos phi of eigenvector
  //C     SPHI   =     sin phi
  //C     NSTAR  =     number of directions
  //C     ISTAR  =     direction index
  //C     CRT    =     cos of rotation transformation
  //C     SRT    =     sin of rotation transformation
  //C     CTH    =     cos rotation
  //C     STH    =     sin rotation
  //C     CHT    =     copy of cos rotation
  //C     XCONT(24,3)  x contur data, max 12 points, 3 contours
  //C     YCONT(24,3)  y contur data, max 12 points, 3 contours
  //C
  //C 151 "aploop.F" 2
  //C#include "declarefl.inc"
  //C      REAL XRP,YRP
  //C      DOUBLE PRECISION DCSPN
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
  //C no additional constraints
  nadfs = 0;
  //C index of profile analysis
  iseca = 0;
  //C initial printout
  if (ipr >= 3) {
    aiprin(cmn, x, vx, 1, iret);
  }
  FEM_DO_SAFE(j, 1, nx) {
    //C save initial X values
    xs(j) = x(j);
    //C reset correction DX
    dx(j) = 0.0e0;
  }
  //C     __________________________________________________________________
  //C     start/restart
  //C count fits
  statement_10:
  nfit++;
  //C      IF(NFIT.GE.10) STOP
  iter = 0;
  cmn.ncst = 0;
  chisq = 0.0e0;
  //C      WRITE(*,*) 'Start/restart fit NFIT ISTATU ',NFIT,ISTATU
  //C      WRITE(*,*) 'NF NX NADFS ',NF,NX,NADFS
  //C      WRITE(*,321) 'Vx',(VX(J),J=1,6)
  //C      WRITE(*,321) 'X ',(X(J),J=1,3)
  //C      WRITE(*,321) 'XS',(XS(J),J=1,3)
  //C      WRITE(*,321) 'DX',(DX(J),J=1,3)
  //C      WRITE(*,*)   'CONSTR IPF',CONSTR,IPF
  //C      WRITE(*,321) 'F',(F(J),J=1,NF)
  //C      WRITE(*,*) 'INDST ',INDST
  //C      WRITE(*,321) 'ST',(AUX(INDST+J),J=1,3)
  //C     __________________________________________________________________
  //C     add constraints for profile analysis
  //C add constraints for profiles
  statement_20:
  if (nadfs == 1) {
    //C 1-parameter constraint
    fadd(1) = cmn.constr - x(ipf);
    //C         WRITE(*,*) 'Extra constraint',FADD(1),CONSTR,X(IPF),IPF
  }
  else if (nadfs == 2) {
    //C 2-parameter constraint
    fadd(1) = cmn.constx - x(ipfx);
    fadd(2) = cmn.consty - x(ipfy);
  }
  FEM_DO_SAFE(j, 1, nf) {
    if (nadfs == 1 && j == nf) {
      fj = fadd(1);
    }
    else if (nadfs == 2 && j + 1 >= nf) {
      fj = fadd(j - nf + 2);
    }
    else {
      fj = f(j);
    }
    //C extended F
    fex(j) = fj;
  }
  //C     __________________________________________________________________
  //C     constraint test summary
  //C      IF(NFIT.GT.1) WRITE(*,*) 'ISTATU=',ISTATU
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
    if (nadfs == 1 && j == nf) {
      fj = fadd(1);
    }
    else if (nadfs == 2 && j + 1 >= nf) {
      fj = fadd(j - nf + 2);
    }
    else {
      fj = f(j);
    }
    //C copy constraint vector
    fcopy(j) = fj;
    //C sum absolute values
    ftest += fem::abs(fj);
    //C sum squares
    frms += fem::pow2(fj);
  }
  //C      WRITE(*,*) 'ISTATU=',ISTATU
  //C      WRITE(*,321) 'X    ',(X    (J),J=1,NX)
  //C      WRITE(*,321) 'FCOPY',(FCOPY(J),J=1,NF)
  //C321  FORMAT(A/(4F15.6))
  //C average |F|
  ftest = fem::max(1.0e-16, ftest / fem::ffloat(nf));
  //C LS mean
  frms = fem::sqrt(frms / fem::ffloat(nf) + 1.0e-32);
  if (iter == 0) {
    aiprin(cmn, x, vx, 3, iret);
  }
  //C      IF(NFIT.GT.1) WRITE(*,*) 'ISTATU=',ISTATU,FTEST
  if (istatu == 1) {
    goto statement_60;
  }
  //C     __________________________________________________________________
  //C     start numerical derivatives
  //C!!
  statement_30:
  istatu = -1;
  //C      IF(NFIT.GT.1) WRITE(*,*) 'start num'
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
  //C      IF(NFIT.GT.1) WRITE(*,*) 'do num',JRET
  iret = -1;
  //C...for constraint calculation
  if (jret < 0) {
    return;
  }
  //C     __________________________________________________________________
  //C     add elements of Jacobian for profile analysis
  //C 1-parameter derivative matrix
  statement_50:
  if (nadfs == 1) {
    FEM_DO_SAFE(j, 1, nx) {
      aux(j + nx * (nf - 1)) = 0.0e0;
    }
    aux(ipf + nx * (nf - 1)) = -1.0e0;
    //C 2-parameter derivative matrix
  }
  else if (nadfs == 2) {
    FEM_DO_SAFE(j, 1, nx) {
      aux(j + nx * (nf - 2)) = 0.0e0;
      aux(j + nx * (nf - 1)) = 0.0e0;
    }
    aux(ipfx + nx * (nf - 2)) = -1.0e0;
    aux(ipfy + nx * (nf - 1)) = -1.0e0;
  }
  //C     __________________________________________________________________
  //C     next iteration
  aniter(cmn, x, vx, fcopy, aux, xp, rh, aux(indwm + 1), dx);
  //C      WRITE(*,*) 'ANITER ',(AUX(INDST+J),J=1,NX)
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
  //C      WRITE(*,*) 'statement 80 in IPLOOP, NFIT=',NFIT
  if (nfit != 1) {
    goto statement_82;
  }
  //C      WRITE(*,*) 'INDDX,INDAS,INDWM,INDPU',INDDX,INDAS,INDWM,INDPU
  acopxv(cmn, x, vx, aux(cmn.inddx + 1), aux(indas + 1), aux(indwm + 1),
    aux(indpu + 1));
  //C      WRITE(*,*) 'vor AIPRIN'
  //C final print
  aiprin(cmn, x, vx, 4, iret);
  //C      WRITE(*,*) '>>>NSECA = ',NSECA
  if (nseca == 0) {
    return;
  }
  //C     __________________________________________________________________
  //C     start of profile analysis
  //C index of save Vx result
  indvs = indpu + nx;
  //C total number of words (so far)
  ndtot = indvs + (nx * nx + nx) / 2;
  ndende = ndtot;
  cmn.ndenda = ndende;
  FEM_DO_SAFE(j, 1, (nx * nx + nx) / 2) {
    //C saved primary cov. matrix
    aux(indvs + j) = vx(j);
    //C store back original VX
    vx(j) = aux(indwm + j);
  }
  chiref = chisq;
  iprsav = ipr;
  //C      WRITE(*,*) 'Original VX '
  //C      WRITE(*,*) (VX(J),J=1,6)
  //C      WRITE(*,*) 'Copy INDVS ',INDVS,' fitted matrix'
  //C      WRITE(*,*) (AUX(INDVS+J),J=1,6)
  //C     prepare storage for secondary fits
  //C     __________________________________________________________________
  //C     start next profile analysis
  statement_81:
  if (iseca >= nseca) {
    goto statement_99;
  }
  //C next profile analysis
  iseca++;
  if (npsec(2, iseca) == 0) {
    //C 1-parameter profile analysis
    nadfs = 1;
    nf = nfprim + nadfs;
    //C         WRITE(*,*) 'Calling A1PROF'
    //C         CALL A1PROF(X,VX,KRET)
    a1prof(kret);
  }
  else {
    //C 2-parameter profile analysis
    nadfs = 2;
    nf = nfprim + nadfs;
    //C        CALL A2PROF(X,VX,KRET)
    a2prof(kret);
  }
  istatu = 0;
  //C restart fit
  goto statement_10;
  //C     __________________________________________________________________
  //C     secondary fit for profile analysis ended, store result
  statement_82:
  if (nadfs == 1) {
    //C      WRITE(*,321) 'ST B1PROF V',(AUX(INDST+J),J=1,3)
    //C        CALL B1PROF(X,VX,KRET)          ! 1-parameter profile analysis
    //C 1-parameter profile analysis
    b1prof(kret);
    //C          WRITE(*,321) 'ST B1PROF N',(AUX(INDST+J),J=1,3)
  }
  else {
    //C        CALL B2PROF(X,VX,KRET)         ! 2-parameter profile analysis
    //C 2-parameter profile analysis
    b2prof(kret);
  }
  istatu = 0;
  //C next fit
  if (kret < 0) {
    goto statement_10;
  }
  goto statement_81;
  //C store primary fit results
  statement_99:
  chisq = chiref;
  //C printflag
  ipr = iprsav;
  FEM_DO_SAFE(j, 1, nx) {
    x(j) = aux(indas + j);
  }
  FEM_DO_SAFE(j, 1, (nx * nx + nx) / 2) {
    vx(j) = aux(indvs + j);
  }
  //C final print
  aiprin(cmn, x, vx, 4, iret);
  {
    write_loop wloop(cmn, 6, star);
    FEM_DO_SAFE(j, 1, 71) {
      wloop, "_";
    }
  }
  {
    write_loop wloop(cmn, 6, star);
    FEM_DO_SAFE(j, 1, 71) {
      wloop, "_";
    }
  }
  iret = 0;
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
  int& lunsim = cmn.lunsim;
  int& ipr = cmn.ipr;
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
  //C 380 "aploop.F" 2
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 381 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 382 "aploop.F" 2
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
    //C 388 "aploop.F" 2
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
        //C error condition
        if (ipr >= 2) {
          write(lunsim, star), "Variable", i, " is", x(i), " reset to normal";
        }
        //C reset: X(i) has to be positive
        ntvar = 0;
      }
      //C sqrt transformation - check
    }
    else if (ntvar == 5) {
      if (x(i) <= 0.0f) {
        //C error condition
        if (ipr >= 2) {
          write(lunsim, star), "Variable", i, " is", x(i), " reset to normal";
        }
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
    //C       CALL ATETOI  ! transform external to internal variables
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
    //C 463 "aploop.F" 2
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
  common_write write(cmn);
  const int naux = 125000;
  arr_ref<double> aux(cmn.aux, dimension(naux));
  int& nx = cmn.nx;
  int& lunsim = cmn.lunsim;
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
  //C 85 "aploop.F" 2
  //C
  //C 1 "comcfit.inc" 1
  //C     =========================== MACRO ================================
  //C     Parameter statement for basic dimensions
  //C     Definition of common for APLCON/ERRPRP/SIM... subroutines
  //C
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 86 "aploop.F" 2
  //C
  //C 1 "cprofil.inc" 1
  //C
  //C     Common for profile likelihood analysis with APLCON
  //C
  //C     NDEXT +1     is start of save area for primary X(.) and VX(.)
  //C     ILR    =     index of current profile point
  //C     NLR    =     number of profile points
  //C     ISECA  =     index of current profile analysis
  //C     NSECA  =     number of profile analyses
  //C     NFADD  =     max number of additional constraints
  //C     IPF    =     index of profile variable
  //C     NPSEC(2,.)   indices for profile analyses
  //C     XL,XR  =     limits for profile analysis
  //C     XLR(.) =     X values for profile
  //C     YLR(.) =     Y values for profile
  //C     FLR(.) =     chi^2 - reference value for profile
  //C     CONSTR =     fixed value
  //C     CONSTX,Y     fixed values
  //C     CENTER =     fitted value (minimum)
  //C     SIGMAX =     parabolic X error
  //C     SIGMAY =     parabolic Y error
  //C     CENTEX =     X center
  //C     CENTEY =     Y center
  //C     RED    =     reduction factor, e.g. 1/2
  //C     REDL   =     reduction factor left
  //C     REDR   =     reduction factor right
  //C     CPHI   =     cos phi of eigenvector
  //C     SPHI   =     sin phi
  //C     NSTAR  =     number of directions
  //C     ISTAR  =     direction index
  //C     CRT    =     cos of rotation transformation
  //C     SRT    =     sin of rotation transformation
  //C     CTH    =     cos rotation
  //C     STH    =     sin rotation
  //C     CHT    =     copy of cos rotation
  //C     XCONT(24,3)  x contur data, max 12 points, 3 contours
  //C     YCONT(24,3)  y contur data, max 12 points, 3 contours
  //C
  //C 87 "aploop.F" 2
  //C
  if (cmn.ncalls != 0) {
    goto statement_10;
  }
  //C     __________________________________________________________________
  //C     indices/pointer etc at first APLOOP entry
  //C max. number of constraints
  nff = cmn.nf + cmn.nfadd;
  //C total number of fit equations
  nxf = nx + nff;
  //C      WRITE(*,*) 'NFF,NXF=',NFF,NXF
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
  //C      WRITE(*,*) 'NDTOTL,NAUX=',NDTOTL,NAUX
  if (ndtotl > naux) {
    //C         WRITE(*,*) 'NDTOTL,NAUX=',NDTOTL,NAUX
    //C error printout
    aiprin(cmn, x, vx, 1, iret);
    FEM_STOP(0);
  }
  FEM_DO_SAFE(i, indlm + 2 * nx + 1, ndtot) {
    //C reset part of aux, unused so far
    aux(i) = 0.0e0;
  }
  //C initial steps ST(.)
  asteps(cmn, x, vx, aux(1 + cmn.indst));
  if (cmn.ipr != 0) {
    write(6, star), "Printout of initial values and correlation coeffs:";
    ciprv(lunsim, x, vx, nx);
    cfcorr(lunsim, vx, nx);
  }
  //C     __________________________________________________________________
  //C     internal APLOOP
  //C default status is -1 = continue
  statement_10:
  iret = -1;
  iploop(cmn, x, vx, f, aux(indxs + 1), aux(inddx + 1), aux(indcf + 1),
    aux(indxp + 1), aux(indrh + 1), iret);
}

//C
//C transform external to internal variable
void
atetoi(
  common& cmn,
  arr_cref<double> x,
  arr_ref<double> vx)
{
  x(dimension(star));
  vx(dimension(star));
  common_write write(cmn);
  // COMMON simcom
  int& nx = cmn.nx;
  int& ipak = cmn.ipak;
  int& indtr = cmn.indtr;
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 861 "aploop.F" 2
  //C
  //C 1 "nauxfit.inc" 1
  //C
  //C     Auxiliary array for all arrays used in APLCON
  //C 1 Mega byte
  //C
  //C 862 "aploop.F" 2
  //C
  //C 1 "declarefl.inc" 1
  //C
  //C     flags used in packfl.inc, unpackfl.inc
  //C
  //C 863 "aploop.F" 2
  //C     ...
  //C transformation back
  int i = fem::int0;
  int ntrfl = fem::int0;
  int ntvar = fem::int0;
  int ntmes = fem::int0;
  int ntder = fem::int0;
  int ntine = fem::int0;
  int ntlim = fem::int0;
  int ntprf = fem::int0;
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
    //C 867 "aploop.F" 2
    //C      transform covariance matrix for transformed variables
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
  }
  return;
  //C
  //C transform internal to external variable
  // UNHANDLED: ENTRY atitoe(x,vx)
  //C transformation back
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
    //C 888 "aploop.F" 2
    if (ntvar == 4) {
      write(6, star), "cov matrix back with ", x(i);
      //C log-normal
      FEM_DO_SAFE(j, 1, nx) {
        vx(ijsym(i, j)) = vx(ijsym(i, j)) * x(i);
        if (i == j) {
          vx(ijsym(i, j)) = vx(ijsym(i, j)) * x(i);
        }
      }
    }
    else if (ntvar == 5) {
      //C sqrt
      FEM_DO_SAFE(j, 1, nx) {
        vx(ijsym(i, j)) = vx(ijsym(i, j)) * 2.0e0 * fem::sqrt(x(i));
        if (i == j) {
          vx(ijsym(i, j)) = vx(ijsym(i, j)) * 2.0e0 * fem::sqrt(x(i));
        }
      }
    }
  }
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
  //C     NADFS  = 0   initial (primary) fit (set by APLCON)
  //C            = 1   1-parameter analysis  (set by BPLCON)
  //C            = 2   2-parameter analysis      "..."
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
  //C     LUNSIM   = printout unit (see parameter statement)
  //C     IDEBUG   = debug flag value
  //C                = 0   no printout
  //C                = 1   only warnings are printed (default)
  //C                = 2,3 more and more printout
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
  //C 1072 "aploop.F" 2
  //C     ...
  //C chi^square
  chi2 = chisq;
  //C number of degrees of freedom
  nd = ndf;
  //C p-value
  pval = chprob(chisq, ndf);
}

} // namespace placeholder_please_replace
