#ifndef APLCON_wrapper_h
#define APLCON_wrapper_h

#include <vector>

// main routines
/**
 * @brief Initialize APLCON
 * @param NVAR number of variables
 * @param MCST number of constraints
 */
void c_aplcon_aplcon(const int NVAR, const int MCST);
/**
 * @brief Loop function for fit
 * @param X current variable values
 * @param VX current covariance matrix (symmetrized)
 * @param F current constraint values
 * @param IRET status of fit iteration
 */
void c_aplcon_aploop(std::vector<double>& X, std::vector<double>& VX, double F[], int* IRET);

// routines to obtain results
/**
 * @brief Obtain ChiSquare, NDoF, Probability
 * @see c_aplcon_apstat
 * @param CHI2 ChiSquare as float
 * @param ND Number of degrees of freedom
 * @param PVAL Probability
 */
void c_aplcon_chndpv(double* CHI2, int* ND, double* PVAL);
/**
 * @brief Obtain information after fit
 * @param FOPT ChiSqare as double
 * @param NFUN number of function calls
 * @param NITER number of iterations
 */
void c_aplcon_apstat(double* FOPT, int* NFUN, int* NITER);
/**
 * @brief Obtain pulls
 * @param PULLS Array of pulls for each variable in X
 */
void c_aplcon_appull(double* PULLS);

// setup/config routines
/**
 * @brief Setup constraint accuracy
 * @param EPSF Constraint accuracy
 */
void c_aplcon_apdeps(const double EPSF);
/**
 * @brief Setup stepsize for variable with index I
 * @param I index of variable (starting from 1)
 * @param STEP stepsize of variable, 0 for fixed
 */
void c_aplcon_apstep(const int I, const double STEP);
/**
 * @brief Setup variable I to be poissionian distributed
 * @param I index of variable
 */
void c_aplcon_apoiss(const int I);

#endif
