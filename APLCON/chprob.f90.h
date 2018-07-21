#include <cmath>

using namespace std;

inline double dgamml_(double x) {
    /* Initialized data */

    constexpr double cof[6] = {76.18009172947146,   -86.50532032941677,
                                24.01409824083091,   -1.231739572450155,
                                .001208650973866179, -5.395239384953e-7};
    constexpr double stp = 2.5066282746310005;

    /* ln[Gamma(x)] */
    /*     ... */
    double yy = x;
    double tmp = (x + .5) * log(x + 5.5) - x - 5.5;
    double ser = 1.000000000190015;
    for (int j = 1; j <= 6; ++j) {
        yy += 1.;
        ser += cof[j - 1] / yy;
    }
    return tmp + log(stp * ser / x);
} /* dgamml_ */

inline double dgamin_(double a, double x) {
    /*     incomplete gamma function P(a,x) */
    /*     returns -1.0 for x < 0 or a =< 0 */
    /*     ... */
    if (x < 0. || a <= 0.) {
        return -1.;
        /* error */
    } else if (x == 0.f) {
        return 0.;
    } else if (x < a + 1.) {
        /* series representation */
        double gln = dgamml_(a);
        /* ln[Gamma(a)] */
        double ap = a;
        double sum = 1. / a;
        double del = sum;
        for (int n = 1; n <= 300; ++n) {
            ap += 1.f;
            del = del * x / ap;
            sum += del;
            if (abs(del) < abs(sum) * 3e-7) {
                goto L10;
            }
        }
L10:
        return sum * exp(-(x) + a * log(x) - gln);
    } else {
        /* continued fraction representation */
        double gln = dgamml_(a);
        /* ln[Gamma(a)] */
        double b = x + 1. - a;
        double c = 9.9999999999999988e29;
        double d = 1. / b;
        double h = d;
        for (int i = 1; i <= 300; ++i) {
            double an = -i * (i - a);
            b += 2.;
            d = an * d + b;
            if (abs(d) < 1e-30) {
                d = 1e-30;
            }
            c = b + an / c;
            if (abs(c) < 1e-30) {
                c = 1e-30;
            }
            d = 1. / d;
            double del = d * c;
            h *= del;
            if (abs(del - 1.) < 3e-7) {
                goto L20;
            }
        }
        return 0;
L20:
        return 1. - exp(-(x) + a * log(x) - gln) * h;
    }
} /* dgamin_ */

inline double chprob_(double chisq, int n) {
    /*     chi square probability for N degrees of freedom at CHISQ */
    /*     =integral from 0 ... CHISQ */
    /* prob from chisquare */
    /*     ... */
    if (chisq <= 0.f) {
        return 1.;
    } else {
        return 1. - dgamin_(n * .5, chisq * .5);
    }
} /* chprob_ */




