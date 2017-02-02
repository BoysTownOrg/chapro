// db.c - convert magnitudes to and from decibels
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"
#include "cha_ff.h"

#undef EXACT // define to check exact values

/***********************************************************/

#ifndef EXACT

static __inline float
pow_ap(float x) // approximate 2^x
{
    float p, q;
    int   n;
    static float a0 = -206059.514;
    static float a1 = -72102.2578;
    static float a2 = -11240.0288;
    static float a3 = -989.027847;
    static float b0 = -206059.514;
    static float b1 = 70727.3131;
    static float b2 = -10763.5093;
    static float b3 = -14.5169712;

    // assume: x > -0.5 and x < 0.5
    n = (x < 0);
    if (n) x = -x;
    p = (((a3 * x + a2) * x + a1) * x + a0);
    q = (((b3 * x + b2) * x + b1) * x + b0);
    if (n) return (q / p);
    return (p / q);
}

static __inline float
log_ap(float x) // approximate ln(x)
{
    float p, q, z, z2;
    static float a0 = 75.1518561;
    static float a1 = -134.730400;
    static float a2 = 74.2011014;
    static float b0 = 37.5759281;
    static float b1 = -79.8905092;
    static float b2 = 56.2155348;

    // assume: x > sqrt(1/2) and x < sqrt(2)
    z = (x - 1) / (x + 1);
    z2 = z * z;
    p = ((a2 * z2 + a1) * z2 + a0);
    q = ((b2 * z2 + b1) * z2 + b0);
    return (z * p / q);
}

#endif // EXACT

/***********************************************************/

FUNC(float)
cha_db1(float x) // 10 * log10(x)
{
    float m, ln;
    int e;
    static float c1 = 1e-38;
    static float c2 = -380;
    static float c3 = 1e38;
    static float c4 = 380;
    static float c5 = 4.34294462;  // 10 / log(10)
    static float c6 = 1.41421354;  // sqrt(2)
    static float c7 = 0.693147182; // log(2);

    if (x < c1) return (c2);
    if (x > c3) return (c4);
    m = frexpf(x, &e);
    if (m < c6) {
        m *= 2;
	e--;
    }
#ifdef EXACT
    ln = logf(m); // exact
#else
    ln = log_ap(m);
#endif
    return (c5 * ln + c7 * e);
}

FUNC(float)
cha_undb1(float x) // 10 ^ (x / 10)
{
    float m, p2;
    int   e, n;
    static float c1 = 1e-38;
    static float c2 = -380;
    static float c3 = 1e38;
    static float c4 = 380;
    static float c5 = 0.332192808;  // log(10) / (10 * log(2))

    if (x < c2) return (c1);
    if (x > c4) return (c3);
    n = (x < 0);
    if (n) x = -x;
    x *= c5;
    e = (int) x;
    m = x - e;
    if (n) {
	e = -e;
	m = -m;
    }
#ifdef EXACT
    p2 = powf(2, m); // exact
#else
    p2 = pow_ap(m); 
#endif
    return (ldexpf(p2, e));
}

/***********************************************************/

FUNC(float)
cha_db2(float x) // 20 * log10(x)
{
    float m, ln;
    int e;
    static float c1 = 1e-38;
    static float c2 = -760;
    static float c3 = 1e38;
    static float c4 = 760;
    static float c5 = 8.68588924;  // 20 / log(10)
    static float c6 = 1.41421354;  // sqrt(2)
    static float c7 = 0.693147182; // log(2);

    if (x < c1) return (c2);
    if (x > c3) return (c4);
    m = frexpf(x, &e);
    if (m < c6) {
        m *= 2;
	e--;
    }
#ifdef EXACT
    ln = logf(m); // exact
#else
    ln = log_ap(m);
#endif
    return (c5 * ln + c7 * e);
}

FUNC(float)
cha_undb2(float x) // 10 ^ (x / 20)
{
    float m, p2;
    int   e, n;
    static float c1 = 1e-38;
    static float c2 = -760;
    static float c3 = 1e38;
    static float c4 = 760;
    static float c5 = 0.166096404;  // log(10) / (20 * log(2))

    if (x < c2) return (c1);
    if (x > c4) return (c3);
    n = (x < 0);
    if (n) x = -x;
    x *= c5;
    e = (int) x;
    m = x - e;
    if (n) {
	e = -e;
	m = -m;
    }
#ifdef EXACT
    p2 = powf(2, m); // exact
#else
    p2 = pow_ap(m); 
#endif
    return (ldexpf(p2, e));
}
