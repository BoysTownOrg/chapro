// iirfb_process.c - IIR-filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"
#include "cha_if.h"

/***********************************************************/

static __inline void
filter_tf(float *x, float *y, int cs, float *coef, float *hist, int op)
{
    float  *b, *a, *xx, *yy, y0;
    int i, j;

    b = coef;
    a = b + op;
    xx = hist;
    yy = xx + op;
    /* loop over time */
    for (i = 0; i < cs; i++) {
        y0 = b[0] * x[i];
        for (j = 1; j < op; j++) {
            y0 += b[j] * xx[j - 1] - a[j] * yy[j - 1];
        }
        fmove(xx + 1, xx, op - 2);
        fmove(yy + 1, yy, op - 2);
        xx[0] = x[i];
        yy[0] = y[i] = y0;
    }
}

/***********************************************************/

// IIR-filterbank analysis
FUNC(void)
cha_iirfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
{
    int     dk, i, k, kk, m, nc, ns, op, *dd;
    float   *bb, *zz, *yk, *yd;

    bb = (float *) cp[_bb];
    dd = (int *) cp[_dd];
    zz = (float *) cp[_zz];
    yd = (float *) cp[_yd];
    nc = CHA_IVAR[_nc];
    op = CHA_IVAR[_op];
    ns = CHA_IVAR[_ns];
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++) {
        kk = k * op * 2;
        yk = y + k * cs;
        filter_tf(x, yk, cs, bb + kk, zz + kk, op);
       /* delay */
        m = k * ns;
        dk = dd[k];
        for (i = 0; i < cs; i++) {
            fmove(yd + m + 1, yd + m, dk);
            yd[m] = yk[i];
            yk[i] = yd[m + dk];
        }
    }
}

// IIR-filterbank synthesis
FUNC(void)
cha_iirfb_synthesize(CHA_PTR cp, float *x, float *y, int cs)
{
    float   xsum;
    int     i, k, nc;

    nc = CHA_IVAR[_nc];
    for (i = 0; i < cs; i++) {
        xsum = 0;
        // loop over filterbank channel
        for (k = 0; k < nc; k++) {
            xsum += x[i + k * cs];
        }
        y[i] = xsum;
    }
}
