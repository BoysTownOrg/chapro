// iirfb_process.c - IIR-filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"
#include "cha_if.h"

#define SOS

/***********************************************************/

#ifdef SOS
static __inline void
filter_sos(float *x, float *y, int cs, float *coef, float *hist, int nsos)
{
#ifdef __arm__
    struct {int nsos; float *hist; float *coef;} S = {nsos, hist, coef};
    arm_biquad_cascade_df1_f32(&S, x, y, cs);
#else
    float b0, b1, b2, a1, a2, x0, x1, x2, y0, y1, y2;
    int i, j;

    for (j = 0; j < nsos; j++) {
        b0 = *coef++;
        b1 = *coef++;
        b2 = *coef++;
        a1 = *coef++;
        a2 = *coef++;
        x1 = hist[0];
        x2 = hist[1];
        y1 = hist[2];
        y2 = hist[3];
        for (i = 0; i < cs; i++) {
            x0 = j ? y[i] : x[i];
            y0 = b0 * x0 + b1 * x1 + b2 * x2 + a1 * y1 + a2 * y2;
            x2 = x1;
            x1 = x0;
            y2 = y1;
            y1 = y0;
            y[i] = y0;
        }
        *hist++ = x1;
        *hist++ = x2;
        *hist++ = y1;
        *hist++ = y2;
    }
#endif // __arm__
}
#else // SOS
static __inline void
filter_tf(float *x, float *y, int cs, float *coef, float *hist, int nz)
{
    float  *b, *a, *xx, *yy, y0;
    int i, j;

    b = coef;
    a = b + nz + 1;
    xx = hist;
    yy = xx + nz;
    /* loop over time */
    for (i = 0; i < cs; i++) {
        y0 = b[0] * x[i];
        for (j = 1; j <= nz; j++) {
            y0 += b[j] * xx[j - 1] - a[j] * yy[j - 1];
        }
        fmove(xx + 1, xx, nz - 1);
        fmove(yy + 1, yy, nz - 1);
        xx[0] = x[i];
        yy[0] = y[i] = y0;
    }
}
#endif // SOS

/***********************************************************/

// IIR-filterbank analysis
FUNC(void)
cha_iirfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
{
    float   *bb, *zz, *yk, *yd, *bk, *zk;
    int     dk, i, k, m, nc, ns, nz, op, ncoef, nhist, *dd;

    bb = (float *) cp[_bb];
    dd = (int *) cp[_dd];
    zz = (float *) cp[_zz];
    yd = (float *) cp[_yd];
    nc = CHA_IVAR[_nc];
    op = CHA_IVAR[_op];
    ns = CHA_IVAR[_ns];
    nz = op - 1;
    nhist = 2 * nz;
#ifdef SOS
    ncoef = 5 * (nz / 2);
#else
    ncoef = 2 * op;
#endif
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++) {
        yk = y + k * cs;
        bk = bb + k * ncoef;
        zk = zz + k * nhist;
#ifdef SOS
        filter_sos(x, yk, cs, bk, zk, nz / 2);
#else
        filter_tf(x, yk, cs, bk, zk, nz);
#endif
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
