// firfb_process.c - FIR-filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

// complex multiply: z = x * y
static __inline void
cmul(float *z, float *x, float *y, int n)
{
    int      i, ir, ii;

    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = i * 2 + 1;
        z[ir] = x[ir] * y[ir] - x[ii] * y[ii];
        z[ii] = x[ir] * y[ii] + x[ii] * y[ir];
    }
}

// FIR-filterbank analysis for short chunk (cs < nw)
static __inline void
firfb_analyze_sc(float *x, float *y, int cs,
    float *hh, float *xx, float *yy, float *zz, int nc, int nw)
{
    float   *hk, *yk, *zk;
    int      i, j, k, nf, ns, nt, nk;

    nk = nw / cs;
    nt = cs * 2;
    nf = cs + 1;
    ns = nf * 2;
    // loop over channels
    for (k = 0; k < nc; k++) {
        fzero(xx, nt);
        fcopy(xx, x, cs);
        cha_fft_rc(xx, nt);
        // loop over sub-window segments
        yk = y + k * cs;
        zk = zz + k * (nw + cs);
        for (j = 0; j < nk; j++) {
            hk = hh + (k * nk + j) * ns;
            cmul(yy, xx, hk, nf);
            cha_fft_cr(yy, nt);
            for (i = 0; i < nt; i++) {
                zk[i + j * cs] += yy[i];
            }
        }
        fcopy(yk, zk, cs);
        fmove(zk, zk + cs, nw);
        fzero(zk + nw, cs);
    }
}

// FIR-filterbank analysis for long chunk (cs >= nw)
static __inline void
firfb_analyze_lc(float *x, float *y, int cs,
    float *hh, float *xx, float *yy, float *zz, int nc, int nw)
{
    float   *hk, *yk, *zk;
    int      i, j, k, nf, nt, ni;

    nt = nw * 2;
    nf = nw + 1;
    // loop over sub-chunk segments
    for (j = 0; j < cs; j += nw) {
        ni = ((cs - j) < nw) ? (cs - j) : nw;
        fzero(xx, nt);
        fcopy(xx, x + j, ni);
        cha_fft_rc(xx, nt);
       // loop over channels
        for (k = 0; k < nc; k++) {
            hk = hh + k * nf * 2;
            cmul(yy, xx, hk, nf);
            cha_fft_cr(yy, nt);
            yk = y + k * cs;
            zk = zz + k * nw;
            for (i = 0; i < ni; i++) {
                yk[i + j] = yy[i] + zk[i];
            }
            fcopy(zk, yy + ni, nw);
        }
    }
}

// FIR-filterbank analysis
FUNC(void)
cha_firfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
{
    float   *hh, *xx, *yy, *zz;
    int      nc, nw;

    nc = CHA_IVAR[_nc];
    nw = CHA_IVAR[_nw];
    hh = (float *) cp[_ffhh];
    xx = (float *) cp[_ffxx];
    yy = (float *) cp[_ffyy];
    zz = (float *) cp[_ffzz];
    if (cs < nw) {
        firfb_analyze_sc(x, y, cs, hh, xx, yy, zz, nc, nw);
    } else {
        firfb_analyze_lc(x, y, cs, hh, xx, yy, zz, nc, nw);
    }
}

// FIR-filterbank synthesis
FUNC(void)
cha_firfb_synthesize(CHA_PTR cp, float *x, float *y, int cs)
{
    float    xsum;
    int      i, k, nc;

    nc = CHA_IVAR[_nc];
    for (i = 0; i < cs; i++) {
        xsum = 0;
        for (k = 0; k < nc; k++) {
            xsum += x[i + k * cs];
        }
        y[i] = xsum;
    }
}
