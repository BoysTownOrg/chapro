// firfb_prepare.c - FIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

// compute FIR-filterbank coefficients
static __inline void
fir_filterbank(float *bb, double *cf, int nc, int nw, int wt, double sr)
{
    double   p, w, a = 0.16, sm = 0;
    float   *ww, *bk, *xx, *yy;
    int      j, k, kk, nt, nf, ns, *be;

    nt = nw * 2;
    nf = nw + 1;
    ns = nf * 2;
    be = (int *) calloc(nc + 1, sizeof(int));
    ww = (float *) calloc(nw, sizeof(float));
    xx = (float *) calloc(ns, sizeof(float));
    yy = (float *) calloc(ns, sizeof(float));
    // window
    for (j = 0; j < nw; j++) {
        p = M_PI * (2.0 * j - nw) / nw;
        if (wt == 0) {
            w = 0.54 + 0.46 * cos(p);                   // Hamming
        } else {
            w = (1 - a + cos(p) + a * cos(2 * p)) / 2;  // Blackman
        }
        sm += w;
        ww[j] = (float) w;
    }
    // frequency bands
    be[0] = 0;
    for (k = 1; k < nc; k++) {
        kk = round(nf * cf[k - 1] * (2 / sr));
        be[k] = (kk > nf) ? nf : kk;
    }
    be[nc] = nf;
    // channel tranfer functions
    fzero(xx, ns);
    xx[nw / 2] = 1;
    cha_fft_rc(xx, nt);
    for (k = 0; k < nc; k++) {
        bk = bb + k * nw;
        fzero(yy, ns);
        fcopy(yy + be[k] * 2, xx + be[k] * 2, (be[k + 1] - be[k]) * 2);
        cha_fft_cr(yy, nt);
	// apply window to iFFT of bandpass
        for (j = 0; j < nw; j++) {
            yy[j] *= ww[j];
        }
	fcopy(bk, yy, nw);
    }
    free(be);
    free(ww);
    free(xx);
    free(yy);
}

// Fourier-transform FIR coefficients for short chunk (cs < nw)
static __inline void
fir_transform_sc(CHA_PTR cp, float *bb, int nc, int nw, int cs)
{
    float *bk, *hk, *hh;
    int j, k, nk, nf, ns, nt;

    nk = nw / cs;
    nt = cs * 2;
    nf = cs + 1;
    ns = nf * 2;
    cha_allocate(cp, nc * ns * nk, sizeof(float), _ffhh);
    hh = cp[_ffhh];
    for (k = 0; k < nc; k++) {
        for (j = 0; j < nk; j++) {
            bk = bb + k * nw + j * cs;
            hk = hh + (k * nk + j) * ns;
            fcopy(hk, bk, cs);
            cha_fft_rc(hk, nt);
        }
    }
}

// Fourier-transform FIR coefficients for long chunk (cs >= nw)
static __inline void
fir_transform_lc(CHA_PTR cp, float *bb, int nc, int nw, int cs)
{
    float *hk, *hh;
    int k, ns, nt;

    nt = nw * 2;
    ns = nt + 2;
    cha_allocate(cp, nc * ns, sizeof(float), _ffhh);
    hh = cp[_ffhh];
    for (k = 0; k < nc; k++) {
        hk = hh + k * ns;
        fcopy(hk, bb + k * nw, nw);
        cha_fft_rc(hk, nt);
    }
}

/***********************************************************/

FUNC(int)
cha_firfb_prepare(CHA_PTR cp, double *cf, int nc, double fs, 
    int nw, int wt, int cs)
{
    float   *bb;
    int      ns, nt;

    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = fs;
    // allocate window buffers
    CHA_IVAR[_nw] = nw;
    CHA_IVAR[_nc] = nc;
    nt = nw * 2;
    ns = nt + 2;
    cha_allocate(cp, ns, sizeof(float), _ffxx);
    cha_allocate(cp, ns, sizeof(float), _ffyy);
    cha_allocate(cp, nc * (nw + cs), sizeof(float), _ffzz);
    // compute FIR-filterbank coefficients
    bb = calloc(nc * nw, sizeof(float));
    fir_filterbank(bb, cf, nc, nw, wt, fs);
    // Fourier-transform FIR coefficients
    if (cs < nw) {  // short chunk
        fir_transform_sc(cp, bb, nc, nw, cs);
    } else {        // long chunk
        fir_transform_lc(cp, bb, nc, nw, cs);
    }
    free(bb);

    return (0);
}
