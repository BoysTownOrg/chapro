// cfirfb_prepare.c - complex FIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

/***********************************************************/

// expand real to complex
static void
rc(float *x, int n)
{
    int i, j, jr, ji;

    for (i = 0; i < n; i++) {
        j = (n - 1) - i;
        jr = 2 * j;
        ji = 2 * j + 1;
        x[jr] = x[j];
        x[ji] = 0;
    }
}

// generate analytic signal
static void
anasig(float *x, int n)
{
    int i;

    rc(x, n);
    cha_fft(x, n);
    for (i = 2; i < n; i++) {
        x[i] *= 2;
        x[i + n] = 0;
    }
    cha_ifft(x, n);
}

/***********************************************************/

// compute FIR-filterbank coefficients
static __inline void
fir_filterbank(float *bb, double *cf, int nc, int nw, int wt, double sr)
{
    double   p, w, a = 0.16;
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
fir_transform_sc(float *bb, float *hh, int nc, int nw, int cs)
{
    float *bk, *hk;
    int j, k, nk, ns, nt;

    nk = nw / cs;
    nt = cs * 2;
    ns = nt * 2;
    bk = (float *) calloc(nw * 2, sizeof(float));
    for (k = 0; k < nc; k++) {
        fcopy(bk, bb + k * nw, nw);
        anasig(bk, nw);
        for (j = 0; j < nk; j++) {
            hk = hh + (k * nk + j) * ns;
            fcopy(hk, bk + j * nt, nt);
            cha_fft(hk, nt);
        }
    }
    free(bk);
}

// Fourier-transform FIR coefficients for long chunk (cs >= nw)
static __inline void
fir_transform_lc(float *bb, float *hh, int nc, int nw, int cs)
{
    float *hk;
    int k, ns, nt;

    nt = nw * 2;
    ns = nt * 2;
    for (k = 0; k < nc; k++) {
        hk = hh + k * ns;
        fcopy(hk, bb + k * nw, nw);
        anasig(hk, nt);
        cha_fft(hk, nt);
    }
}

/***********************************************************/

// save bandwidth
static __inline void
save_bw(CHA_PTR cp, double *cf, int nc)
{
    float *fck, *bwk, pcf;
    int i;

    fck = (float *) cha_allocate(cp, nc, sizeof(float), _fc);
    bwk = (float *) cha_allocate(cp, nc, sizeof(float), _bw);
    pcf = 0;
    for (i = 0; i < nc; i++) {
        fck[i] = (float) (cf[i] + pcf) / 2000;
        bwk[i] = (float) (cf[i] - pcf) / 1000;
        pcf = (float) cf[i];
    }
}

/***********************************************************/

FUNC(int)
cha_cfirfb_prepare(CHA_PTR cp, double *cf, int nc, double sr, 
    int nw, int wt, int cs)
{
    float   *bb, *hh;
    int      ns, nt, nh;

    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = sr / 1000;
    // allocate window buffers
    CHA_IVAR[_nw] = nw;
    CHA_IVAR[_nc] = nc;
    nt = nw * 2;
    ns = nt + 2;
    nh = (cs < nw) ? (cs * 4) * (nw / cs) : (nw * 4);
    cha_allocate(cp, ns * 2, sizeof(float), _ffxx);
    cha_allocate(cp, ns * 2, sizeof(float), _ffyy);
    cha_allocate(cp, nc * (nw + cs) * 2, sizeof(float), _ffzz);
    cha_allocate(cp, nc * nh, sizeof(float), _ffhh);
    // save bandwidths
    save_bw(cp, cf, nc);
    // compute FIR-filterbank coefficients
    bb = calloc(nc * nw, sizeof(float));
    fir_filterbank(bb, cf, nc, nw, wt, sr);
    // Fourier-transform FIR coefficients
    hh = cp[_ffhh];
    if (cs < nw) {  // short chunk
        fir_transform_sc(bb, hh, nc, nw, cs);
    } else {        // long chunk
        fir_transform_lc(bb, hh, nc, nw, cs);
    }
    free(bb);
    // allocate chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);

    return (0);
}
