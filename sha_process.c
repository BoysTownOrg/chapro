// sha_process.c - nonlinear-frequency-compression processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"

static double g0, a1, a2, a3;
static float *AA, *II, *SS;
static int exr, hbw;

/***********************************************************/
#ifdef ARM_DSP

static int fft_initialized = 0;
static arm_rfft_fast_instance_f32 fft_instance;

static __inline void
rfft(float *x, int n)
{
    if (!fft_initialized) {
        arm_rfft_fast_init_f32(&fft_instance, n);
	fft_initialized++;
    }
    arm_rfft_fast_f32(&fft_instance, x, x, 0);
}

static __inline void
rifft(float *x, int n)
{
    if (!fft_initialized) {
        arm_rfft_fast_init_f32(&fft_instance, n);
	fft_initialized++;
    }
    arm_rfft_fast_f32(&fft_instance, x, x, 1);
}

#else // ARM_DSP

#define rfft(x,n)    cha_fft_rc(x,n)
#define rifft(x,n)   cha_fft_cr(x,n)

#endif // ARM_DSP

// apply compression
static __inline void
compress(float *y, float *x, int nf, float *g1)
{
    float g, gg, xr, xi, xx, JJ;
    int e, k, kr, ki, k1, k2, k2mn, k2mx, kk;
    static float spl_ref = 1.1219e-6;
    static float eps = 1e-12;

    gg = 1 / (spl_ref * spl_ref) / 4; // reference intensity [11-Jul-2021]
    // compute amplitude and intensity
    for (k = 0; k < nf; k++) {
        kr = 2 * k;
        ki = kr + 1;
        xr = x[kr];
        xi = x[ki];
        if (g1) {
            xr *= g1[k];
            xi *= g1[k];
        }
        xx = (xr * xr + xi * xi) * gg;
        II[k] = xx;
        AA[k] = (float)sqrt(xx);
    }
    // apply suppression
    if (hbw) {
        for (k1 = 0; k1 < nf; k1++) {
            y[k1] = 0;
            k2mn = k1 - hbw;
            k2mx = k1 + hbw;
            k2mn = (k2mn <  0) ?  0 : k2mn;
            k2mx = (k2mx > nf) ? nf : k2mx;
            for (k2 = k2mn; k2 < k2mx; k2++) {
                kk = k1 + k2 * nf;
                y[k1] += (k2 == k1) ? SS[kk] : SS[kk] * II[k2];
            }
        }
        for (k = 0; k < nf; k++) {
            xx = y[k];
            II[k1] = xx;
            AA[k1] = (float)sqrt(xx);
        }
    }
    // calculate compression and apply to input
    for (k = 0; k < nf; k++) {
        if (exr > 1) {                  // include expansion ??
            JJ = 1 / (II[k] + eps); 
            e = exr;
            while (--e) JJ *= JJ;      // repeat xr times 
        } else {
            JJ = 0;
        }
        g = g0 / (float)sqrt(1 + a1 * AA[k] + a2 * II[k] + a3 * JJ);
        kr = 2 * k;
        ki = kr + 1;
        y[kr] = x[kr] * g; // real
        y[ki] = x[ki] * g; // imag
    }
}

// short-term FFT analyze
static __inline void
short_term_analyze(float *xx, float *XX, float *ww, int nw, int ns, int nf)
{
    int i;

    for (i = 0; i < nw; i++) {
        XX[i] = xx[i] * ww[i];   // apply window to input
    }
    fzero(XX + nw, nw);
    rfft(XX, nf);                // FFT
    fcopy(xx, xx + ns, ns);      // save last half of input window
}

// short-term FFT synthesize
static __inline void
short_term_synthesize(float *yy, float *YY, int ns, int nf)
{
    int i, nn;

    rifft(YY, nf);               // IFFT
    nn = nf - ns;
    fmove(yy, yy + ns, nn);      // shift previous output
    for (i = 0; i < nn; i++) {
        yy[i] += YY[i];          // overlap-add output
    }
    fcopy(yy + nn, YY + nn, ns); // save response tail
}

// nonlinear-frequency-compression short chunk
static __inline void
sha_sc(CHA_PTR cp, float *x, float *y, int cs,
    float *xx, float *yy, float *XX, float *YY, float *ww,
    float *g1, int nw)
{
    int icp, ics, nn, nf, ns, ncs;

    // process chunk
    ncs = CHA_IVAR[_sha_ncs];
    ics = CHA_IVAR[_sha_ics];
    ns = nw / 2;
    nf = nw * 2;
    nn = ics * cs;
    fcopy(xx + nn + ns, x, cs);
    fcopy(y, yy + nn, cs);
    icp = (ics + 1) % ncs;
    if (icp == 0) { // perform frequency-map after every shift
        short_term_analyze(xx, XX, ww, nw, ns, nf);
        compress(YY, XX, nf, g1);
        short_term_synthesize(yy, YY, ns, nf);
    }
    // update chunk count
    CHA_IVAR[_sha_ics] = icp;
}

// nonlinear-frequency-compression long chunk
static __inline void
sha_lc(float *x, float *y, int cs,
    float *xx, float *yy, float *XX, float *YY, float *ww, 
    float *g1, int nw)
{
    int k, nn, nf, ns;

    // process chunk
    ns = nw / 2;
    nf = nw * 2;
    nn = cs / ns;
    for (k = 0; k < nn; k++) {
        fcopy(xx + nn + ns, x + k * ns, ns);
        fcopy(y + k * ns, yy + nn, ns);
        // perform frequency-map after every shift
        short_term_analyze(xx, XX, ww, nw, ns, nf);
        compress(YY, XX, nf, g1);
        short_term_synthesize(yy, YY, ns, nf);
    }
}

/***********************************************************/

// nonlinear-frequency-compression analysis
FUNC(void)
cha_sha_process(CHA_PTR cp, float *x, float *y, int cs)
{
    float *ww, *xx, *yy, *XX, *YY, *g1;
    int nw;

    // copy parameters and pointers from cha_data
    nw = CHA_IVAR[_sha_nw];
    ww = (float *) cp[_sha_ww];
    xx = (float *) cp[_sha_xx];
    yy = (float *) cp[_sha_yy];
    XX = (float *) cp[_sha_XX];
    YY = (float *) cp[_sha_YY];
    g1 = (float *) cp[_sha_g1];
    SS = (float *) cp[_sha_SS];
    AA = (float *) cp[_sha_AA];
    II = (float *) cp[_sha_II];
    g0 = CHA_DVAR[_sha_g0];
    a1 = CHA_DVAR[_sha_a1];
    a2 = CHA_DVAR[_sha_a2];
    a3 = CHA_DVAR[_sha_a3];
    exr = CHA_IVAR[_sha_xr];
    hbw = CHA_IVAR[_sha_hbw];
    if (cs <= (nw / 2)) {       // short chunk ??
        sha_sc(cp, x, y, cs, xx, yy, XX, YY, ww, g1, nw);
    } else {                   // long chunk (not yet implemented)
        sha_lc(x, y, cs, xx, yy, XX, YY, ww, g1, nw);
    }
}
