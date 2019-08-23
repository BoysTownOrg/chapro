// cfirfb_process.c - complex FIR-filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"

/***********************************************************/

#ifdef ARM_DSP
#define cmul(z,x,y,n)   arm_cmplx_mult_cmplx_f32(x,y,z,n)

static int fft_initialized = 0;
static arm_cfft_instance_f32 fft_instance;

static __inline void
cfft_init(arm_cfft_instance_f32 *S, int length)
  switch (length) {
    case 32:
      *S = arm_cfft_sR_f32_len32;
      break;
    case 64:
      *S = arm_cfft_sR_f32_len64;
      break;
    case 128:
      *S = arm_cfft_sR_f32_len128;
      break;
    case 256:
      *S = arm_cfft_sR_f32_len256;
      break;
  }
}

static __inline void
fft(float *x, int n)
{
    if (!fft_initialized) {
        cfft__init(&fft_instance, n);
	fft_initialized++;
    }
    arm_cfft_f32(&fft_instance, x, x, 0, 1);
}

static __inline void
ifft(float *x, int n)
{
    if (!fft_initialized) {
        cfft__init(&fft_instance, n);
	fft_initialized++;
    }
    arm_cfft_f32(&fft_instance, x, x, 1, 1);
}
#else
#define fft(x,n)    cha_fft(x,n)
#define ifft(x,n)   cha_ifft(x,n)

static __inline void
cmul(float *z, float *x, float *y, int n)
{
    register float zr, zi;
    int      i, ir, ii;

// complex multiply: z = x * y
    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = i * 2 + 1;
        zr = x[ir] * y[ir] - x[ii] * y[ii];
        zi = x[ir] * y[ii] + x[ii] * y[ir];
        z[ir] = zr;
        z[ii] = zi;
    }
}
#endif

/***********************************************************/

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

/***********************************************************/

// FIR-filterbank analysis for short chunk (cs < nw)
static __inline void
firfb_analyze_sc(float *x, float *y, int cs,
    float *hh, float *xx, float *yy, float *zz, int nc, int nw)
{
    float   *hk, *yk, *zk;
    int      i, j, k, ns, nt, nk;

    nk = nw / cs;
    nt = cs * 2;
    ns = nt * 2;
    // loop over channels
    for (k = 0; k < nc; k++) {
        fzero(xx, nt);
        fcopy(xx, x, cs);
        rc(xx, nt);
        fft(xx, nt);
        // loop over sub-window segments
        yk = y + k * cs * 2;
        zk = zz + k * (nw + cs) * 2;
        for (j = 0; j < nk; j++) {
            hk = hh + (k * nk + j) * ns;
            cmul(yy, xx, hk, nt);
            ifft(yy, nt);
            for (i = 0; i < nt * 2; i++) {
                zk[i + j * cs * 2] += yy[i];
            }
        }
        fcopy(yk, zk, cs * 2);
        fmove(zk, zk + cs * 2, nw * 2);
        fzero(zk + nw * 2, cs * 2);
    }
}

// FIR-filterbank analysis for long chunk (cs >= nw)
static __inline void
firfb_analyze_lc(float *x, float *y, int cs,
    float *hh, float *xx, float *yy, float *zz, int nc, int nw)
{
    float   *hk, *yk, *zk;
    int      i, j, k, nt, ns, ni;

    nt = nw * 2;
    ns = nt * 2;
    // loop over sub-chunk segments
    for (j = 0; j < cs; j += nw) {
        ni = ((cs - j) < nw) ? (cs - j) : nw;
        fzero(xx, nt);
        fcopy(xx, x + j, ni);
        rc(xx, nt);
        fft(xx, nt);
       // loop over channels
        for (k = 0; k < nc; k++) {
            hk = hh + k * ns;
            cmul(yy, xx, hk, nt);
            ifft(yy, nt);
            yk = y + k * cs * 2;
            zk = zz + k * nw * 2;
            for (i = 0; i < ni * 2; i++) {
                yk[i + j] = yy[i] + zk[i];
            }
            fcopy(zk, yy + ni * 2, nw * 2);
        }
    }
}

/***********************************************************/

// FIR-filterbank analysis
FUNC(void)
cha_cfirfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
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
cha_cfirfb_synthesize(CHA_PTR cp, float *x, float *y, int cs)
{
    float sm;
    int i, k, kr, nc;

    nc = CHA_IVAR[_nc];
    for (i = 0; i < cs; i++) {
        sm = 0;
        /* loop over filterbank channel */
        for (k = 0; k < nc; k++) {
            kr = (i + k * cs) * 2;
            sm += x[kr];
        }
        y[i] = sm;
    }
}
