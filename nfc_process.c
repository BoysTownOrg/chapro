// nfc_process.c - nonlinear-frequency-compression processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"

#include "seriallog.h"

/***********************************************************/
#ifdef ARM_DSP

static int fft_initialized = 0;
static arm_rfft_fast_instance_f32 fft_instance;

static __inline void
rfft(float *x, int n)
{
    if (!fft_initialized)
    {
        arm_rfft_fast_init_f32(&fft_instance, n);
        fft_initialized++;
    }
    arm_rfft_fast_f32(&fft_instance, x, x, 0);
}

static __inline void
rifft(float *x, int n)
{
    if (!fft_initialized)
    {
        arm_rfft_fast_init_f32(&fft_instance, n);
        fft_initialized++;
    }
    arm_rfft_fast_f32(&fft_instance, x, x, 1);
}

#else // ARM_DSP

#define rfft(x, n) cha_fft_rc(x, n)
#define rifft(x, n) cha_fft_cr(x, n)

#endif // ARM_DSP

// map frequencies
static __inline void
fmap(float *y, float *x, int nf, int *mm, int nm)
{
    float *xx, *yy;
    int k, kk, nn;

    seriallog("fmap");

    nn = mm[0] * 2;
    fcopy(y, x, nn);
    fzero(y + nn, nf - nn);

    for (k = 0; k < (nm - 1); k++)
    {
        yy = y + 2 * k + nn;
        for (kk = mm[k]; kk < mm[k + 1]; kk++)
        {
            xx = x + 2 * kk;
            yy[0] += xx[0]; // real
            yy[1] += xx[1]; // imag
        }
    }
}

// NFC analyze
static __inline void
short_term_analyze(float *xx, float *XX, float *ww, int nw, int ns, int nf)
{
    int i;

    for (i = 0; i < nw; i++)
    {
        XX[i] = xx[i] * ww[i]; // apply window to input
    }
    fzero(XX + nw, nw);
    rfft(XX, nf);           // FFT
    fcopy(xx, xx + ns, ns); // save last half of input window
}

// NFC synthesize
static __inline void
short_term_synthesize(float *yy, float *YY, int ns, int nf)
{
    int i, nn;

    rifft(YY, nf); // IFFT
    nn = nf - ns;
    fmove(yy, yy + ns, nn); // shift previous output
    for (i = 0; i < nn; i++)
    {
        yy[i] += YY[i]; // overlap-add output
    }
    fcopy(yy + nn, YY + nn, ns); // save response tail
}

// nonlinear-frequency-compression short chunk
static __inline void
nfc_sc(CHA_PTR cp, float *x, float *y, int cs,
       float *xx, float *yy, float *XX, float *YY, float *ww,
       int *mm, int nm, int nw)
{
    int icp, ics, nn, nf, ns, ncs;

    seriallog("nfc_sc");

    // process chunk
    ncs = CHA_IVAR[_nfc_ncs];
    ics = CHA_IVAR[_nfc_ics];
    ns = nw / 2;
    nf = nw * 2;
    nn = ics * cs;
    fcopy(xx + nn + ns, x, cs);
    fcopy(y, yy + nn, cs);
    icp = (ics + 1) % ncs;
    if (icp == 0)
    { // perform NFC after every shift
        short_term_analyze(xx, XX, ww, nw, ns, nf);
        fmap(YY, XX, nf, mm, nm); // compress frequency range
        short_term_synthesize(yy, YY, ns, nf);
    }
    // update chunk count
    CHA_IVAR[_nfc_ics] = icp;
}

// nonlinear-frequency-compression long chunk
static __inline void
nfc_lc(float *x, float *y, int cs,
       float *xx, float *yy, float *XX, float *YY, float *ww,
       int *mm, int nm, int nw)
{
    int k, nn, nf, ns;

    seriallog("nfc_lc");

    // process chunk
    ns = nw / 2;
    nf = nw * 2;
    nn = cs / ns;
    for (k = 0; k < nn; k++)
    {
        fcopy(xx + nn + ns, x + k * ns, ns);
        fcopy(y + k * ns, yy + nn, ns);
        // perform NFC after every shift
        short_term_analyze(xx, XX, ww, nw, ns, nf);
        fmap(YY, XX, nf, mm, nm); // compress frequency range
        short_term_synthesize(yy, YY, ns, nf);
    }
}

/***********************************************************/

// nonlinear-frequency-compression analysis
FUNC(void)
cha_nfc_process(CHA_PTR cp, float *x, float *y, int cs)
{
    float *ww, *xx, *yy, *XX, *YY;
    int nw, nm, *mm;

    seriallog("cha_nfc_process");

    // copy parameters and pointers from cha_data
    nw = CHA_IVAR[_nfc_nw];
    nm = CHA_IVAR[_nfc_nm];
    mm = (int *)cp[_nfc_mm];
    ww = (float *)cp[_nfc_ww];
    xx = (float *)cp[_nfc_xx];
    yy = (float *)cp[_nfc_yy];
    XX = (float *)cp[_nfc_XX];
    YY = (float *)cp[_nfc_YY];
    if (cs < (nw / 2))
    { // short chunk ??
        nfc_sc(cp, x, y, cs, xx, yy, XX, YY, ww, mm, nm, nw);
    }
    else
    { // long chunk (not yet implemented)
        nfc_lc(x, y, cs, xx, yy, XX, YY, ww, mm, nm, nw);
    }
}
