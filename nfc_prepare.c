// nfc_prepare.c - nonlinear-frequency-compression preparation functions

#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

static int
nfc_map(int nw, double lbf, double ubf, double sr, int *map)
{
    double df, dk, kk;
    int k, n1, n2, nn;

    df  = sr / (2 * nw);
    lbf = fmax( df,fmin(lbf, sr/2));
    ubf = fmax(lbf,fmin(ubf, sr/2));
    n1  = round(lbf / df);
    n2  = round(ubf / df);
    nn  = n2 - n1 + 1;
    if (map) {
        dk = log((double) nw / n1) / log((double) n2 / n1);
        for (k = 0; k < nn; k++) {
            kk = log((double) (k + n1) / n1);
            map[k] = round(n1 * exp(kk * dk));
        }
    }

    return (nn);
}


static void
nfc_window(float *w, int nw, int wt, int nsw)
{
    double   p, sm = 0, a = 0.16, b = 0.46;
    int      j;

    // window
    for (j = 0; j < nw; j++) {
        p = M_PI * (2.0 * j - nw) / nw;
        if (wt == 0) {
            w[j] = (float)((1 - b) + b * cos(p));                   // Hamming
        } else {
            w[j] = (float)((1 - a + cos(p) + a * cos(2 * p)) / 2);  // Blackman
        }
        sm += w[j];
    }
    for (j = 0; j < nw; j++) {
        w[j] *= (float)(nw / sm / nsw);
    }
}

/***********************************************************/

FUNC(int)
cha_nfc_prepare(CHA_PTR cp, CHA_NFC *nfc)
{
    double  sr, lbf, ubf;
    float  *ww, *g1, *g2;
    int      cs, nf, nw, nn, wt, nm, *mm;

    cs = nfc->cs;
    nw = nfc->nw;
    nm = nfc->nm;
    wt = nfc->wt;
    sr = nfc->sr;
    lbf = nfc->lbf;
    ubf = nfc->ubf;
    printf("NFC parameters: nw=%d lbf=%.0f ubf=%.0f nm=%d\n",
        nw, lbf, ubf, nm);
    if (cs <= 0) {
        return (1);
    }
    if (cs % 2 != 0 || nw % 2 != 0)
        return 1;
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = sr / 1000;
    // copy NFC parameters
    CHA_IVAR[_nfc_nw] = nw;
    CHA_IVAR[_nfc_wt] = wt;
    // specify NFC frequency map 
    if (nfc->mm && nfc->nm) { // copy from NFC struct ??
        nm = CHA_IVAR[_nfc_nm] = nfc->nm;
        mm = (int32_t *)cha_allocate(cp, nm, sizeof(int32_t), _nfc_mm);
        memcpy(mm,nfc->mm,nm * sizeof(int32_t));
    } else {                   // compute log-frequency map
        nm = CHA_IVAR[_nfc_nm] = nfc_map(nw, lbf, ubf, sr, 0);
        mm = (int32_t *)cha_allocate(cp, nm, sizeof(int32_t), _nfc_mm);
        nfc_map(nw, lbf, ubf, sr, mm);
    }
    // compute chunks per shift
    CHA_IVAR[_nfc_ics] = 0;
    CHA_IVAR[_nfc_ncs] = (nw / 2) / cs;
    // compute window
    ww = cha_allocate(cp, nw, sizeof(float), _nfc_ww);
    nfc_window(ww, nw, wt, 2);
    // allocate NFC buffers
    nf = nw * 2 + 1;
    cha_allocate(cp, nw, sizeof(float), _nfc_xx);
    cha_allocate(cp, nf * 2, sizeof(float), _nfc_yy);
    cha_allocate(cp, nf * 2, sizeof(float), _nfc_XX);
    cha_allocate(cp, nf * 2, sizeof(float), _nfc_YY);
    // allocate gain buffers
    if (nfc->g1) {
        g1 = cha_allocate(cp, nw, sizeof(float), _nfc_g1);
        fcopy(g1, nfc->g1, nw);
    }
    if (nfc->g2) {
        nn = round(2 * ubf / sr);
        g2 = cha_allocate(cp, nn, sizeof(float), _nfc_g2);
        fcopy(g2, nfc->g2, nn);
    }

    return (0);
}
