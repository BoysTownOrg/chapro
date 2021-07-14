// sha_prepare.c - nonlinear-frequency-compression preparation functions

#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

static void
sha_window(float *w, int nw, int wt, int nsw)
{
    double   p, sm = 0, a = 0.16, b = 0.46;
    int      j;

    // window
    for (j = 0; j < nw; j++) {
        p = M_PI * (2.0 * j - nw) / nw;
        if (wt == 0) {
            w[j] = (float)((1 - b) + b * cos(p));                  // Hamming
        } else {
            w[j] = (float)((1 - a + cos(p) + a * cos(2 * p)) / 2); // Blackman
        }
        sm += w[j];
    }
    for (j = 0; j < nw; j++) {
        w[j] *= (float)(nw / sm / nsw);
    }
}

/***********************************************************/

FUNC(int)
cha_sha_prepare(CHA_PTR cp, CHA_SHA *sha)
{
    double  sr, Gmax, Lmax, Lckp, Lekp;
    double  g0, a1, a2, a3, gg, aa;
    float  *ww, *g1, *supp;
    int     cs, nf, nw, wt, xr, hbw;

    cs = sha->cs;
    nw = sha->nw;
    wt = sha->wt;
    xr = sha->xr;
    sr = sha->sr;
    g1 = sha->g1;
    hbw = sha->hbw;
    Gmax = sha->Gmax;
    Lmax = sha->Lmax;
    Lckp = sha->Lckp;
    Lekp = sha->Lekp;
    printf("SHA parameters: nw=%d hbw=%d ", nw, hbw);
    printf("Gmax=%.0f Lmax=%.0f Lckp=%.0f Lekp=%.0f\n", Gmax, Lmax, Lckp, Lekp);
    if (cs <= 0) {
        return (1);
    }
    if (cs % 2 != 0 || nw % 2 != 0)
        return 1;
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = sr / 1000;
    // copy SHA parameters
    CHA_IVAR[_sha_nw] = nw;
    CHA_IVAR[_sha_wt] = wt;
    // compute chunks per shift
    CHA_IVAR[_sha_ics] = 0;
    CHA_IVAR[_sha_ncs] = (nw / 2) / cs;
    // compute window
    ww = cha_allocate(cp, nw, sizeof(float), _sha_ww);
    sha_window(ww, nw, wt, 2);
    // allocate SHA buffers
    nf = nw + 1;
    cha_allocate(cp, nw, sizeof(float), _sha_xx);
    cha_allocate(cp, nf * 2, sizeof(float), _sha_yy);
    cha_allocate(cp, nf * 2, sizeof(float), _sha_XX);
    cha_allocate(cp, nf * 2, sizeof(float), _sha_YY);
    // allocate gain buffers
    if (sha->g1) {
        g1 = cha_allocate(cp, nw, sizeof(float), _sha_g1);
        fcopy(g1, sha->g1, nw);
    }
    // prepare compression parameters
    g0 = pow(10, Gmax / 20);
    a1 = pow(10, -Lckp / 20);
    a2 = pow(10,(Gmax - Lmax) / 10);
    aa = (1e12 * a2)/(1+a1 * 1e6+a2 * 1e12);
    a1 = a1 * aa;
    a2 = a2 * aa;
    a3 = (Lckp <= 0) ? 0 : pow(10, xr * Lekp / 10);
    gg = sha->ref * (nw / 4) * sqrt(2); // reference amplitude
    gg = 1 / (gg*gg);                  // intensity reference
    CHA_DVAR[_sha_g0] = g0;
    CHA_DVAR[_sha_a1] = a1;
    CHA_DVAR[_sha_a2] = a2;
    CHA_DVAR[_sha_a3] = a3;
    CHA_DVAR[_sha_gg] = gg;
    CHA_IVAR[_sha_xr] = xr;
    CHA_IVAR[_sha_hbw] = hbw;
    // allocate compression buffers
    cha_allocate(cp, nf, sizeof(float), _sha_AA);
    cha_allocate(cp, nf, sizeof(float), _sha_II);
    // copy SHA suppression
    if (sha->supp) { // copy suppressive-influence matrix
        nf = nw + 1;
        supp = (float *)cha_allocate(cp, nf * nf, sizeof(float), _sha_SS);
        fcopy(supp, sha->supp, nf * nf);
    }

    return (0);
}
