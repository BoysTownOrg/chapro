// feedback_process.c - feedback-management processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

FUNC(void)
cha_afc_input(CHA_PTR cp, float *x, float *y, int cs)
{
    float *rng0, *rng1, *rng2, *rng3, *efbp, *sfbp, *wfrp, *ffrp, *merr;
    float ye, yy, mum, dif, fbm, dm, xx, ss, ee, uu, ef, uf, mu, rho, eps, pwr;
    int i, ii, ij, j, afl, wfl, ffl, fbl, nqm, rsz, rhd, mask;

    mu  = (float) CHA_DVAR[_mu];
    rho = (float) CHA_DVAR[_rho];
    eps = (float) CHA_DVAR[_eps];
    pwr = (float) CHA_DVAR[_pwr];
    fbm = (float) CHA_DVAR[_fbm];
    rsz = CHA_IVAR[_rsz];
    rhd = CHA_IVAR[_rhd];
    afl = CHA_IVAR[_afl];
    wfl = CHA_IVAR[_wfl];
    ffl = CHA_IVAR[_ffl];
    fbl = CHA_IVAR[_fbl];
    nqm = CHA_IVAR[_nqm];
    efbp = (float *) cp[_efbp];
    sfbp = (float *) cp[_sfbp];
    wfrp = (float *) cp[_wfrp];
    ffrp = (float *) cp[_ffrp];
    merr = (float *) cp[_merr];
    rng0 = (float *) cp[_rng0];
    rng1 = (float *) cp[_rng1];
    rng2 = (float *) cp[_rng2];
    rng3 = (float *) cp[_rng3];
    if (ffl <= 0) rng3 = rng0; // bypass rng3
    if (wfl <= 0) rng2 = rng3; // bypass rng2
    mask = rsz - 1;
    // subtract estimated feedback signal
    for (i = 0; i < cs; i++) {
        xx = x[i];
        ii = (rhd + i) & mask;
        // simulate feedback
        yy = 0;
        for (j = 0; j < fbl; j++) {
            ij = (ii - j + rsz) & mask;
            ss = rng0[ij];
            yy += sfbp[j] * ss;
        }
        // apply fixed-feedback filter
        if (ffl > 0) {
            uu = 0;
            for (j = 0; j < ffl; j++) {
                ij = (ii - j + rsz) & mask;
                ss = rng0[ij];
                uu += ffrp[j] * ss;
            }
            rng3[ii] = uu;
        }
        // estimate feedback
        ye = 0;
        for (j = 0; j < afl; j++) {
            ij = (ii - j + rsz) & mask;
            uu = rng3[ij];
            ye += efbp[j] * uu;
        }
        // apply feedback to input signal
        ee = xx + yy - ye;
        // apply signal-whitening filter
        if (wfl > 0) {
            rng1[ii] = ee;
            ef = uf = 0;
            for (j = 0; j < wfl; j++) {
                ij = (ii - j + rsz) & mask;
                ef += rng1[ij] * wfrp[j];
                uf += rng3[ij] * wfrp[j];
            }
            rng2[ii] = uf;
        } else {
            ef = ee;
        }
        // update adaptive feedback coefficients
        ss = rng0[ii];
        pwr = rho * pwr + ee * ee + ss * ss;
        mum = mu / (eps + pwr);  // modified mu
        for (j = 0; j < afl; j++) {
            ij = (ii - j + rsz) & mask;
            uf = rng2[ij];
            efbp[j] += mum * ef * uf;
        }
        // save quality metrics
        if (nqm > 0) {
            dm = 0;
            for (j = 0; j < nqm; j++) {
                dif = sfbp[j] - efbp[j];
                dm += dif * dif;
            }
            merr[i] = dm / fbm;
        }
        // copy AFC signal to output
        y[i] = ee;
    }
    // save power estimate
    CHA_DVAR[_pwr] = pwr;
}

FUNC(void)
cha_afc_output(CHA_PTR cp, float *x, int cs)
{
    float *ring;
    int i, j, rsz, rhd, rtl, mask;

    rsz = CHA_IVAR[_rsz];
    rhd = CHA_IVAR[_rhd];
    rtl = CHA_IVAR[_rtl];
    ring = (float *) cp[_rng0];
    mask = rsz - 1;
    rhd = rtl;
    // copy chunk to ring buffer
    for (i = 0; i < cs; i++) {
        j = (rhd + i) & mask;
        ring[j] = x[i];
    }
    rtl = (rhd + cs) % rsz;
    CHA_IVAR[_rhd] = rhd;
    CHA_IVAR[_rtl] = rtl;
}
