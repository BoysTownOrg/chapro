// feedback_process.c - feedback-management processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

static int rhd = 0;

/***********************************************************/

FUNC(void)
cha_afc_input(CHA_PTR cp, float *x, float *y, int cs)
{
    float ye, yy, mum, dif, dm, xx, ss, ee, uu, ef, uf;
    int i, ii, ij, j;
    static float *rng0, *rng1, *rng2, *rng3, *efbp, *sfbp, *wfrp, *ffrp, *merr;
    static float mu, rho, eps, fbm;
    static float pwr = 0;
    static int rsz, mask, afl, wfl, ffl, fbl, nqm; 
    static int first_time = 1;

    if (first_time) {
        efbp = (float *) cp[_efbp];
        sfbp = (float *) cp[_sfbp];
        wfrp = (float *) cp[_wfrp];
        ffrp = (float *) cp[_ffrp];
        merr = (float *) cp[_merr];
        rng0 = (float *) cp[_rng0];
        rng1 = (float *) cp[_rng1];
        rng2 = (float *) cp[_rng2];
        rng3 = (float *) cp[_rng3];
        mu  = (float) CHA_DVAR[_mu];
        rho = (float) CHA_DVAR[_rho];
        eps = (float) CHA_DVAR[_eps];
        fbm = (float) CHA_DVAR[_fbm];
        rsz = CHA_IVAR[_rsz];
        afl = CHA_IVAR[_afl];
        wfl = CHA_IVAR[_wfl];
        ffl = CHA_IVAR[_ffl];
        fbl = CHA_IVAR[_fbl];
        nqm = CHA_IVAR[_nqm];
        if (ffl <= 0) rng3 = rng0; // bypass rng3
        if (wfl <= 0) rng2 = rng3; // bypass rng2
        mask = rsz - 1;
        first_time = 0;
    }
    // subtract estimated feedback signal
    for (i = 0; i < cs; i++) {
        xx = x[i];
        ii = (rhd + i) & mask;
        // simulate feedback
        yy = 0;
        for (j = 0; j < fbl; j++) {
            ij = (ii - j + rsz) & mask;
            yy += sfbp[j] * rng0[ij];
        }
        // apply fixed-feedback filter
        if (ffl > 0) {
            uu = 0;
            for (j = 0; j < ffl; j++) {
                ij = (ii - j + rsz) & mask;
                uu += ffrp[j] * rng0[ij];
            }
            rng3[ii] = uu;
        }
        // estimate feedback
        if (afl > 0) {
            ye = 0;
            for (j = 0; j < afl; j++) {
                ij = (ii - j + rsz) & mask;
                ye += efbp[j] * rng3[ij];
            }
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
}

FUNC(void)
cha_afc_output(CHA_PTR cp, float *x, int cs)
{
    int i, j;
    static float *rng0;
    static int rsz, mask;
    static int rtl = 0;
    static int first_time = 1;

    if (first_time) {
        rng0 = (float *) cp[_rng0];
        rsz = CHA_IVAR[_rsz];
        mask = rsz - 1;
        first_time = 0;
    }
    // copy chunk to ring buffer
    rhd = rtl;
    for (i = 0; i < cs; i++) {
        j = (rhd + i) & mask;
        rng0[j] = x[i];
    }
    rtl = (rhd + cs) % rsz;
}
