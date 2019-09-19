// afc_process.c - adaptive-feedback-cancelation processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

/***********************************************************/

FUNC(void)
cha_afc_input(CHA_PTR cp, float *x, float *y, int cs)
{
    float ye, yy, mmu, dif, dm, xx, ss, ee, uu, ef, uf;
    int i, ih, ij, is, id, j, n, iqm;
    static float *rng0, *rng1, *rng2, *rng3;
    static float *efbp, *sfbp, *wfrp, *ffrp, *qm;
    static float mu, rho, eps, fbm;
    static float pwr = 0;
    static int rhd, rsz, mask, afl, wfl, pfl, fbl, nqm, hdel, *iqmp; 

    if (CHA_IVAR[_mxl] == 0) { // if no AFC, do nothing
        return;
    }
    if (CHA_IVAR[_in1] == 0) {
        iqmp = (int *)   cp[_iqmp];
        qm   = (float *) cp[_qm];
        efbp = (float *) cp[_efbp];
        sfbp = (float *) cp[_sfbp];
        wfrp = (float *) cp[_wfrp];
        ffrp = (float *) cp[_ffrp];
        rng0 = (float *) cp[_rng0];
        rng3 = (float *) cp[_rng3];
        rng2 = (float *) cp[_rng2];
        rng1 = (float *) cp[_rng1];
        mu   = (float) CHA_DVAR[_mu];
        rho  = (float) CHA_DVAR[_rho];
        eps  = (float) CHA_DVAR[_eps];
        fbm  = (float) CHA_DVAR[_fbm];
        rsz  = CHA_IVAR[_rsz];
        afl  = CHA_IVAR[_afl];
        wfl  = CHA_IVAR[_wfl];
        pfl  = CHA_IVAR[_pfl];
        fbl  = CHA_IVAR[_fbl];
        nqm  = CHA_IVAR[_nqm];
        hdel = CHA_IVAR[_hdel];
        if (pfl <= 0) rng1 = rng0; // bypass rng1
        if (wfl <= 0) rng2 = rng1; // bypass rng2
        mask = rsz - 1;
    }
    // ss -> rng0
    // uu -> rng1
    // uf -> rng2
    // ee -> rng3
    rhd = CHA_IVAR[_rhd];
    iqm = iqmp[0];
    // subtract estimated feedback signal
    for (i = 0; i < cs; i++) {
        xx = x[i];
        ih = (rhd + i) & mask;
        is = ih + rsz;
        id = is - hdel;
        // simulate feedback
        yy = 0;
        for (j = 0; j < fbl; j++) {
            ij = (id - j) & mask;
            yy += sfbp[j] * rng0[ij];
        }
        // apply persistent-feedback filter
        ss = rng0[ih];
        if (pfl > 0) {
            uu = 0;
            for (j = 0; j < pfl; j++) {
                ij = (is - j) & mask;
                uu += ffrp[j] * rng0[ij];
            }
            rng1[ih] = uu;
        }
        // estimate feedback
        ye = 0;
        if (afl > 0) {
            for (j = 0; j < afl; j++) {
                ij = (id - j) & mask;
                ye += efbp[j] * rng1[ij];
            }
        }
        // apply feedback to input signal
        ee = xx + yy - ye;
        // apply signal-whitening filter
        if (wfl > 0) {
            rng3[ih] = ee;
            ef = uf = 0;
            for (j = 0; j < wfl; j++) {
                ij = (is - j) & mask;
                ef += rng3[ij] * wfrp[j];
                uf += rng1[ij] * wfrp[j];
            }
            rng2[ih] = uf;
        } else {
            ef = ee;
        }
        // update adaptive feedback coefficients
        if (afl > 0) {
            uf = rng2[id & mask];
            pwr = rho * pwr + ef * ef + uf * uf;
            mmu = mu / (eps + pwr);  // modified mu
            for (j = 0; j < afl; j++) {
                ij = (id - j) & mask;
                uf = rng2[ij];
                efbp[j] += mmu * ef * uf;
            }
        }
        // save quality metrics
        if (nqm > 0) {
            dm = 0;
            for (j = 0; j < afl; j++) {
                dif = sfbp[j] - efbp[j];
                dm += dif * dif;
            }
            qm[iqm++] = dm / fbm;
        }
        // copy AFC signal to output
        y[i] = ee;
    }
    iqmp[0] = iqm;
}

FUNC(void)
cha_afc_output(CHA_PTR cp, float *x, int cs)
{
    int i, j, rhd;
    static float *rng0;
    static int rsz, mask;
    static int rtl = 0;

    if (CHA_IVAR[_in1] == 0) {
        rng0 = (float *) cp[_rng0];
        rsz = CHA_IVAR[_rsz];
        mask = rsz - 1;
        CHA_IVAR[_in1] = CHA_IVAR[_in2];
    }
        rsz = CHA_IVAR[_rsz];
    // copy chunk to ring buffer
    rhd = rtl;
    for (i = 0; i < cs; i++) {
        j = (rhd + i) & mask;
        rng0[j] = x[i];
    }
    rtl = (rhd + cs) % rsz;
    CHA_IVAR[_rhd] = rhd;
}
