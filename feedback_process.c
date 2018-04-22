// feedback_process.c - feedback-management processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

#define db1(x)          (10*log10f(x))

/***********************************************************/

FUNC(void)
cha_afc_input(CHA_PTR cp, float *x, float *y, int cs)
{
	double mu, rho, pwr, eps = 1e-9;
	float *ring, *filt, *sfbp, *merr, fbe, fbs, mum, dif, fbp, dm, fm, s0, s1;
	int i, ii, ij, j, afl, fbl, nqm, rsz, rhd;

    mu  = CHA_DVAR[_mu];
    rho = CHA_DVAR[_rho];
    pwr = CHA_DVAR[_pwr];
    rsz = CHA_IVAR[_rsz];
    rhd = CHA_IVAR[_rhd];
    afl = CHA_IVAR[_afl];
    fbl = CHA_IVAR[_fbl];
    nqm = CHA_IVAR[_nqm];
    ring = (float *) cp[_ring];
    filt = (float *) cp[_filt];
    sfbp = (float *) cp[_sfbp];
    merr = (float *) cp[_merr];
	// subtract estimated feedback signal
	for (i = 0; i < cs; i++) {
        s0 = x[i];
		ii = rhd + i;
		// simulate feedback
        fbs = 0;
		for (j = 0; j < fbl; j++) {
            ij = ii - j;
            if (ij >= rsz) {
                ij -= rsz;
			} else if (ij < 0) {
				ij += rsz;
			}
            fbs += ring[ij] * sfbp[j];
		}
		// estimate feedback
        fbe = 0;
		for (j = 0; j < afl; j++) {
            ij = ii - j;
            if (ij >= rsz) {
                ij -= rsz;
			} else if (ij < 0) {
				ij += rsz;
			}
            fbe += ring[ij] * filt[j];
		}
		// apply feedback to input signal
		s1 = s0 + fbs - fbe;
		// update adaptive feedback coefficients
		pwr = rho * pwr + s0 * s0 + s1 * s1;
        mum = (float) (mu / (eps + pwr));  // modified mu
		for (j = 0; j < afl; j++) {
            ij = ii - j;
            if (ij >= rsz) {
                ij -= rsz;
			} else if (ij < 0) {
				ij += rsz;
			}
            filt[j] += mum * ring[ij] * s1;
		}
		// save quality metrics
		if (nqm > 0) {
            dm = fm = 0;
            for (j = 0; j < nqm; j++) {
                dif = sfbp[j] - filt[j];
                fbp = sfbp[j];
                dm += dif * dif;
                fm += fbp * fbp;
            }
            merr[i] = db1(dm / fm);
        }
        // copy AFC signal to output
        y[i] = s1;
	}
    // save power estimate
    CHA_DVAR[_pwr] = pwr;
}

FUNC(void)
cha_afc_output(CHA_PTR cp, float *x, int cs)
{
	float *ring;
	int i, j, rsz, rhd, rtl;

    rsz = CHA_IVAR[_rsz];
    rhd = CHA_IVAR[_rhd];
    rtl = CHA_IVAR[_rtl];
    ring = (float *) cp[_ring];
	rhd = rtl;
	// copy chunk to ring buffer
	for (i = 0; i < cs; i++) {
		j = rhd + i;
		if (j >= rsz) {
			j -= rsz;
		}
		ring[j] = x[i];
	}
	rtl = (rhd + cs) % rsz;
    CHA_IVAR[_rhd] = rhd;
    CHA_IVAR[_rtl] = rtl;
}
