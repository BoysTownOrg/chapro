// feedback_prepare.c - feedback-management preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

// ITE feedback path from Kosgen (2006)

static float ite_fbp[100] = {
     0.000000,-0.001715, 0.000306, 0.007936,-0.014126, 0.001207, 0.001541, 0.040498, 0.077863, 0.069577,
    -0.008025,-0.106044,-0.151131,-0.124264,-0.055777, 0.013607, 0.075156, 0.082472, 0.048781, 0.004444,
    -0.033910,-0.033370,-0.013614, 0.011497, 0.031987, 0.040751, 0.031424, 0.014976,-0.001471,-0.006256,
    -0.008080,-0.003218,-0.000933, 0.001967,-0.000208,-0.003719,-0.009216,-0.012280,-0.014076,-0.013594,
    -0.013487,-0.010067,-0.007326,-0.004848,-0.003921,-0.003026,-0.002302,-0.001824,-0.001345,-0.000225,
     0.001573, 0.003371, 0.003689, 0.003543, 0.002901, 0.002056, 0.001375, 0.000886, 0.000398,-0.000091,
    -0.000580,-0.001069,-0.001569,-0.002153,-0.002737,-0.003321,-0.003905,-0.004489,-0.005073,-0.004653,
    -0.004184,-0.003715,-0.003246,-0.002778,-0.002309,-0.001965,-0.001722,-0.001479,-0.001236,-0.000993,
    -0.000750,-0.000819,-0.000969,-0.001118,-0.001268,-0.001418,-0.001567,-0.001717,-0.001867,-0.002016,
	-0.002166,-0.002200,-0.002200,-0.002200,-0.002200,-0.002200,-0.002200,-0.002200,-0.002200,-0.002200};

/***********************************************************/

FUNC(int)
cha_afc_prepare(CHA_PTR cp, double mu, double rho, int afl, double fbg, int sqm)
{
    double fbm = 0;
	float *sfbp, *efbp, *merr, gn;
	int i, cs, fbl = 0, nqm = 0, rsz = 32;

    cha_prepare(cp);
    // allocate ring buffer
	while (rsz < afl) rsz *= 2;
    cha_allocate(cp, rsz, sizeof(float), _ring);
    CHA_IVAR[_rsz] = rsz;
    CHA_IVAR[_rhd] = 0;
    CHA_IVAR[_rtl] = 0;
    // allocate afc-filter buffers
    if (afl > 0) {
        cha_allocate(cp, afl, sizeof(float), _efbp);
        efbp = (float *) cp[_efbp];
        fzero(efbp, afl);
    }
    CHA_DVAR[_mu]  = mu;
    CHA_DVAR[_rho] = rho;
    CHA_IVAR[_afl] = afl;
    // initialize simulated feedback path
	if (fbg > 0) {
        fbl = 100;
        cha_allocate(cp, fbl, sizeof(float), _sfbp);
        sfbp = (float *) cp[_sfbp];
        gn = (float) fbg;
        fbm = 0;
        for (i = 0; i < fbl; i++) {
            sfbp[i] = ite_fbp[i] * gn;
            fbm += sfbp[i] * sfbp[i];
        }
	}
    CHA_IVAR[_fbl] = fbl;
    CHA_DVAR[_fbm] = fbm;
    // initialize quality metrics
	if (sqm && (afl > 0) && (fbl > 0)) {
        nqm = (fbl < afl) ? fbl : afl;
        cs = CHA_IVAR[_cs];
        cha_allocate(cp,  cs, sizeof(float), _merr);
        merr = (float *) cp[_merr];
        fzero(merr, cs);
	}
    CHA_IVAR[_nqm] = nqm;
    // initialize instantaneous-power estimate
    CHA_DVAR[_pwr] = 0;
    return (0);
}
