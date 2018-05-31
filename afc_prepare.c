// feedback_prepare.c - feedback-management preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

// ITE feedback path from Kosgen (2006)

static double ite_fbp[100] = {
     0.001764, 0.000049, 0.002070, 0.009700,-0.012362, 0.002971, 0.003305, 0.042262, 0.079627, 0.071341,
    -0.006261,-0.104280,-0.149367,-0.122500,-0.054013, 0.015371, 0.076920, 0.084236, 0.050545, 0.006208,
    -0.032146,-0.031606,-0.011850, 0.013261, 0.033751, 0.042515, 0.033188, 0.016740, 0.000293,-0.004492,
    -0.006316,-0.001454, 0.000831, 0.003731, 0.001556,-0.001955,-0.007452,-0.010516,-0.012312,-0.011830,
    -0.011723,-0.008303,-0.005562,-0.003084,-0.002157,-0.001262,-0.000538,-0.000060, 0.000419, 0.001539,
     0.003337, 0.005135, 0.005453, 0.005307, 0.004665, 0.003820, 0.003139, 0.002650, 0.002162, 0.001673,
     0.001184, 0.000695, 0.000195,-0.000389,-0.000973,-0.001557,-0.002141,-0.002725,-0.003309,-0.002889,
    -0.002420,-0.001951,-0.001482,-0.001014,-0.000545,-0.000201, 0.000042, 0.000285, 0.000528, 0.000771,
     0.001014, 0.000945, 0.000795, 0.000646, 0.000496, 0.000346, 0.000197, 0.000047,-0.000103,-0.000252,
    -0.000402,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436};

/***********************************************************/

FUNC(int)
cha_afc_prepare(CHA_PTR cp, double mu, double rho, double eps, int afl, 
                int wfl, int ffl, double fbg, int sqm)
{
    double fbm = 0;
    float *sfbp;
    int i, cs, fbl = 0, nqm = 0, rsz = 32;

    cha_prepare(cp);
    // allocate HA-output ring buffer
    while (rsz < afl) rsz *= 2;
    while (rsz < wfl) rsz *= 2;
    while (rsz < ffl) rsz *= 2;
    cha_allocate(cp, rsz, sizeof(float), _rng0);
    CHA_IVAR[_rsz] = rsz;
    // allocate AFC-filter buffer
    if (afl > 0) {
        cha_allocate(cp, afl, sizeof(float), _efbp);
    }
    CHA_DVAR[_mu]  = mu;
    CHA_DVAR[_rho] = rho;
    CHA_DVAR[_eps] = eps;
    CHA_IVAR[_afl] = afl;
    // initialize signal-whitening filter
    if (wfl > 0) {
        cha_allocate(cp, wfl, sizeof(float), _wfrp);
        cha_allocate(cp, rsz, sizeof(float), _rng1);
        cha_allocate(cp, rsz, sizeof(float), _rng2);
    }
    CHA_IVAR[_wfl] = wfl;
    // initialize fixed-feedback filter
    if (ffl > 0) {
        cha_allocate(cp, ffl, sizeof(float), _ffrp);
        cha_allocate(cp, rsz, sizeof(float), _rng3);
    }
    CHA_IVAR[_ffl] = ffl;
    // initialize simulated feedback path
    if (fbg > 0) {
        fbl = 100;
        cha_allocate(cp, fbl, sizeof(float), _sfbp);
        sfbp = (float *) cp[_sfbp];
        fbm = 0;
        for (i = 0; i < fbl; i++) {
            sfbp[i] = (float) (ite_fbp[i] * fbg);
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
    }
    CHA_IVAR[_nqm] = nqm;
    return (0);
}
