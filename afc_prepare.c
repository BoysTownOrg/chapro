// afc_prepare.c - adaptive-feedback-cancelation preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "ite_fb.h"
#include "iltass.h"

static float bndlmt[100] = {
     0.92015965,-0.15211262,-0.14604275,-0.13782281,-0.12793407,-0.11682282,-0.10489615,
    -0.09251921,-0.08001370,-0.06765745,-0.05568504,-0.04428918,-0.03362277,-0.02380146,
    -0.01490669,-0.00698888,-0.00007095, 0.00584826, 0.01079044, 0.01479437, 0.01791271,
     0.02020885, 0.02175406, 0.02262487, 0.02290074, 0.02266197, 0.02198791, 0.02095546,
     0.01963786, 0.01810364, 0.01641592, 0.01463189, 0.01280246, 0.01097212, 0.00917895,
     0.00745474, 0.00582527, 0.00431062, 0.00292560, 0.00168017, 0.00057997,-0.00037318,
    -0.00118079,-0.00184716,-0.00237889,-0.00278442,-0.00307360,-0.00325725,-0.00334682,
    -0.00335406,-0.00329071,-0.00316831,-0.00299792,-0.00279001,-0.00255432,-0.00229976,
    -0.00203432,-0.00176508,-0.00149818,-0.00123881,-0.00099126,-0.00075895,-0.00054449,
    -0.00034975,-0.00017591,-0.00002356, 0.00010725, 0.00021692, 0.00030624, 0.00037630,
     0.00042843, 0.00046414, 0.00048506, 0.00049291, 0.00048940, 0.00047625, 0.00045513,
     0.00042763, 0.00039524, 0.00035934, 0.00032117, 0.00028187, 0.00024241, 0.00020363,
     0.00016623, 0.00013079, 0.00009777, 0.00006749, 0.00004019, 0.00001601,-0.00000501,
    -0.00002289,-0.00003771,-0.00004960,-0.00005875,-0.00006535,-0.00006965,-0.00007188,
    -0.00007230,-0.00007116};

/***********************************************************/

static void
white_filt(float *h, int n)
{
    int i, m;

    if (n < 3) {
        h[0] = 1;
    } else {
        m = (n - 1) / 2;
        if (m >= WFSZ) m = WFSZ - 1;
        h[m] = iltass[0];
        for (i = 1; i <= m; i++) {
            h[m - i] = h[m + i] = iltass[i];
        }
    }
}

/***********************************************************/

FUNC(int)
cha_afc_prepare(CHA_PTR cp, CHA_AFC *afc)
{
    double mu, rho, eps, alf, fbg;
    int afl, wfl, pfl, pup, hdel, sqm;
    double fbm = 0;
    float *sfbp, *wfrp, *ffrp;
    int i, cs, fbl, mxl = 0, nqm, rsz = 32;

    cha_prepare(cp);
    mu  = afc->mu;
    rho = afc->rho;
    eps = afc->eps;
    alf = afc->alf;
    fbg = afc->fbg;
    afl = afc->afl;
    wfl = afc->wfl;
    pfl = afc->pfl;
    pup = afc->pup;
    sqm = afc->sqm;
    nqm = afc->nqm;
    hdel = afc->hdel;
    // allocate HA-output ring buffer
    cs = CHA_IVAR[_cs];
    fbl = (fbg > 0) ? FBSZ : 0;
    mxl = (fbl > afl) ? fbl : (afl > wfl) ? afl : (wfl > pfl) ? wfl : pfl;
    while (rsz < (mxl + hdel + cs))  rsz *= 2;
    CHA_IVAR[_rsz] = rsz;
    CHA_IVAR[_mxl] = mxl;
    if (mxl == 0) { // if no AFC, do nothing
        return (0);
    }
    CHA_IVAR[_in1] = 0;
    CHA_IVAR[_in2] = 1;
    CHA_IVAR[_rhd] = 0;
    cha_allocate(cp, rsz, sizeof(float), _rng0);
    // allocate AFC-filter buffer
    if (afl > 0) {
        cha_allocate(cp, afl, sizeof(float), _efbp);
    }
    CHA_IVAR[_afl] = afl;
    CHA_DVAR[_mu]  = mu;
    CHA_DVAR[_rho] = rho;
    CHA_DVAR[_eps] = eps;
    CHA_IVAR[_puc] = 0;
    // initialize whiten filter
    if (wfl > 0) {
        cha_allocate(cp, rsz, sizeof(float), _rng3); // ee -> rng3
        cha_allocate(cp, rsz, sizeof(float), _rng2); // uf -> rng2
        wfrp = cha_allocate(cp, wfl, sizeof(float), _wfrp);
        white_filt(wfrp, wfl);
    }
    CHA_IVAR[_wfl] = wfl;
    // initialize band-limit filter
    if (pfl > 0) {
        cha_allocate(cp, rsz, sizeof(float), _rng1); // uu -> rng1
        ffrp = cha_allocate(cp, pfl, sizeof(float), _ffrp);
        for (i = 0; i < pfl; i++) {
            ffrp[i] = (float)alf * bndlmt[i];
        }
    } else {
        pup = 0;
    }
    CHA_IVAR[_pfl] = pfl;
    CHA_IVAR[_pup] = pup;
    CHA_DVAR[_alf] = alf;
    // initialize simulated feedback path
    if (fbl > 0) {
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
    if (sqm && (afl > 0) && (fbl >= afl)) {
        cha_allocate(cp, nqm, sizeof(float), _qm);
        cha_allocate(cp, 1, sizeof(int), _iqmp);
        CHA_IVAR[_nqm] = nqm;
    } else {
        CHA_IVAR[_iqmp] = 0;
        CHA_IVAR[_nqm] = 0;
    }
    // initialize hardware delay
    CHA_IVAR[_hdel] = hdel; // should this be 38 ???
    // copy filter info back to CHA_AFC
    afc->pcp = NULL;
    cha_afc_filters(cp, afc);

    return (0);
}

