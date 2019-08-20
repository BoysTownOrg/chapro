// compressor_prepare.c - instantaneous-compression preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_gf.h"

/***********************************************************/

static void
compr_ds_expand (
    CHA_PTR cp, 
    CHA_CLS *cls,
    int ds
) {
    double fs, *bw;
    int *dsm, *dso;
    int *mf, *gc, *sc;
    int i, j, k, jk, nc, nf, ns, tc, jmx = 0, jmn = 0, scmx = 0, scmn = 0;
    static int ms = 5; /* mimimun samples per envelope cycle */

    ns = ds;
    fs = CHA_DVAR[_fs];
    nc = CHA_IVAR[_nc];
    bw = cls->bw;
    dsm = (int *) cha_allocate(cp, nc, sizeof(int), _dsm);    
    dso = (int *) cha_allocate(cp, nc, sizeof(int), _dso);
    gc = (int *) calloc(nc * ns, sizeof(int));    
    sc = (int *) calloc(ns, sizeof(int));    
    mf = (int *) calloc(nc, sizeof(int));
    nf = 0;
    for (k = ds; k  > 0; k--) {
        if ((ds % k) == 0) {
            mf[nf++] = k;
        }
    }
    for (k = 0; k < nc; k++) {
        dsm[k] = (int) floor(2000 * fs / bw[k] / ms);
        if (dsm[k] > ds) {
            dsm[k] = ds;
        }
    }
    for (k = 0; k < nc; k++) {
        for (j = 0; j < (nf-1); j++) {
            if ((dsm[k] < mf[j]) && (dsm[k] > mf[j+1])) {
                dsm[k] = mf[j+1];
            }
        }
    }
    for (k = 0; k < nc; k++) {
        dso[k] = k;
    }
    for (i = 0; i < ns; i++) {
        for (j = 0; j < ns; j++) {
            sc[j] = 0;
            for (k = 0; k < nc; k++) {
                jk = j * nc + k;
                tc = (((j - dso[k]) % dsm[k]) == 0);
                gc[jk] = tc;
                sc[j] += tc;
            }
            if (scmx < sc[j] || j == 0) {
                scmx = sc[j];
                jmx = j;
            }
            if (scmn > sc[j] || j == 0) {
                scmn = sc[j];
                jmn = j;
            }
        }
        if ((sc[jmx] - sc[jmn]) <= 1) break;
        for (k = 0; k < nc; k++) {
            jk = jmx * nc + k;
            if (gc[jk]) {
                dso[k] = ns + jmn;
                break;
            }
        }
    }
    for (k = 0; k < nc; k++) {
        dso[k] %= dsm[k];
    }
    free(gc);
    free(sc);
    free(mf);
}

/***********************************************************/

static void
compr_levl_config(CHA_PTR cp, CHA_CLS *cls, double lr)
{
    int k, nc;
    float *Lcs, *Lcm, *Lce, *Lmx;

    /* allocate arrays */
    nc = cls->nc;
    Lcs = (float *) cha_allocate(cp, nc, sizeof(float), _Lcs);
    Lcm = (float *) cha_allocate(cp, nc, sizeof(float), _Lcm);
    Lce = (float *) cha_allocate(cp, nc, sizeof(float), _Lce);
    Lmx = (float *) cha_allocate(cp, nc, sizeof(float), _Lmx);
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++) {
        Lcs[k] = (float) cls->Lcs[k];
        Lcm[k] = (float) cls->Lcm[k];
        Lce[k] = (float) cls->Lce[k];
        Lmx[k] = (float) cls->Lmx[k];
    }
    CHA_DVAR[_lrpk] = lr * M_SQRT2;
}

static void
compr_gain_config(CHA_PTR cp, CHA_CLS *cls)
{
    float *Gcs, *Gcm, *Gce, *Gmx;
    int k, nc;

    CHA_IVAR[_cm] = 1; // set compression mode
    /* allocate arrays */
    nc = cls->nc;
    Gcs = (float *) cha_allocate(cp, nc, sizeof(float), _Gcs);
    Gcm = (float *) cha_allocate(cp, nc, sizeof(float), _Gcm);
    Gce = (float *) cha_allocate(cp, nc, sizeof(float), _Gce);
    Gmx = (float *) cha_allocate(cp, nc, sizeof(float), _Gmx);
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++) {
        Gcs[k] = (float) cls->Gcs[k];
        Gcm[k] = (float) cls->Gcm[k];
        Gce[k] = (float) cls->Gce[k];
        Gmx[k] = (float) cls->Gmx[k];
    }
}

static void
compressor_prep(CHA_PTR cp, int nc)
{
    float *Gmn, *Gcs, *Gce;
    int *dsm, k;

    /* allocate arrays */
    dsm = (int *) cp[_dsm];
    /* initialize Gmn */
    Gmn = (float *) cha_allocate(cp, nc, sizeof(float), _Gmn);
    Gcs = (float *) cp[_Gcs];
    Gce = (float *) cp[_Gce];
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++) {
       Gmn[k] = fmin(Gcs[k], Gce[k]);
    }
    cha_allocate(cp, nc, sizeof(float), _Lfb);
    cha_allocate(cp, nc, sizeof(float), _Gsup);
    cha_allocate(cp, nc, sizeof(float), _Gpre);
    cha_allocate(cp, nc, sizeof(float), _gsup);
    cha_allocate(cp, nc, sizeof(float), _ginc);
    cha_allocate(cp, dsm[0] * nc * 2, sizeof(float), _zdr);
}

/***********************************************************/

FUNC(int)
cha_icmp_prepare(CHA_PTR cp, CHA_CLS *cls, double lr, int ds)
{
    cha_prepare(cp);
    compr_ds_expand(cp, cls, ds);
    compr_levl_config(cp, cls, lr);
    compr_gain_config(cp, cls);
    compressor_prep(cp, cls->nc);

    return (0);
}

