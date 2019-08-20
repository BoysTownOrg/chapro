// compressor_process.c - instantaneous-compression processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_gf.h"

/***********************************************************/

FUNC(void)
cha_icmp_process(CHA_PTR cp, float *x, float *y, int cs)
{
    float *Lcs, *Lcm, *Lce, *Lmx, *Gcs, *Gcm, *Gce, *Gmx;
    float *Gmn, *Lfb, *Gsup, *Gpre, *gsup, *ginc, *zdr, *zdi;
    float Lsup = 0, rnge, head;
    float agtf, lrpk;
    int *dsm, *dso;
    int i, k, kr, ki, nc, dsmx, id, kd, dsmo, dsi, im, cm, jm;
    static float eps = 1e-20f;

    // compression
    nc = CHA_IVAR[_nc];
    cm = CHA_IVAR[_cm];
    lrpk = (float) CHA_DVAR[_lrpk];
    dsm = (int *) cp[_dsm];
    dso = (int *) cp[_dso];
    Lcs = (float *) cp[_Lcs];
    Lcm = (float *) cp[_Lcm];
    Lce = (float *) cp[_Lce];
    Lmx = (float *) cp[_Lmx];
    Gcs = (float *) cp[_Gcs];
    Gcm = (float *) cp[_Gcm];
    Gce = (float *) cp[_Gce];
    Gmx = (float *) cp[_Gmx];
    Gmn = (float *) cp[_Gmn];
    Lfb = (float *) cp[_Lfb];
    Gsup = (float *) cp[_Gsup];
    Gpre = (float *) cp[_Gpre];
    gsup = (float *) cp[_gsup];
    ginc = (float *) cp[_ginc];
    zdr = (float *) cp[_zdr];
    zdi = zdr + dsm[0] * nc;
    for (i = 0; i < cs; i++) {
        dsmx = dsm[0];
        dsmo = dsmx - 1;
        dsi = (dsmx > 1);
        /* loop over filterbank channel */
        #pragma omp parallel for
        for (k = 0; k < nc; k++) {
            /* apply down-sample delay */
            for (id = dsmo; id > 0; id--) {
                kd = k * dsmx + id;
                zdr[kd] = zdr[kd - 1];
                zdi[kd] = zdi[kd - 1];
            }
            kd = k * dsmx;
            /* retrieve input */
            kr = (i + k * cs) * 2;
            ki = kr + 1;
            zdr[kd] = x[kr];
            zdi[kd] = x[ki];
            im = dsi ? ((i - dso[k]) % dsm[k]) : 0;
            if (im == 0) {              /* calculate channel levels ?? */
                agtf = (float) _hypot(zdr[kd], zdi[kd]);
                if (agtf < eps) {
                    agtf = eps;
                }
                Lfb[k] = db2(agtf / lrpk);
            }
        }
        /* loop over suppressor channel */
        #pragma omp parallel for
        for (k = 0; k < nc; k++) {
            im = dsi ? ((i - dso[k]) % dsm[k]) : 0;
            if (im == 0) {              /* calculate channel gains ?? */
                if (cm) {               /* compress mode */
                    Lsup = Lfb[k];
                }
                /* calculate compressive gain */
                if (Lsup <= Lcs[k]) {
                    Gsup[k] = Gcs[k];
                } else if (Lsup <= Lcm[k]) {
                    rnge = (Lcm[k] - Lsup) / (Lcm[k] - Lcs[k]);
                    Gsup[k] = Gcm[k] + (Gcs[k] - Gcm[k]) * rnge;
                } else if (Lsup <= Lce[k]) {
                    rnge = (Lce[k] - Lsup) / (Lce[k] - Lcm[k]);
                    Gsup[k] = Gce[k] + (Gcm[k] - Gce[k]) * rnge;
                } else {
                    Gsup[k] = Gce[k];
                }
                head = Lmx[k] - Lfb[k];
                if (Gsup[k] > head) {
                    Gsup[k] = head;
                } else if (Gsup[k] > Gmx[k]) {
                    Gsup[k] = Gmx[k];
                } else if (Gsup[k] < Gmn[k]) {
                    Gsup[k] = Gmn[k];
                }
            }
            /* apply compressive gain to filterbank outputs */
            jm = dsi ? (im + 1) % dsm[k] : 0;
            if (jm == 0) {                /* next time-step calculates gain */ 
                gsup[k] = undb2(Gsup[k]);
                Gpre[k] = Gsup[k];
            } else {                      /* interpolate gain */
                if (im == 0) {
                    ginc[k] = undb2((Gsup[k] - Gpre[k]) / dsm[k]);
                }
                gsup[k] *= ginc[k];
            }
            /* store outputs */
            kd = k * dsmx + dsmo;
            kr = (i + k * cs) * 2;
            ki = kr + 1;
            y[kr] = gsup[k] * zdr[kd];
            y[ki] = gsup[k] * zdi[kd];
        }
    }
}
