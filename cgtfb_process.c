// filterbank.c - gammatone-filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_gf.h"

#define fmove(x,y,n)    memmove(x,y,(n)*sizeof(float))
#define GTFO            4

/***********************************************************/

FUNC(void)
cha_cgtfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
{
    double  *br, *bi, *ar, *ai, *zr, *zi;
    double  *bkr, *bki, *akr, *aki, *zkr, *zki;
    double   xx, yyr, yyi;
    float   *ydr, *ydi, *zg;
    int      i, d, m, k, kr, ki, n, nc, ns, op, no = GTFO, *dn;

    dn = (int *) cp[_dn];
    nc = CHA_IVAR[_nc];
    ns = CHA_IVAR[_ns];
    ydr = (float *) cp[_ydr];
    ydi = ydr + ns * nc;
    op = no + 1;
    br = (double *) cp[_br];
    ar = (double *) cp[_ar];
    zr = (double *) cp[_zr];
    zg = (float *) cp[_zg];
    bi = br + nc * op;
    ai = ar + nc * op;
    zi = zr + nc * op;
    for (i = 0; i < cs; i++) {
        /* initialize i/o */
        xx = x[i];
        /* loop over filterbank channel */
        #pragma omp parallel for
        for (k = 0; k < nc; k++) {
            bkr = br + k * op;
            bki = bi + k * op;
            akr = ar + k * op;
            aki = ai + k * op;
            zkr = zr + k * op;
            zki = zi + k * op;
            /* assume:
               bbr[0] = bbi[0] = 0;
               bbr[4] = bbi[4] = 0;
            */
            yyr = zkr[0];
            yyi = zki[0];
            zkr[0] = (bkr[1] * xx) - (akr[1] * yyr - aki[1] * yyi) + zkr[1];
            zki[0] = (bki[1] * xx) - (akr[1] * yyi + aki[1] * yyr) + zki[1];
            zkr[1] = (bkr[2] * xx) - (akr[2] * yyr - aki[2] * yyi) + zkr[2];
            zki[1] = (bki[2] * xx) - (akr[2] * yyi + aki[2] * yyr) + zki[2];
            zkr[2] = (bkr[3] * xx) - (akr[3] * yyr - aki[3] * yyi) + zkr[3];
            zki[2] = (bki[3] * xx) - (akr[3] * yyi + aki[3] * yyr) + zki[3];
            zkr[3] =               - (akr[4] * yyr - aki[4] * yyi);
            zki[3] =               - (akr[4] * yyi + aki[4] * yyr);
            /* peak shift */
            if (dn) {
                d = dn[k];
                m = k * ns;
                n = m + d;
                fmove(ydr + m + 1, ydr + m, d);
                fmove(ydi + m + 1, ydi + m, d);
                ydr[m] = (float) yyr;
                ydi[m] = (float) yyi;
                yyr = ydr[n];
                yyi = ydi[n];
            }
            /* channel out */
            kr = (i + k * cs) * 2;
            ki = kr + 1;
            y[kr] = ((float) yyr) * zg[k];
            y[ki] = ((float) yyi) * zg[k];
        }
    }
}

FUNC(void)
cha_cgtfb_synthesize(CHA_PTR cp, float *x, float *y, int cs)
{
    float sm;
    int i, k, kr, nc;

    nc = CHA_IVAR[_nc];
    for (i = 0; i < cs; i++) {
        sm = 0;
        /* loop over filterbank channel */
        for (k = 0; k < nc; k++) {
            kr = (i + k * cs) * 2;
            sm += x[kr];
        }
        y[i] = sm;
    }
}
