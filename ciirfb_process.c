// ciirfb_process.c - complex filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

#define fmove(x,y,n)    memmove(x,y,(n)*sizeof(float))
#define GTFO            4

/***********************************************************/

static __inline void
filter_tf(float *x, float *y, int cs, float *coef, float *hist, int op)
{
    float xx, yyr, yyi, *zr, *zi, *br, *bi, *ar, *ai;
    int i, ir, ii, j, no;

    no = op - 1;
    br = coef;
    bi = br + op;
    ar = bi + op;
    ai = ar + op;
    zr = hist;
    zi = zr + op;
    /* loop over time */
    for (i = 0; i < cs; i++) {
        xx = x[i];
        yyr = (br[0] * xx) + zr[0];
        yyi = (bi[0] * xx) + zi[0];
         for (j = 1; j < no; j++) {
             zr[j - 1] = (br[j] * xx) - (ar[j] * yyr - ai[j] * yyi) + zr[j];
             zi[j - 1] = (bi[j] * xx) - (ar[j] * yyi + ai[j] * yyr) + zi[j];
         }
        zr[no - 1] = (br[no] * xx) - (ar[no] * yyr - ai[no] * yyi);
        zi[no - 1] = (bi[no] * xx) - (ar[no] * yyi + ai[no] * yyr);
        /* channel out */
        ir = i * 2;
        ii = ir + 1;
        y[ir] = yyr;
        y[ii] = yyi;
    }
}

/***********************************************************/

FUNC(void)
cha_ciirfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
{
    float   *br, *zr, *bkr, *zkr;
    float   *ydr, *ydi, *yk;
    int      i, ir, ii, d, m, k, n, nc, ns, op, no = GTFO, *dn;

    dn = (int *) cp[_dn];
    nc = CHA_IVAR[_nc];
    ns = CHA_IVAR[_ns];
    ydr = (float *) cp[_ydr];
    ydi = ydr + ns * nc;
    op = no + 1;
    br = (float *) cp[_br];
    zr = (float *) cp[_zr];
    /* loop over filterbank channel */
    #pragma omp parallel for
    for (k = 0; k < nc; k++) {
        yk = y + k * cs * 2;
        bkr = br + k * op * 4;
        zkr = zr + k * op * 4;
        filter_tf(x, yk, cs, bkr, zkr, op);
        for (i = 0; i < cs; i++) {
            ir = i * 2;
            ii = ir + 1;
            /* peak shift */
            if (dn) {
                d = dn[k];
                m = k * ns;
                n = m + d;
                fmove(ydr + m + 1, ydr + m, d);
                fmove(ydi + m + 1, ydi + m, d);
                ydr[m] = yk[ir];
                ydi[m] = yk[ii];
                yk[ir] = ydr[n];
                yk[ii] = ydi[n];
            }
        }
    }
}

FUNC(void)
cha_ciirfb_synthesize(CHA_PTR cp, float *x, float *y, int cs)
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
