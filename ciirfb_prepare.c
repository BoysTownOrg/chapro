// ciirfb_prepare.c - complex filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_gf.h"

/***********************************************************/

static void
root2poly(
    double *yr, double *yi,
    float  *x,
    float  *g,
    int n
) {
    double ytemp;
    int i, j, jr, ji;

    dzero(yr, n + 1);
    dzero(yi, n + 1);
    yr[0] = 1;
    for (j = 0; j < n; j++) {
        jr = j * 2;
        ji = jr + 1;
        for (i = j; i >= 0; i--) {
            yr[i + 1] -= yr[i] * x[jr] - yi[i] * x[ji];
            yi[i + 1] -= yi[i] * x[jr] + yr[i] * x[ji];
        }
    }
    if (g) {
        for (i = 0; i <= n; i++) {
            ytemp = yr[i] * g[0] - yi[i] * g[1];
            yi[i] = yi[i] * g[0] + yr[i] * g[1];
            yr[i] = ytemp;
        }
    }
}

static void
zp2tf(
    double *coef,
    float *z, float *p,
    int no, int nf
) {
    double *br, *bi, *ar, *ai;
    float  *zz, *pp;
    int     j, op;

    op = no + 1;
    for (j = 0; j < nf; j++) {
        zz = z + j * no * 2;
        pp = p + j * no * 2;
        br = coef + j * op * 4;
        bi = br + op;
        ar = bi + op;
        ai = ar + op;
        root2poly(br, bi, zz, 0, no);
        root2poly(ar, ai, pp, 0, no);
    }
}

/***********************************************************/

static void
filterbank_prepare(CHA_PTR cp, float *z, float *p, float *g, int *d, 
                   int nc, int no, double sr, int cs)
{
    double  fs, yyr, yyi, gnr, gni, *bkr, *bki, *br, *coef;
    int     i, j, k, mxd, ns, op, *dn;

    dn = (int *) cha_allocate(cp, nc, sizeof(int), _dn);
    for (i = 0; i < nc; i++) {
        dn[i] = d[i];
    }
    //-----------------------------
    op = no + 1;
    coef = (double *) calloc(nc * op * 4, sizeof(double));
    br = (double *) cha_allocate(cp, nc * op * 4, sizeof(double), _br);
    zp2tf(coef, z, p, no, nc);
    for (i = 0; i < nc * op * 4; i++) {
        br[i] = coef[i];
    }
    free(coef);
    //-----------------------------
    dn = (int *) cp[_dn];
    mxd = 0;
    if (dn) {      /* find maximum delay */
        for (k = 0; k < nc; k++) {
            if (mxd < dn[k]) {
                mxd = dn[k];
            }
        }
    }
    ns = mxd + 1;
    CHA_IVAR[_ns] = ns;
    if (g) {     /* apply complex gain to filter coefficients */
        for (k = 0; k < nc; k++) {
            gnr = g[k * 2];
            gni = g[k * 2 + 1];
            bkr = br + k * op * 4;
            bki = bkr + op;
            for (j = 0; j < op; j++) {
                yyr = gnr * bkr[j] - gni * bki[j];
                yyi = gnr * bki[j] + gni * bkr[j];
                bkr[j] = yyr;
                bki[j] = yyi;
            }
        }
    }
    cha_allocate(cp, nc * op * 4, sizeof(double), _zr);
    cha_allocate(cp, ns * nc * 2, sizeof(float), _ydr);
    //-----------------------------
    fs = sr / 1000;
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = fs;
    CHA_IVAR[_nc] = nc;
}

/***********************************************************/

FUNC(int)
cha_ciirfb_prepare(CHA_PTR cp, float *z, float *p, float *g, int *d, 
                   int nc, int no, double sr, int cs)
{
    cha_prepare(cp);
    filterbank_prepare(cp, z, p, g, d, nc, no, sr, cs);

    return (0);
}
