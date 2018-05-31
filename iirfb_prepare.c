// iirfb_prepare.c - IIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_if.h"

/***********************************************************/


static void
root2poly(float *r, double *p, int n)
{
    double *pp, *qq;
    int i, ir, ii, j, jr, ji;

    pp = (double *) calloc((n + 1) * 2, sizeof(double));
    qq = (double *) calloc((n + 1) * 2, sizeof(double));
    dzero(pp, (n + 1) * 2);
    dzero(qq, (n + 1) * 2);
    pp[0] = qq[0] = 1;
    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = i * 2 + 1;
        qq[2] = pp[2] - r[ir];
        qq[3] = pp[3] - r[ii];
        for (j = 0; j < i; j++) {
            jr = j * 2;
            ji = j * 2 + 1;
            qq[jr + 4] = pp[jr + 4] - (pp[jr + 2] * r[ir] - pp[ji + 2] * r[ii]);
            qq[ji + 4] = pp[ji + 4] - (pp[ji + 2] * r[ir] + pp[jr + 2] * r[ii]);
        }
        dcopy(pp, qq, (n + 1) * 2);
    }
    // return real part of product-polynomial coefficients
    for (i = 0; i < (n + 1); i++) {
        p[i] = pp[i * 2];
    }
    free(pp);
    free(qq);
}

static void
zp2tf(float *z, float *p, int nz, int nb, double *b, double *a)
{
    double *bk, *ak;
    float  *zk, *pk;
    int k;

    for (k = 0; k < nb; k++) {
        zk = z + k * nz * 2;
        pk = p + k * nz * 2;
        bk = b + k * (nz + 1);
        ak = a + k * (nz + 1);
        root2poly(zk, bk, nz);
        root2poly(pk, ak, nz);
    }
}

/***********************************************************/

// compute IIR-filterbank coefficients
static __inline int
iir_filterbank(CHA_PTR cp, float *z, float *p, float *g, int *d, int nc, int op, double fs)
{
    double *b, *a;
    float *bb, *aa;
    int i, j, mxd, *dd;

    // convert zeros & poles to IIR coefficients
    b = (double *) calloc(nc * op, sizeof(double));
    a = (double *) calloc(nc * op, sizeof(double));
    zp2tf(z, p, op - 1, nc, b, a);
    // copy IIR coefficients
    bb = (float *) cp[_bb];
    dd = (int *) cp[_dd];
    aa = bb + op;
    mxd = 0;
    for (i = 0; i < nc; i++) {
        // copy coefficients
        for (j = 0; j < op; j++) {
            bb[i * op * 2 + j] = (float) b[i * op + j] * g[i];
            aa[i * op * 2 + j] = (float) a[i * op + j];
        }
        // copy delay & save maximum
        dd[i] = d[i];
        if (mxd < d[i]) {
            mxd = d[i];
        }
    }
    free(b);
    free(a);
    return (mxd + 1); // size of peak-delay buffer
}

/***********************************************************/

FUNC(int)
cha_iirfb_prepare(CHA_PTR cp, float *z, float *p, float *g, int *d, int nc, int nz, double fs, int cs)
{
    int      ns, op;

    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_DVAR[_fs] = fs;
    CHA_IVAR[_cs] = cs;
    // allocate filter-coefficient buffers
    op = nz + 1;
    CHA_IVAR[_nc] = nc;
    CHA_IVAR[_op] = op;
    cha_allocate(cp, nc * op * 2, sizeof(float), _bb);
    cha_allocate(cp, nc, sizeof(int), _dd);
    // save IIR-filterbank coefficients
    ns = iir_filterbank(cp, z, p, g, d, nc, op, fs);
    CHA_IVAR[_ns] = ns;
    cha_allocate(cp, nc * op * 2, sizeof(float), _zz);
    cha_allocate(cp, nc * ns, sizeof(float), _yd);

    return (0);
}
