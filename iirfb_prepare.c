// iirfb_prepare.c - IIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_if.h"

/***********************************************************/


static void
root2poly(double *r, double *p, int n)
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
zp2ba(double *z, double *p, int nz, int nb, double *b, double *a)
{
    double *zk, *pk, *bk, *ak;
    int k;

    if ((nz > 0) && (nb > 0)) {
        for (k = 0; k < nb; k++) {
            zk = z + k * nz * 2;
            pk = p + k * nz * 2;
            bk = b + k * (nz + 1);
            ak = a + k * (nz + 1);
            root2poly(zk, bk, nz);
            root2poly(pk, ak, nz);
        }
    }
}

/***********************************************************/

// compute IIR-filterbank coefficients
static __inline void
iir_filterbank(CHA_PTR cp, double *b, double *a, double *g, double *d, int nc, int op, double fs)
{
    float *bb, *aa, *gg, *dd;
    int i, mxd;

    // copy IIR coefficients
    bb = (float *) cp[_bb];
    aa = (float *) cp[_aa];
    gg = (float *) cp[_gg];
    dd = (float *) cp[_dd];
    for (i = 0; i < (nc * op); i++) {
        bb[i] = (float) b[i];
        aa[i] = (float) a[i];
    }
    for (i = 0; i < nc; i++) {
        gg[i] = (float) g[i];
        dd[i] = (float) d[i];
    }
    // find maximum delay
    mxd = 0;
    for (i = 0; i < nc; i++) {
        if (mxd < (int) dd[i]) {
            mxd = (int) dd[i];
        }
    }
    CHA_IVAR[_ns] = mxd + 1;
}

/***********************************************************/

FUNC(int)
cha_iirfb_prepare(CHA_PTR cp, double *z, double *p, double *g, double *d, int nc, int nz, double fs, int cs)
{
    double *b, *a;
    int      ns, op;

    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = fs;
    // covert zeros & poles to IIR coefficients
    op = nz + 1;
    b = (double *) calloc(nc * op, sizeof(double));
    a = (double *) calloc(nc * op, sizeof(double));
    zp2ba(z, p, nz, nc, b, a);
    // allocate filter-coefficient buffers
    CHA_IVAR[_nc] = nc;
    CHA_IVAR[_op] = op;
    cha_allocate(cp, nc * op, sizeof(float), _bb);
    cha_allocate(cp, nc * op, sizeof(float), _aa);
    cha_allocate(cp, nc, sizeof(float), _gg);
    cha_allocate(cp, nc, sizeof(float), _dd);
    // save IIR-filterbank coefficients
    iir_filterbank(cp, b, a, g, d, nc, op, fs);
    ns = CHA_IVAR[_ns];
    cha_allocate(cp, nc * op, sizeof(float), _zz);
    cha_allocate(cp, nc * ns, sizeof(float), _yd);
    free(b);
    free(a);

    return (0);
}
