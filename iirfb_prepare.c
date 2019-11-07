// iirfb_prepare.c - IIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

/***********************************************************/

static void
root2poly(float *r, float *p, int n)
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
        p[i] = (float) pp[i * 2];
    }
    free(pp);
    free(qq);
}

static void
zp2sos(float *z, float *p, float g, int nsos, float *c)
{
    float b[3], a[3], *cc = c;
    int j;

    for (j = 0; j < nsos; j++) {
        root2poly(z + j * 4, b, 2);
        root2poly(p + j * 4, a, 2);
        *cc++ = b[0];
        *cc++ = b[1];
        *cc++ = b[2];
        *cc++ = -a[1];
        *cc++ = -a[2];
    }
    for (j = 0; j < 3; j++) {
        c[j] *= g;
    }
}

/***********************************************************/

// compute IIR-filterbank coefficients
static __inline int
iir_filterbank(CHA_PTR cp, float *z, float *p, float *g, int *d, int nc, int nz, double sr)
{
    float *bb, *cc, *zz, *pp;
    int i, op, mxd, *dd;

    // copy IIR coefficients
    bb = (float *) cp[_bb];
    dd = (int *) cp[_dd];
    op = nz + 1;
    mxd = 0;
    for (i = 0; i < nc; i++) {
        zz = z + i * nz * 2;
        pp = p + i * nz * 2;
        cc = bb + i * op * 2;
        zp2sos(zz, pp, g[i], nz / 2, cc);
        dd[i] = d[i];
        if (mxd < d[i]) {
            mxd = d[i];
        }
    }
    return (mxd + 1); // size of peak-delay buffer
}

/***********************************************************/

FUNC(int)
cha_iirfb_prepare(CHA_PTR cp, float *z, float *p, float *g, int *d, int nc, int nz, double sr, int cs)
{
    int nn, op, ncoef, nhist;

    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_DVAR[_fs] = sr / 1000;
    CHA_IVAR[_cs] = cs;
    // allocate filter-coefficient buffers
    op = nz + 1;
    CHA_IVAR[_nc] = nc;
    CHA_IVAR[_op] = op;
    nhist = 2 * nz;
    ncoef = 5 * (nz / 2);
    cha_allocate(cp, nc * ncoef, sizeof(float), _bb);
    cha_allocate(cp, nc * nhist, sizeof(float), _zz);
    cha_allocate(cp, nc, sizeof(int), _dd);
    // save IIR-filterbank coefficients
    nn = iir_filterbank(cp, z, p, g, d, nc, nz, sr);
    CHA_IVAR[_nn] = nn;
    cha_allocate(cp, nc * nn, sizeof(float), _yd);
    // allocate input & chunk buffers
    cha_allocate(cp, nc * cs, sizeof(float), _xx);
    cha_allocate(cp, nc * cs, sizeof(float), _cc);

    return (0);
}
