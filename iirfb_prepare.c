// iirfb_prepare.c - IIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_if.h"

/***********************************************************/

// compute IIR-filterbank coefficients
static __inline void
iir_filterbank(CHA_PTR cp, double *b, double *a, double *g, double *d, int nc, int op, double fs)
{
	float *bb, *aa, *gg, *dd;
	int i, mxd;

	bb = (float *) cp[_bb];
	aa = (float *) cp[_aa];
	gg = (float *) cp[_gg];
	dd = (float *) cp[_dd];
	// copy IIR coefficients
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
cha_iirfb_prepare(CHA_PTR cp, double *b, double *a, double *g, double *d, int nc, int op, double fs, int cs)
{
    int      ns;

    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    CHA_DVAR[_fs] = fs;
    // allocate filter-coefficient buffers
    CHA_IVAR[_nc] = nc;
    CHA_IVAR[_op] = op;
    cha_allocate(cp, nc * op, sizeof(float), _bb);
    cha_allocate(cp, nc * op, sizeof(float), _aa);
    cha_allocate(cp, nc, sizeof(float), _gg);
    cha_allocate(cp, nc, sizeof(float), _dd);
    // compute IIR-filterbank coefficients
    iir_filterbank(cp, b, a, g, d, nc, op, fs);
    ns = CHA_IVAR[_ns];
    cha_allocate(cp, nc * op, sizeof(float), _zz);
    cha_allocate(cp, nc * ns, sizeof(float), _yd);

    return (0);
}
