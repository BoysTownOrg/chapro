// iirfb_process.c - IIR-filterbank processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "chapro.h"
#include "cha_if.h"

#define fmove(x,y,n)    memmove(x,y,(n)*sizeof(float))

// IIR-filterbank analysis
static __inline void
analyze(CHA_PTR cp, float *x, float *y, int cs)
{
	int     dk, i, j, k, m, nc=8, ns, op=5;
	float   xx, yy, *bb, *aa, *zz, *gg, *dd, *bk, *ak, *zk, *yd;

    bb = (float *) cp[_bb];
    aa = (float *) cp[_aa];
	gg = (float *) cp[_gg];
	dd = (float *) cp[_dd];
    zz = (float *) cp[_zz];
    yd = (float *) cp[_yd];
    nc = CHA_IVAR[_nc];
    ns = CHA_IVAR[_ns];
    for (i = 0; i < cs; i++) {
        /* initialize i/o */
        xx = x[i];
        /* loop over filterbank channel */
        #pragma omp parallel for
        for (k = 0; k < nc; k++) {
			bk = bb + k * op;
			ak = aa + k * op;
			zk = zz + k * op;
            yy = (bk[0] * xx) + zk[0];
			for (j = 1; j < op; j++) {
                zk[j - 1] = (bk[j] * xx) - (ak[j] * yy) + zk[j];
			}
            /* peak-shift delay */
            m = k * ns;
			dk = (int) dd[k];
            fmove(yd + m + 1, yd + m, dk);
            yd[m] = yy;
            yy = yd[m + dk];
            /* channel out */
            y[i + k * cs] = yy * gg[k];
		}
	}
}

/***********************************************************/

// IIR-filterbank analysis
FUNC(void)
cha_iirfb_analyze(CHA_PTR cp, float *x, float *y, int cs)
{
	int     dk, i, j, k, m, nc, ns, op;
	float   xx, yy, *bb, *aa, *zz, *gg, *dd, *bk, *ak, *zk, *yd;

    bb = (float *) cp[_bb];
    aa = (float *) cp[_aa];
	gg = (float *) cp[_gg];
	dd = (float *) cp[_dd];
    zz = (float *) cp[_zz];
    yd = (float *) cp[_yd];
    nc = CHA_IVAR[_nc];
    op = CHA_IVAR[_op];
    ns = CHA_IVAR[_ns];
    for (i = 0; i < cs; i++) {
        /* initialize i/o */
        xx = x[i];
        /* loop over filterbank channel */
        #pragma omp parallel for
        for (k = 0; k < nc; k++) {
			bk = bb + k * op;
			ak = aa + k * op;
			zk = zz + k * op;
            yy = (bk[0] * xx) + zk[0];
			for (j = 1; j < op; j++) {
                zk[j - 1] = (bk[j] * xx) - (ak[j] * yy) + zk[j];
			}
            /* peak-shift delay */
            m = k * ns;
			dk = (int) dd[k];
            fmove(yd + m + 1, yd + m, dk);
            yd[m] = yy;
            yy = yd[m + dk];
            /* channel out */
            y[i + k * cs] = yy * gg[k];
		}
	}
}

// IIR-filterbank synthesis
FUNC(void)
cha_iirfb_synthesize(CHA_PTR cp, float *x, float *y, int cs)
{
    float   xsum;
    int     i, k, nc;

    nc = CHA_IVAR[_nc];
    for (i = 0; i < cs; i++) {
        xsum = 0;
        // loop over filterbank channel
        for (k = 0; k < nc; k++) {
            xsum += x[i + k * cs];
        }
        y[i] = xsum;
    }
}
