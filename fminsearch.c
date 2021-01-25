// sp_fminsearch.c - parameter fitting routine using simplex method

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

#define MXNV	21		/* maximum number of parameters + 1 */
#define ALFA	1.0		/* reflection coefficient */
#define BETA	0.5		/* contraction coefficient */
#define GAMA	2.0		/* expansion coefficient */

static double  (*variance) (float *pv, void *pp);
static int     (*early_exit) ();
static void    (*report) (float *pv);
static void    *parm;

static int sfconvergence();
static void sfinflate();
static void sfhilo();
static void sfcentroid();
static void sfreflect();
static void sfaccept();
static void sfsavelo();
static void sfexpansion();
static void sfcontraction();
static void sfshrink();

static double lores, hires, lstres;
static double icons = 0.01;      /* initial constant */
static double ifrac = 0.01;      /* initial fraction */
static double ffrac = 0.0001;    /* final fraction */
static int npar, nval, hiidx, loidx, lastii;
static float simp[MXNV][MXNV];

/****************************************************************************
 *
 *
 * usage:
 *          void simpfit(iniv, npv, mxiter, mniter, pvar, prep, peex)
 *          float iniv[npv];
 *          int npv, mxiter, mniter;
 *
 * where:
 *          iniv   - array of initial parameter values
 *                   (final values will be returned here also)
 *          npv    - number of values in the parameter array
 *          mxiter - maximum number of iterations before convergence
 *          mniter - minimum number of iterations between reports
 *          pvar - pointer to function returning variance
 *          prep - pointer to function producing reports
 *          peex - pointer to function checking for "early exit"
 *          ppar - pointer to user-defined paramter data
 *
 * The user must also supply a function called "variance" that
 * will compute the residual variance for a given set of parameters
 * and a function called "report" to handle reports of intermediate
 * parameter values prior to convergence:
 *
 *         float pv[npv];
 *         void *pp;
 *         double variance(pv, pp)
 *         void report(pv)
 *
 ****************************************************************************
*/

static void 
simpfit(
    float  *iniv,
    int     npv, 
    int     mxiter, 
    int     mniter,
    double  (*pvar) (float *, void *),
    void    (*prep) (float *),
    int     (*peex) (),
    void    *ppar
)
{
    int     i;
    float   cent[MXNV], nxtv[MXNV], nxres;

    if (npv >= MXNV)
	return;
    variance = pvar;
    report = prep;
    early_exit = peex;
    parm = ppar;
    lores = hires = lstres = 0;
    lastii = 0;
    npar = npv;
    nval = npv + 1;
    sfinflate(iniv);
    sfhilo(0, 0);
    for (i = 0; i < mxiter; i++) {
        if (early_exit && early_exit())
            break;
	sfhilo(i, mniter);
	sfcentroid(cent);
	sfreflect(cent, nxtv);
	if (nxtv[npar] < hires) {
	    sfaccept(nxtv);
	    if (nxtv[npar] < lores) {
		nxres = nxtv[npar];
		sfexpansion(cent, nxtv);
		if (nxtv[npar] < nxres)
		    sfaccept(nxtv);
	    }
	} else {
	    sfcontraction(cent, nxtv);
	    if (nxtv[npar] < hires)
		sfaccept(nxtv);
	    else
		sfshrink();
	    if (sfconvergence()!= 0)
		break;
	}
    }
    sfhilo(i, 0);
    sfsavelo(iniv);

    return;
}

static void 
sfinflate(iniv)
float   iniv[MXNV];
{
    int     i, j;
    float   d;

    for (i = 0; i < npar; i++)
	simp[0][i] = iniv[i];
    simp[0][npar] = (float) variance(simp[0], parm);
    if (report) {
        report(simp[0]);
    }

    for (j = 1; j < nval; j++) {
	for (i = 0; i < npar; i++)
	    simp[j][i] = simp[0][i];
	d = simp[0][j - 1];
	simp[j][j - 1] = (float) ((d == 0) ? icons : d * (1 + ifrac));
	simp[j][npar] = (float) variance(simp[j], parm);
    }

    return;
}

static void 
sfhilo(ii, ir)
int     ii, ir;
{
    int     i;

    hiidx = 0;
    hires = simp[hiidx][npar];
    loidx = 0;
    lores = simp[loidx][npar];
    for (i = 0; i < nval; i++) {
	if (hires < simp[i][npar]) {
	    hires = simp[i][npar];
	    hiidx = i;
	}
	if (lores > simp[i][npar]) {
	    lores = simp[i][npar];
	    loidx = i;
	}
    }

    if (ir > 0) {
	if (ii < lastii + ir)
	    return;
	if (lores >= lstres)
	    return;
    }
    lastii = ii;
    lstres = lores;
    if (report) {
        report(simp[loidx]);
    }

    return;
}

static void 
sfcentroid(cent)
float   cent[MXNV];
{
    int     i, j;

    for (i = 0; i < npar; i++) {
	cent[i] = 0;
	for (j = 0; j < nval; j++)
	    if (j != hiidx)
		cent[i] = cent[i] + simp[j][i];
	cent[i] = cent[i] / npar;
    }

    return;
}

static void 
sfreflect(cent, nxtv)
float   cent[MXNV], nxtv[MXNV];
{
    int     i;

    for (i = 0; i < npar; i++)
	nxtv[i] = (float) ((1 + ALFA) * cent[i] - ALFA * simp[hiidx][i]);
    nxtv[npar] = (float) variance(nxtv, parm);

    return;
}

static void 
sfaccept(nxtv)
float   nxtv[MXNV];
{
    int     i;

    for (i = 0; i < nval; i++)
	simp[hiidx][i] = nxtv[i];

    return;
}

static void 
sfsavelo(iniv)
float   iniv[MXNV];
{
    int     i;

    for (i = 0; i < npar; i++)
	iniv[i] = simp[loidx][i];

    return;
}

static void 
sfexpansion(cent, nxtv)
float   cent[MXNV], nxtv[MXNV];
{
    int     i;

    for (i = 0; i < npar; i++)
	nxtv[i] = (float) ((1 - GAMA) * cent[i] + GAMA * simp[hiidx][i]);
    nxtv[npar] = (float) variance(nxtv, parm);

    return;
}

static void 
sfcontraction(cent, nxtv)
float   cent[MXNV], nxtv[MXNV];
{
    int     i;

    for (i = 0; i < npar; i++)
	nxtv[i] = (float) ((1 - BETA) * cent[i] + BETA * simp[hiidx][i]);
    nxtv[npar] = (float) variance(nxtv, parm);

    return;
}

static void 
sfshrink()
{
    int     i, j;

    for (j = 0; j < nval; j++) {
	for (i = 0; i < npar; i++)
	    simp[j][i] = (float) ((1 - BETA) * simp[loidx][i] + BETA * simp[j][i]);
	simp[j][npar] = (float) variance(simp[j], parm);
    }

    return;
}

static int 
sfconvergence()
{
    int     i, j;
    double  error, hi[MXNV], lo[MXNV];

    for (i = 0; i < nval; i++) {
	hi[i] = fabs(simp[0][i]);
	lo[i] = fabs(simp[0][i]);
    }
    for (j = 1; j < nval; j++) {
	for (i = 0; i < nval; i++) {
	    if (hi[i] < fabs(simp[j][i]))
		hi[i] = fabs(simp[j][i]);
	    if (lo[i] > fabs(simp[j][i]))
		lo[i] = fabs(simp[j][i]);
	}
    }
    for (i = 0; i < npar; i++) {
	error = (hi[i] > lo[i]) ? (1 - lo[i] / hi[i]) : 0;
	if (error > ffrac)
	    return (0);
    }

    return (1);
}

/*************************************************************************/

// sp_fminsearch - search for parameters with minimum error
FUNC(int)
sp_fminsearch(
    float  *iniv,
    int     npv, 
    double  (*pvar) (float *, void *),
    OPT *o,
    void *ppar
)
{
    int     mxiter, mniter;
    void    (*prep) (float *) = NULL;
    int     (*peex) () = NULL;

    if (npv >= MXNV) { // number of parameters must be
        return (1);    // less than simplex array size
    }
    icons = 0.00025;
    ifrac = 0.05;
    mxiter = 1000;
    mniter = npv;
    prep = NULL;
    peex = NULL;
    if (o) {
        if (o->icons > 0) {
            icons = o->icons;
        }
        if (o->ifrac > 0) {
            ifrac = o->ifrac;
        }
        if (o->ffrac > 0) {
            ffrac = o->ffrac;
        }
        if (o->maxiter > 0) {
            mxiter = o->maxiter;
        }
        if (o->miniter > 0) {
            mniter = o->miniter;
        }
        if (o->report) {
            prep = o->report;
        }
        if (o->escape) {
            peex = o->escape;
        }
    }
    simpfit(iniv, npv, mxiter, mniter, pvar, prep, peex, ppar);

    return (0);
}
