// iirfb.c - IIR-filterbank parameter computation

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/***********************************************************/

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880 
#endif

#define fmin(x,y)       ((x<y)?(x):(y))
#define fmove(x,y,n)    memmove(x,y,(n)*sizeof(float))
#define fcopy(x,y,n)    memcpy(x,y,(n)*sizeof(float))
#define fzero(x,n)      memset(x,0,(n)*sizeof(float))
#define dcopy(x,y,n)    memcpy(x,y,(n)*sizeof(double))
#define round(x)        ((int)floorf((x)+0.5))
#define log2(x)         (logf(x)/M_LN2)

/***********************************************************/

// p2ss - transform poles to state-space (assumes np is even)
static void
p2ss(float *p, int np, float *a, float *b, float *c, float *d)
{
    int k, kd, kr, ki;

    fzero(a, np*np);
    fzero(b, np);
    fzero(c, np);
    for (k = 1; k < np; k++) {
        a[k * (np + 1) - 1] = 1;
    }
    for (k = 0; k < (np / 2); k++) {
        kr = k * 2;
        ki = kr + 1;
        kd = 2 * k * (np + 1); 
        a[kd] = p[ki] * ((p[ki] < 0) ? 2 : -2);
        a[kd + 1] = -1;
    }
    b[0] = 1;
    c[np - 1] = 1;
    d[0] = 0;
}

// mxinv - in-place matrix inversion
static void
mxinv(float *a, int n)
{
    float at, *b;
    int i, j, jj, k;
    
    // create identity matrix
    b = (float *) calloc(n * n, sizeof(float));
    // work down from top
    for (j = 0; j < n; j++) {
        k = j * (n + 1);  // index of diagonal element
        at = 1 / a[k];     // assume diagonal element is not close to zero
        for (jj = 0; jj < n; jj++) {
            a[j * n + jj] *= at;
        }
        a[k] = at;
        if (j < (n - 1)) {
            for (i = (j + 1); i < n; i++) {
                k = i * n + j;
                at = -a[k];
                a[k] = 0;
                for (jj = 0; jj < n; jj++) {
                    a[i * n + jj] += a[j * n + jj] * at;
                }
            }
        }
    }
    // copy lower-triangle of a to b
    fzero(b, n * n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < (i + 1); j++) {
            k = i * n + j;
            b[k] = a[k];
        }
    }
    // copy work up from bottom
    for (i = (n - 2); i >= 0; i--) {
        for (jj = (i + 1); jj < n; jj++) {
            for (j = 0; j < n; j++) {
                b[i * n + j] -= a[i * n + jj] * b[jj * n + j];
            }
        }
    }
    // copy b to a
    fcopy(a, b, n * n);
    free(b);
}

// mxmul - matrix multiply
static void
mxmul(float *a, float *b, float *c, int m, int o, int n, double d)
{
    float *e, s, t;
    int i, j, k;

    t = (float) d;
    e = calloc(m * n, sizeof(float));
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            s = 0;
            for (k = 0; k < o; k++) {
                s += a[i * o + k] * b[k * n + j];
            }
            e[i * n + j] = s * t;
        }
    }
    fcopy(c, e, m * n);
    free(e);
}

// mxadd - matrix add
static void
mxadd(float *a, float *b, float *c, int m, int n)
{
    int i, j, k;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
			k = i * n + j;
            c[k] = a[k] + b[k];
        }
    }
}

// ssblt - state-space bilinear transformaion
static void
ssblt(float *a, float *b, float *c, float *d, int n)
{
    float *t1, *t2, *e, f[1];
    int i, j, k;

    e = (float *) calloc(n, sizeof(float));
    fcopy(e, b, n);
    fcopy(f, d, 1);
    t1 = (float *) calloc(n * n, sizeof(float));
    t2 = (float *) calloc(n * n, sizeof(float));
    for (i = 0; i < n; i++) {
        k = i * (n + 1);
        t1[k] = 1;
        t2[k] = 1;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            k = i * n + j;
            t1[k] += a[k] / 4;
            t2[k] -= a[k] / 4;
        }
    }
    mxinv(t2, n);
    mxmul(t2, t1, a, n, n, n, 1);
    mxmul(t2,  b, b, n, n, 1, M_SQRT2 / 2);
    mxmul( c, t2, c, 1, n, n, M_SQRT2 / 2);
    mxmul( c,  e, d, 1, n, 1, M_SQRT2 / 4);
    mxadd( d,  f, d, 1, 1);
    free(e);
    free(t1);
    free(t2);
}

static float
trace(float *a, int n)
{
	float s;
	int   i;

	s = 0;
	for (i = 0; i < n; i++) {
		s += a[i * (n + 1)];
	}
	return s;
}

static void
poly(float *a, float *p, int n)
{
	int i, j;
	float *b, *e;

	b = (float *) calloc(n * n, sizeof(float));
	e = (float *) calloc(n * n, sizeof(float));
	fcopy(b, a, n * n);
	fzero(e, n * n);
	p[0] = 1;
	for (i = 1; i <= n; i++) {
		if (i > 1) {
			for (j = 0; j < n; j++) {
				e[j * (n + 1)] = p[i - 1];
			}
            mxadd(b, e, e, n, n);
            mxmul(a, e, b, n, n, n, 1);
		}
		p[i] = -trace(b, n) / i;
	}
	free(b);
	free(e);
}

// ss2ba - state-space to IIR
static void
ss2ba(float *a, float *b, float *c, float *d, int n, float *bb, float *aa)
{
	float *e, *f, *g;
	int i;

	e = calloc(n * n, sizeof(float));
	f = calloc(n * n, sizeof(float));
	g = calloc(n + 1, sizeof(float));
	// denominator
	poly(a, aa, n);
	// numerator
    mxmul(b, c, e, n, 1, n, -1);
    mxadd(a, e, f, n, n);
	poly(f, g, n);
	for (i = 0; i <= n; i++) {
		bb[i] = g[i] + (d[0] - 1) * aa[i];
	}
	free(e);
	free(f);
	free(g);
}

/***********************************************************/

// ss2lp - convert normalized state-space to low-pass
static void
ss2lp(float *wn, int n, float *a, float *b, float *c, float *d)
{
    int k;

    for (k = 0; k < (n * n); k++) {
        a[k] *= wn[0];
    }
    for (k = 0; k < n; k++) {
        b[k] *= wn[0];
    }
}

// ss2hp - convert normalized state-space to high-pass
static void
ss2hp(float *wn, int n, float *a, float *b, float *c, float *d)
{
	float *e, *f;
    int k;

	e = (float *) calloc(n * n, sizeof(float));
	f = (float *) calloc(n, sizeof(float));
	fcopy(e, a, n * n);
	fcopy(f, b, n);
	mxinv(e, n);
    for (k = 0; k < (n * n); k++) {
        a[k] = wn[0] * e[k];
    }
	mxmul(e, b, b, n, n, 1, -wn[0]);
	mxmul(c, e, c, 1, n, n, 1);
	mxmul(c, f, e, 1, n, 1, 1); 
	d[0] -= e[0];
	free(e);
}

// ss2bp - convert normalized state-space to band-pass
static void
ss2bp(float *wn, int np, float *a, float *b, float *c, float *d)
{
}

// ss2bs - convert normalized state-space to band-stop
static void
ss2bs(float *wn, int np, float *a, float *b, float *c, float *d)
{
}

// lp2ba - transform normalized all-pole low-pass to IIR
static void
lp2ba(float *p, int np, float *b, float *a, float *wn, int ft)
{
    float *ssa, *ssb, *ssc, ssd[1], u[2];

    // transform poles to state-space
    ssa = (float *) calloc(np * np, sizeof(float));
    ssb = (float *) calloc(np, sizeof(float));
    ssc = (float *) calloc(np, sizeof(float));
    p2ss(p,np,ssa,ssb,ssc,ssd);
    if ((ft == 0) || (ft == 1)) {
        u[0] = (float) (4 * tan(M_PI * wn[0] / 2));
	} else if ((ft == -1) || (ft == 2)) {
        u[0] = (float) (4 * tan(M_PI * wn[0] / 2));
        u[1] = (float) (4 * tan(M_PI * wn[1] / 2));
	} else {
		fprintf(stderr, "*** Unknown filter type = %d.\n", ft);
		exit(1);
	}
    if (ft < 0) {
        ss2bs(u,np,ssa,ssb,ssc,ssd);
	} else if (ft == 0) {
        ss2lp(u,np,ssa,ssb,ssc,ssd);
	} else if (ft == 1) {
        ss2hp(u,np,ssa,ssb,ssc,ssd);
	} else if (ft == 2) {
        u[0] = (float) (4 * tan(M_PI * wn[0] / 2));
        ss2bp(u,np,ssa,ssb,ssc,ssd);
    }
    ssblt(ssa,ssb,ssc,ssd, np);
    ss2ba(ssa,ssb,ssc,ssd, np, b, a);
    free(ssa);
    free(ssb);
    free(ssc);
}

/**********************************************************/

// butterp - normalized low-pass Butterworth filter
static void
butterp(float *p, int n)
{
    double aa;
    int i, ir, ii, m;

    if (n < 1) {
        return;
    }
    m = n / 2;
    for (i = 0; i < m; i++) {
        ir = 2 * i;
        ii = ir + 1;
        aa = ii * M_PI / (2 * n);
        p[ir] = (float) (-sin(aa));
        p[ii] = (float) (cos(aa));
    }
    if (n % 2) {
        p[n - 1] = -1;
    }
}

static void
butter(	            // Butterworth filter design
    float *b, 		// input coeffcients
    float *a, 		// output coeffcients
    int n,		    // filter order
    float *wn,		// cutoff frequency
    double fs,		// sampling frequency
    int ft		    // filter type 
)
{
    float *p;

    p = (float *) calloc(n, sizeof(float));
    butterp(p, n);
    lp2ba(p, n, b, a, wn, ft);
    free(p);
}

/***********************************************************/

// compute IIR-filterbank coefficients
void
iirfb(double *b, double *a, double *g, double *d, int *nc, int *op)
{
	double *s, cos=9;
	float *bb, *aa, *gg, *dd, *wn, *bk, *ak;
	int i, k, mxd;
    static double _b[8*5] = {
    2.39658456e-06, 9.58633825e-06, 1.43795074e-05, 9.58633825e-06, 2.39658456e-06,
    4.68723083e-04, 0.00000000e+00,-9.37446165e-04, 0.00000000e+00, 4.68723083e-04,
    1.24747171e-03, 0.00000000e+00,-2.49494342e-03, 0.00000000e+00, 1.24747171e-03,
    3.19250430e-03, 0.00000000e+00,-6.38500859e-03, 0.00000000e+00, 3.19250430e-03,
    7.88234386e-03, 0.00000000e+00,-1.57646877e-02, 0.00000000e+00, 7.88234386e-03,
    1.87984229e-02, 0.00000000e+00,-3.75968459e-02, 0.00000000e+00, 1.87984229e-02,
    4.31011542e-02, 0.00000000e+00,-8.62023084e-02, 0.00000000e+00, 4.31011542e-02,
	1.49131039e-01,-5.96524156e-01, 8.94786233e-01,-5.96524156e-01, 1.49131039e-01};
	static double _a[8*5] = {
    1.00000000e+00,-3.78901370e+00, 5.38901420e+00,-3.40967042e+00, 8.09708267e-01,
    1.00000000e+00,-3.91609799e+00, 5.77268885e+00,-3.79618186e+00, 9.39709137e-01,
    1.00000000e+00,-3.84380066e+00, 5.59364362e+00,-3.65174519e+00, 9.02627286e-01,
    1.00000000e+00,-3.70143260e+00, 5.26455707e+00,-3.40542867e+00, 8.46700113e-01,
    1.00000000e+00,-3.41182853e+00, 4.65527387e+00,-2.98260118e+00, 7.65145005e-01,
    1.00000000e+00,-2.81545747e+00, 3.57693518e+00,-2.26695910e+00, 6.51625525e-01,
    1.00000000e+00,-1.62765417e+00, 2.00985702e+00,-1.14453709e+00, 5.05359476e-01,
	1.00000000e+00,-6.16229925e-01, 6.06292132e-01,-1.38330252e-01, 2.52443143e-02};
    static double _g[8] = {1,0.51286,0.73348,0.77050,0.67926,0.78124,0.87297,-1};
    static double _d[8] = {23,1,22,36,44,50,53,58};
	static int m=8, n=5, no=2;
    static double fs = 24000;
    static double cf[7] = {317.1666,502.9734,797.6319,1264.9,2005.9,3181.1,5044.7};

	// copy precomputed coefficients
	dcopy(b,_b,m*n);
	dcopy(a,_a,m*n);
	dcopy(g,_g,m);
	dcopy(d,_d,m);
	*nc = m;
	*op = n;
	// allocate arrays
	bb = (float *) calloc(m * n, sizeof(float));
	aa = (float *) calloc(m * n, sizeof(float));
	gg = (float *) calloc(m, sizeof(float));
	dd = (float *) calloc(m, sizeof(float));
	// copy in
	for (i = 0; i < (m * n); i++) {
		bb[i] = (float) b[i];
		aa[i] = (float) a[i];
	}
	// fetch Butterworth coefficients
	s = (double *) calloc(m, sizeof(double));
	wn = (float *) calloc(m, sizeof(float));
	for (i = 0; i < m; i++) {
		s[i] = 1 + cos / cf[i];
	}
	// low-pass filter
	bk = bb;
	ak = aa;
	wn[0] = (float) ((cf[0] / s[0]) * (2 / fs));
	butter(bk, ak, 2 * no, wn, fs, 0);
	// band-pass filters
	for (k = 1; k < (m - 1); k++) {
        bk = bb + k * n;
        ak = aa + k * n;
        wn[0] = (float) ((cf[k-1] * s[k-1]) * (2 / fs));
        wn[1] = (float) ((cf[k] / s[k]) * (2 / fs));
        //butter(bk, ak, no, wn, fs, 2);
	}
	// high-pass filter
    k = m - 1;
    bk = bb + k * n;
    ak = aa + k * n;
    wn[0] = (float) ((cf[k-1] * s[k-1]) * (2 / fs));
    //butter(bk, ak, 2 * no, wn, fs, 1);
	// copy out
	for (i = 0; i < (m * n); i++) {
		b[i] = bb[i];
		a[i] = aa[i];
	}
}
