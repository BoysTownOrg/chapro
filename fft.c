// fft.c - complex FFT and IFFT for chapro firmware

#include <stdlib.h>
#include <math.h>
#include "chapro.h"

//-----------------------------------------------------------

static int
ilog2(int n)
{
    int     m;

    for (m = 1; m < 32; m++)
	if (n == (1 << m))
	    return (m);
    return (-1);
}

static int
cfft2(float *x, int m, int d)
{
    double  c1, c2, z;
    float  *y, tx, ty, t1, t2, u1, u2;
    int     i, i1, j, k, i2, l, l1, l2, n;
    int     ii, jj;

    n = 1 << m;

// re-order

    i2 = n >> 1;
    j = 0;
    y = x + 1;
    for (i = 0; i < n - 1; i++) {
	if (i < j) {
            ii = i * 2;
            jj = j * 2;
	    tx = x[ii];
	    ty = y[ii];
	    x[ii] = x[jj];
	    y[ii] = y[jj];
	    x[jj] = tx;
	    y[jj] = ty;
	}
	k = i2;
	while (k <= j) {
	    j -= k;
	    k >>= 1;
	}
	j += k;
    }

// radix-2 Fourier transform

    c1 = -1;
    c2 = 0;
    l2 = 1;
    for (l = 0; l < m; l++) {
	l1 = l2;
	l2 <<= 1;
	u1 = 1;
	u2 = 0;
	for (j = 0; j < l1; j++) {
	    for (i = j; i < n; i += l2) {
		i1 = i + l1;
		ii = i * 2;
		jj = i1 * 2;
		t1 = u1 * x[jj] - u2 * y[jj];
		t2 = u1 * y[jj] + u2 * x[jj];
		x[jj] = x[ii] - t1;
		y[jj] = y[ii] - t2;
		x[ii] += t1;
		y[ii] += t2;
	    }
	    z = u1 * c1 - u2 * c2;
	    u2 = (float) (u1 * c2 + u2 * c1);
	    u1 = (float) z;
	}
	c2 = sqrt((1.0 - c1) / 2.0);
	c1 = sqrt((1.0 + c1) / 2.0);
	if (d)
	    c2 = -c2;
    }

    return (0);
}

//-----------------------------------------------------------

static int
cdft(float *x, int n, int d)
{
    double a, ur, ui, vr, vi, wr, wi, xr, xi, yr, yi, zz;
    float *x0, *x1, *y0, *y1, *y;
    int i, ii, j, jj, m;
    static double tpi = 2 * M_PI;

    m = n * 2;
    y = (float *) calloc(m, sizeof(float));
    x0 = x;
    y0 = y;
    y1 = y + 1;
    x1 = x + 1;
    a = -tpi / n;
    if (d)
	a = -a;
    ur = cos(a);
    ui = sin(a);
    for (i = 0; i < n; i++) {
	if (i == 0) {
	    vr = 1;
	    vi = 0;
	} else {
	    zz = ur * vr - ui * vi;
	    vi = ur * vi + ui * vr;
	    vr = zz;
	}
        xr = x0[0];
        xi = x1[0];
	yr = xr * vr - xi * vi;
	yi = xr * vi + xi * vr;
	wr = vr;
	wi = vi;
	for (j = 1; j < n; j++) {
	    zz = vr * wr - vi * wi;
	    wi = vr * wi + vi * wr;
	    wr = zz;
	    jj = j * 2;
	    xr = x0[jj];
	    xi = x1[jj];
	    yr += xr * wr - xi * wi;
	    yi += xr * wi + xi * wr;
	}
	ii = i * 2;
	y0[ii] = (float) yr;
	y1[ii] = (float) yi;
    }
    for (i = 0; i < n; i++) {
	x0[i] = y0[i];
	x1[i] = y1[i];
    }
    free(y);

    return (0);
}

//-----------------------------------------------------------

FUNC(int) cha_fft(
    float *x, int n
) {
    int m, err;

    if (n < 2)
	return (0);

    m = ilog2(n);
    if (m > 0)
        err = cfft2(x, m, 0);
    else
        err = cdft(x, n, 0);

    return (err);
}

FUNC(int) cha_ifft(
    float *x, int n
) {
    int m, err, i;

    if (n < 2)
	return (0);

    m = ilog2(n);
    if (m > 0)
        err = cfft2(x, m, 1);
    else
        err = cdft(x, n, 1);

// scale inverse by 1/n

    for (i = 0; i < n * 2; i++) {
        x[i] /= n;
    }

    return (err);
}
