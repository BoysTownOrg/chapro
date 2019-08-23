// ifft.c - inverse FFT

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

/***********************************************************/

// complex-to-complex FFT

FUNC(void)
cha_fft(float *x, int n)
{
    int m;

    m = ilog2(n);
    cfft2(x, m, 1);
}

// complex-to-complex inverse FFT

FUNC(void)
cha_ifft(float *x, int n)
{
    int i, m;

    m = ilog2(n);
    cfft2(x, m, 0);
    // scale inverse by 1/n
    for (i = 0; i < n * 2; i++) {
        x[i] /= n;
    }
}
