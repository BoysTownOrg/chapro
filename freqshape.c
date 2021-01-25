// freqshape.c

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

FUNC(int) sp_freqshape(
    float *x0, 
    float *x1, float *y1, int n1,
    float *x2, float *y2, int n2
) {
    float   a, b, s;
    int     i, j, k, m;

    if (n2 < 2)
	return (1);		// table too small
    for (j = 1; j < n2; j++) {
	if (x2[j - 1] > x2[j])
	    return (2);		// table nonmontonic
    }
    if (x1 != y1)
	sp_copy(x1, y1, n1);
    sp_rcfft(y1, n1);
    j = 0;
    m = n1 / 2 + 1;
    for (i = 0; i < m; i++) {
	if ((x0[i] < x2[0]) || (x0[i] > x2[n2 - 1])) {
	    return (3);		// outside table range
	}
	while ((j < (n2 - 1)) && (x0[i] >= x2[j + 1])) {
	    j++;
	}
	if (x0[i] == x2[j]) {
	    a = y2[j];
	} else if (x0[i] == x2[j + 1]) {
	    a = y2[j + 1];
	} else {
	    b = (x2[j + 1] - x0[i]) / (x2[j + 1] - x2[j]);
	    a = b * y2[j] + (1 - b) * y2[j + 1];
	}
	s = (float) pow(10, -a / 20);
	k = i * 2;
	y1[k + 0] *= s;
	y1[k + 1] *= s;
    }
    sp_crfft(y1, n1);

    return (0);
}

FUNC(int) sp_firdb(
    float *b, 
    int n, 
    double fs, 
    float *fr, 
    float *at, 
    int nt
) {
    double a, s, f, df;
    float  *y;
    int     i, j, k, m, n1, nf, err = 0;
    static int wt = SP_WT_HAMMING;

    nf = sp_nxtpow2(n) * 2;
    df = fs / nf;
    n1 = nf / 2 + 1;
    y = (float *) calloc(nf + 2, sizeof(float));
    j = 0;
    for (i = 0; i < n1; i++) {
	f = i * df;
	if ((f < fr[0]) || (f > fr[nt - 1])) {
	    return (3);		// outside table range
	}
	while ((j < (nt - 1)) && (f >= fr[j + 1])) {
	    j++;
	}
	if (f == fr[j]) {
	    a = at[j];
	} else if (f == fr[j + 1]) {
	    a = at[j + 1];
	} else {
	    s = (fr[j + 1] - f) / (fr[j + 1] - fr[j]);
	    a = s * at[j] + (1 - s) * at[j + 1];
	}
	y[i * 2] = (float) pow(10, -a / 20);
    }
    err = sp_crfft(y, nf);
    if (err)
	return (err);
    m = n / 2;
    for (i = 0; i < m; i++) {
	j = i + m;
	k = i - m + nf;
	b[i] = y[k];
	b[j] = y[i];
    }
    err = sp_window(y, n, wt);
    for (i = 0; i < n; i++) {
	b[i] *= y[i];
    }
    free(y);

    return (err);
}

FUNC(int) sp_frqshp(
    float *x, float *y, int n, 
    int nb, double fs, 
    float * fr, float *at, int nt,
    int wrap
) {
    float  *f, *b, *z;
    int     i, m, n1, err = 0;

    if (wrap && (n == nb) && (n == sp_nxtpow2(n))) {
        n1 = n / 2 + 1;
	f = (float *) calloc(n1, sizeof(float));
	sp_linspace(f, n1, 0, fs / 2);
        err = sp_freqshape(f, x, y, n, fr, at, nt); 
	free(f);
    } else {
	b = (float *) calloc(nb, sizeof(float));
	z = (float *) calloc(nb, sizeof(float));
	err = sp_firdb(b, nb, fs, fr, at, nt);
	err = sp_fftfiltz(b, nb, x, y, n, z);
	if (!err) {
	    if (wrap) {
		for (i = 0; i < nb; i++) {
		    y[i] += z[i];
		    z[i] = y[i];
		}
	    }
	    m = nb / 2;
	    for (i = 0; i < (n - m); i++) {
		y[i] = y[i + m];
	    }
	    for (i = 0; i < m; i++) {
		y[i + n - m] = z[i];
	    }
	}
	free(b);
	free(z);
    }

    return (err);
}
