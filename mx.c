// mx.c - matrix functions
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/***********************************************************/

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

