// tst_cpx - test writing & reading complex data in MAT files

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "sigpro.h"

void
print_complex(float complex *a, char *sa, float complex *b, char *sb, int n)
{
    int i;

    printf("%s = ", sa);
    for (i = 0; i < n; i++) {
	printf("%s %.1f + %.1fi", i ? "," : "", creal(a[i]), cimag(a[i]));
    }
    printf("    ");
    printf("%s = ", sb);
    for (i = 0; i < n; i++) {
	printf("%s %.1f + %.1fi", i ? "," : "", creal(b[i]), cimag(b[i]));
    }
    printf("\n");
}

void
verify_index(int i, char *vn, char *fn)
{
    if (i < 0) {
	printf("variable %s not found in %s\n", vn, fn);
        exit (1);
    }
}

void
verify_load(VAR *v, char *fn)
{
    if (v == NULL) {
	printf("Unable to load file: %s\n", fn);
        exit (1);
    }
}

int
main(int ac, char **av)
{
    float complex a[2], *c;
    float complex b[2], *d;
    int i;
    VAR *v, *w;
    static char *an = "aaaa";           // variable name 1
    static char *bn = "bbbb";           // variable name 2
    static char *fn = "cdata.mat";      // file name

    // initialize complex arrays
    a[0] = 1 + 2 * I;
    a[1] = 3 + 4 * I;
    b[0] = 5 + 6 * I;
    b[1] = 7 + 8 * I;
    print_complex(a, "a", b, "b", 2);

    // save complex variables in MATLAB data file
    v = sp_var_alloc(2);
    sp_var_set(v + 0, an, &a, 2, 1, "f8c");
    sp_var_set(v + 1, bn, &b, 2, 1, "f4c");
    sp_mat_save(fn, v);

    // load complex variables from MATLAB data file
    w = sp_mat_load(fn);
    verify_load(w, fn);
    i = sp_var_find(w, an);
    verify_index(i, an, fn);
    c = w[i].data;
    i = sp_var_find(w, bn);
    verify_index(i, bn, fn);
    d = w[i].data;
    print_complex(c, "c", d, "d", 2);

    sp_var_clear_all();
    return (0);
}

