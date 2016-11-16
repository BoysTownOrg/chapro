// cha_core.c - CHAPRO core functions

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "chapro.h"
#include "version.h"

#define free_null(p)    if(p){free(p);p=NULL;}

/***********************************************************/

FUNC(char *) 
cha_version(void)
{
    return (VER);                  
};

FUNC(void) 
cha_prepare(CHA_PTR cp)
{
    if (cp[_size] == NULL) {
        cp[_size] = calloc(NPTR, sizeof(int));
        cp[_ivar] = calloc(NVAR, sizeof(int));
        cp[_dvar] = calloc(NVAR, sizeof(double));
        ((int *)cp[_size])[_size] = NPTR * sizeof(int);
        ((int *)cp[_size])[_ivar] = NVAR * sizeof(int);
        ((int *)cp[_size])[_dvar] = NVAR * sizeof(double);
    }
};

FUNC(void *) 
cha_allocate(CHA_PTR cp, int cnt, int siz, int idx)
{
    cha_prepare(cp);
    assert(idx < NPTR);
    free_null(cp[idx]);
    cp[idx] = calloc(cnt, siz);
    ((int *)cp[_size])[idx] = cnt * siz;

    return (cp[idx]);
};

FUNC(void) 
cha_cleanup(CHA_PTR cp)
{
    int i;

    for (i = 0; i < NPTR; i++) {
        free_null(cp[i]);
    }
};

FUNC(int)
cha_data_gen(CHA_PTR cp, char *fn)
{
    int ptsiz, arsiz, arlen, i, j, *cpsiz;
    unsigned long *ulptr;
    FILE *fp;
    static char *head[] = {
        "#ifndef CHA_DATA_H",
        "#define CHA_DATA_H",
        ""
    };
    static char *tail[] = {
        "",
        "#endif // CHA_DATA_H"
    };
    static int hdsz = sizeof(head) / sizeof(char *);
    static int tlsz = sizeof(tail) / sizeof(char *);
    static int ptpl = 16;
    static int arpl = 8;
    static int ivpl = 8;
    static int dvpl = 5;

    fp = fopen(fn, "wt");
    if (fp == NULL) {
        return (1);
    }
    cpsiz = (int *) cp[_size];
    if (cpsiz == NULL) {
        fclose(fp);
        return (2);
    }
    // print header
    arsiz = 0;
    for (i = 0; i < NPTR; i++) {
        arsiz += cpsiz[i];
    }
    if (arsiz == 0) {
        fclose(fp);
        return (3);
    }
    fprintf(fp, "// cha_data.h - array size = %d bytes\n", arsiz);
    for (i = 0; i < hdsz; i++) {
        fprintf(fp, "%s\n", head[i]);
    }
    // initialize ptr arrays
    ptsiz = 0;
    for (i = 0; i < NPTR; i++) {
        if (cp[i]) ptsiz = i + 1;
    }
    for (i = 0; i < ptsiz; i++) {
        if (i == _size) {
            arlen = cpsiz[i] / sizeof(long);
            arsiz = 0;
            ulptr = (unsigned long *) cp[i];
            if (ulptr) {
                for (j = 0; j < arlen; j++) {
                    if (ulptr[j]) arsiz = j + 1;
                }
            }
            fprintf(fp, "static CHA_DATA p%02d[%8d] = { // _size\n", i, arlen);
            for (j = 0; j < arsiz; j++) {
                if ((j % arpl) == 0) fprintf(fp, "        ");
                fprintf(fp, "%10lu", ulptr[j]);
                if (j < (arsiz - 1)) fprintf(fp, ",");
                if ((j % arpl) == (arpl - 1)) fprintf(fp, "\n");
            }
            if ((j % arpl) != 0) fprintf(fp, "\n");
            fprintf(fp, "};\n");
        } else if (i == _ivar) {
            arlen = cpsiz[i] / sizeof(int);
            arsiz = 1;
            for (j = 1; j < arlen; j++) {
                if (CHA_IVAR[j]) arsiz = j + 1;
            }
            fprintf(fp, "static CHA_DATA p%02d[%8d] = { // _ivar\n", i, arlen);
            for (j = 0; j < arsiz; j++) {
                if ((j % ivpl) == 0) fprintf(fp, "        ");
                fprintf(fp, "%10d", CHA_IVAR[j]);
                if (j < (arsiz - 1)) fprintf(fp, ",");
                if ((j % ivpl) == (ivpl - 1)) fprintf(fp, "\n");
            }
            if ((j % ivpl) != 0) fprintf(fp, "\n");
            fprintf(fp, "};\n");
        } else if (i == _dvar) {
            arlen = cpsiz[i] / sizeof(double);
            arsiz = 1;
            for (j = 1; j < arlen; j++) {
                if (CHA_DVAR[j]) arsiz = j + 1;
            }
            fprintf(fp, "static double   p%02d[%8d] = { // _dvar\n", i, arlen);
            for (j = 0; j < arsiz; j++) {
                if ((j % dvpl) == 0) fprintf(fp, "        ");
                fprintf(fp, "%15.9g", CHA_DVAR[j]);
                if (j < (arsiz - 1)) fprintf(fp, ",");
                if ((j % dvpl) == (dvpl - 1)) fprintf(fp, "\n");
            }
            if ((j % dvpl) != 0) fprintf(fp, "\n");
            fprintf(fp, "};\n");
        } else if (cpsiz[i] == 0) {
            fprintf(fp, "// empty array ->     p%02d\n", i);
        } else if ((cpsiz[i] % sizeof(long)) == 0) {
            arlen = cpsiz[i] / sizeof(long);
            arsiz = 0;
            ulptr = (unsigned long *) cp[i];
            if (ulptr) {
                for (j = 0; j < arlen; j++) {
                    if (ulptr[j]) arsiz = j + 1;
                }
            }
            if (arsiz < 2) {
                fprintf(fp, "static CHA_DATA p%02d[%8d] = {%10lu};\n",
                    i, arlen, ulptr[0]);
            } else {
                fprintf(fp, "static CHA_DATA p%02d[%8d] = {\n", i, arlen);
                for (j = 0; j < arsiz; j++) {
                    if ((j % arpl) == 0) fprintf(fp, "        ");
                    fprintf(fp, "0x%08lX", ulptr[j]);
                    if (j < (arsiz - 1)) fprintf(fp, ",");
                    if ((j % arpl) == (arpl - 1)) fprintf(fp, "\n");
                }
                if ((j % arpl) != 0) fprintf(fp, "\n");
                fprintf(fp, "};\n");
            }
        } else if ((cpsiz[i] % sizeof(short)) == 0) {
            arlen = cpsiz[i] / sizeof(short);
            fprintf(fp, "// NOTE: Only zero data implemented unless size%%4==0.\n");
            fprintf(fp, "static unsigned short p%02d[%8d] = {0};\n", i, arlen);
        } else {
            arlen = cpsiz[i];
            fprintf(fp, "// NOTEL Only zero data implemented unless size%%4==0.\n");
            fprintf(fp, "static unsigned char  p%02d[%8d] = {0};\n", i, arlen);
        }
    }
    fprintf(fp, "\n");
    // initialize ptr
    if (ptsiz < 1) {
        fprintf(fp, "static CHA_DATA *cha_data[NPTR] = {0};\n");
    } else {
        fprintf(fp, "static CHA_DATA *cha_data[NPTR] = {\n");
        fprintf(fp, "    ");
        for (i = 0; i < _reserve; i++) {
            fprintf(fp, "(CHA_DATA *)p%02d,", i);
        }
        fprintf(fp, "\n");
        for (i = _reserve; i < ptsiz; i++) {
            j = i - _reserve;
            if ((j % ptpl) == 0) fprintf(fp, "    ");
            if (cpsiz[i] == 0) {
                fprintf(fp, "NULL");
            } else {
                fprintf(fp, " p%02d", i);
            }
           if (i < (ptsiz - 1)) fprintf(fp, ",");
           if ((j % ptpl) == (ptpl - 1)) fprintf(fp, "\n");
        }
        if ((j % ptpl) != (ptpl - 1)) fprintf(fp, "\n");
        fprintf(fp, "};\n");
    }
    // print trailer
    for (i = 0; i < tlsz; i++) {
        fprintf(fp, "%s\n", tail[i]);
    }
    fclose(fp);

    return (0);
};
