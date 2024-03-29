// cha_core.c - CHAPRO core functions

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "chapro.h"
#include "version.h"

#define free_null(p)    if(p){free(p);p=NULL;}

/***********************************************************/

FUNC(char *) 
cha_version(void)
{
    return (VER);                  
}

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
}

FUNC(void) 
cha_chunk_size(CHA_PTR cp, int cs)
{
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
}

FUNC(void *) 
cha_allocate(CHA_PTR cp, int cnt, int siz, int idx)
{
    cha_prepare(cp);
    assert(idx < NPTR);
    free_null(cp[idx]);
    cp[idx] = calloc(cnt, siz);
    ((int *)cp[_size])[idx] = cnt * siz;

    return (cp[idx]);
}

FUNC(void) 
cha_cleanup(CHA_PTR cp)
{
    int i;

    for (i = 0; i < NPTR; i++) {
        free_null(cp[i]);
    }
}

/***********************************************************/

static int32_t *
data_plan(CHA_PTR cp, int *psz, int *asz)
{
    int i, ptsiz, arsiz;
    int32_t *cpsiz;

    // count pointers
    ptsiz = 0;
    for (i = 0; i < NPTR; i++) {
        if (cp[i]) ptsiz = i + 1;
    }
    // sum array sizes
    cpsiz = (int32_t *) cp[_size]; //altered 11/17/2021 from (int *) to (int32_t *) to better match data type and avoid compiler warnings
    arsiz = 0;
    if (cpsiz) {
        for (i = 0; i < NPTR; i++) {
            if (cpsiz[i] < 0) {
                arsiz = 0;
		break;
	    }
            arsiz += cpsiz[i];
        }
    }
    *psz = ptsiz;
    *asz = arsiz;
    return(cpsiz);
}

static void
state_make(CHA_STA *state, int32_t *cpsiz, void **cp, int ptsiz, int arsiz)
{
    char *data;
    int i;
    void *mkdata;
    void **mkcp;

    // allocate memory
    mkcp = (void **)calloc(NPTR, sizeof(void *));
    mkdata = calloc(arsiz, 1);
    // copy pointers
    memcpy(mkcp, cp, NPTR * sizeof(void *));
    // copy data
    data = (char *)mkdata;
    for (i = 0; i < NPTR; i++) {
        if (cp[i]) {
            memcpy(data, cp[i], cpsiz[i]);
            mkcp[i] = data;
            data += cpsiz[i];
        }
    }
    // return state variables
    state->ptsiz = ptsiz;
    state->arsiz = arsiz;
    state->cp    = mkcp;
    state->data  = mkdata;
    state->sr    = CHA_DVAR[_fs] * 1000;
    state->cs    = CHA_IVAR[_cs];
}

/***********************************************************/

FUNC(int)
cha_data_gen(CHA_PTR cp, char *fn)
{
    int ptsiz, arsiz, arlen, i, j = 0;
    int32_t *cpsiz;
    uint32_t *ulptr;
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
    static int ptpl = 15;
    static int szpl = 8;
    static int ivpl = 8;
    static int dvpl = 4;
    static int arpl = 6;
    static int cptr = 4; // number of cast pointers

    cpsiz = data_plan(cp, &ptsiz, &arsiz);
    if (arsiz == 0) {
        return (1);
    }
    // open file
    fp = fopen(fn, "wt");
    if (fp == NULL) {
        return (2);
    }
    // print header
    fprintf(fp, "// cha_data.h - data size = %d bytes\n", arsiz);
    for (i = 0; i < hdsz; i++) {
        fprintf(fp, "%s\n", head[i]);
    }
    // initialize magic number
    if (ptsiz > 0) {
        fprintf(fp, "static CHA_DATA cha_magic[4] = ");
        fprintf(fp, "{0x55530,0x68131,%d,%d};\n", ptsiz, arsiz);
    }
    // initialize ptr arrays
    for (i = 0; i < ptsiz; i++) {
        if (i == _size) {
            arlen = cpsiz[i] / sizeof(int32_t);
            arsiz = 0;
            ulptr = (uint32_t *) cp[i];
            if (ulptr) {
                for (j = 0; j < arlen; j++) {
                    if (ulptr[j]) arsiz = j + 1;
                }
            }
            fprintf(fp, "static CHA_DATA p%02d[%8d] = { // _size\n", i, arlen);
            for (j = 0; j < arsiz; j++) {
                if ((j % szpl) == 0) fprintf(fp, "    ");
                fprintf(fp, "%8u", ulptr[j]);
                if (j < (arsiz - 1)) fprintf(fp, ",");
                if ((j % szpl) == (szpl - 1)) fprintf(fp, "\n");
            }
            if ((j % szpl) != 0) fprintf(fp, "\n");
            fprintf(fp, "};\n");
        } else if (i == _ivar) {
            arlen = cpsiz[i] / sizeof(int);
            arsiz = 1;
            for (j = 1; j < arlen; j++) {
                if (CHA_IVAR[j]) arsiz = j + 1;
            }
            fprintf(fp, "static CHA_DATA p%02d[%8d] = { // _ivar\n", i, arlen);
            for (j = 0; j < arsiz; j++) {
                if ((j % ivpl) == 0) fprintf(fp, "    ");
                fprintf(fp, "%8d", CHA_IVAR[j]);
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
                if ((j % dvpl) == 0) fprintf(fp, "    ");
                fprintf(fp, "%15.9g", CHA_DVAR[j]);
                if (j < (arsiz - 1)) fprintf(fp, ",");
                if ((j % dvpl) == (dvpl - 1)) fprintf(fp, "\n");
            }
            if ((j % dvpl) != 0) fprintf(fp, "\n");
            fprintf(fp, "};\n");
        } else if (cpsiz[i] == 0) {
            fprintf(fp, "// empty array: p%02d\n", i);
        } else if ((cpsiz[i] % sizeof(int32_t)) == 0) {
            arlen = cpsiz[i] / sizeof(int32_t);
            arsiz = 0;
            ulptr = (uint32_t *) cp[i];
            if (ulptr) {
                for (j = 0; j < arlen; j++) {
                    if (ulptr[j]) arsiz = j + 1;
                }
            }
            if (arsiz < 2) {
                fprintf(fp, "static CHA_DATA p%02d[%8d] = {%10u};\n",
                    i, arlen, ulptr[0]);
            } else {
                fprintf(fp, "static CHA_DATA p%02d[%8d] = {\n", i, arlen);
                for (j = 0; j < arsiz; j++) {
                    if ((j % arpl) == 0) fprintf(fp, "    ");
                    fprintf(fp, "0x%08X", ulptr[j]);
                    if (j < (arsiz - 1)) fprintf(fp, ",");
                    if ((j % arpl) == (arpl - 1)) fprintf(fp, "\n");
                }
                if ((j % arpl) != 0) fprintf(fp, "\n");
                fprintf(fp, "};\n");
            }
        } else if ((cpsiz[i] % sizeof(int16_t)) == 0) {
            arlen = cpsiz[i] / sizeof(int16_t);
            fprintf(fp, "// NOTE: Only zero data implemented unless size%%4==0.\n");
            fprintf(fp, "static uint16_t p%02d[%8d] = {0};\n", i, arlen);
        } else {
            arlen = cpsiz[i];
            fprintf(fp, "// NOTEL Only zero data implemented unless size%%4==0.\n");
            fprintf(fp, "static unsigned char  p%02d[%8d] = {0};\n", i, arlen);
        }
    }
    fprintf(fp, "\n");
    // initialize ptr
    if (ptsiz < 1) {
        fprintf(fp, "static CHA_DATA *cha_data[1] = {cha_magic};\n");
    } else {
        fprintf(fp, "static CHA_DATA *cha_data[%d] = {\n", NPTR);
        fprintf(fp, "    ");
        for (i = 0; i < cptr; i++) {
            fprintf(fp, "(CHA_DATA *)p%02d,", i);
        }
        fprintf(fp, "\n");
        for (i = cptr; i < ptsiz; i++) {
            j = i - cptr;
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
}

FUNC(int)
cha_data_save(CHA_PTR cp, char *fn)
{
    int ptsiz, arsiz, dtsiz, i;
    int32_t *cpsiz;
    CHA_DATA cha_magic[4];
    FILE *fp;

    cpsiz = data_plan(cp, &ptsiz, &arsiz);
    if (arsiz == 0) {
        return (1);
    }
    // open file
    fp = fopen(fn, "wb");
    if (fp == NULL) {
        return (2);
    }
    dtsiz = sizeof(CHA_DATA);
    // write magic number
    cha_magic[0] = 0x55530;
    cha_magic[1] = 0x68131;
    cha_magic[2] = ptsiz;
    cha_magic[3] = arsiz;
    fwrite(cha_magic, dtsiz, 4, fp);
    // write ptr arrays
    for (i = 0; i < ptsiz; i++) {
        if (cpsiz[i]) {
            arsiz = cpsiz[i] / dtsiz;
            fwrite(cp[i], dtsiz, arsiz, fp);
        }
    }
    fclose(fp);

    return (0);
}

FUNC(int)
cha_data_load(CHA_PTR cp, char *fn)
{
    int ptsiz, arsiz, dtsiz, i, rv;
    CHA_DATA cha_magic[4], file_magic[4], *file_size, *file_data;
    FILE *fp;

    (void) data_plan(cp, &ptsiz, &arsiz);
    if (arsiz == 0) {
        return (1);
    }
    // open file
    fp = fopen(fn, "rb");
    if (fp == NULL) {
        return (2);
    }
    dtsiz = sizeof(CHA_DATA);
    rv = 0; // default return value
    // read magic number
    cha_magic[0] = 0x55530;
    cha_magic[1] = 0x68131;
    cha_magic[2] = ptsiz;
    cha_magic[3] = arsiz;
    rv = fread(file_magic, dtsiz, 4, fp);
    // check magic number
    for (i = 0; i < 4; i++) {
        if (file_magic[i] != cha_magic[i]) {
            break;
        }
    }
    if (i < 4) {
        fclose(fp);
        return (4);
    }
    // read  & check size array
    ptsiz = file_magic[2];
    arsiz = NPTR * dtsiz;
    file_size = (CHA_DATA *) cha_allocate(cp, NPTR, dtsiz, 0);
    rv = fread(file_size, NPTR, dtsiz, fp);
    if (file_size[0] != arsiz) {
        return (4);
    }
    // read data arrays
    for (i = 1; i < ptsiz; i++) {
        if (file_size[i]) {
            arsiz = file_size[i] / dtsiz;
            file_data = (CHA_DATA *) cha_allocate(cp, arsiz, dtsiz, i);
            rv = fread(file_data, arsiz, dtsiz, fp);
        }
    }
    fclose(fp);

    return (rv);
}

FUNC(int)
cha_state_save(CHA_PTR cp, CHA_STA *state)
{
    int ptsiz, arsiz;
    int32_t *cpsiz;

    cpsiz = data_plan(cp, &ptsiz, &arsiz);
    if (arsiz == 0) {
        return (1);
    }
    state_make(state, cpsiz, cp, ptsiz, arsiz);

    return (0);
}

FUNC(int)
cha_state_copy(CHA_STA *new_state, CHA_STA *old_state)
{
    int ptsiz, arsiz;
    int32_t *cpsiz;
    void **cp;

    cp = old_state->cp;
    cpsiz = (int32_t *) cp[_size]; //altered 11/17/2021 from (int *) to (int32_t *) to better match data type and avoid compiler warnings
    ptsiz = old_state->ptsiz;
    arsiz = old_state->arsiz;
    state_make(new_state, cpsiz, cp, ptsiz, arsiz);
    return (0);
}

FUNC(int)
cha_state_free(CHA_STA *state)
{
    free(state->cp);
    free(state->data);
    state->cp = NULL;
    state->data = NULL;
    return (0);
}

