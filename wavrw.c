// wavrw.c - read & write WAV-format files.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "sigpro.h"
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#define _access access
#define _close  close
#define _lseek  lseek
#define _open   open
#define _read   read
#define _write  write 
#define _hypot  hypot
#define _strdup strdup
#endif

#define limit(a,b,c)  (((b)<(a))?(a):((b)>(c))?(c):(b))
#define free_null(x)  if(x){free(x);x=NULL;}
#define MAXNCH        256
#ifndef O_BINARY
#define O_BINARY      0
#endif
#define PMODE         (S_IREAD|S_IWRITE)
#define AFLAG         (O_RDWR|O_BINARY)
#define OFLAG         (O_RDWR|O_CREAT|O_BINARY)
#define OPMOD         (O_RDONLY|O_BINARY)

/*..........................................................................*/

#define MINWAV 36

typedef struct {
    int16_t   format;
    int16_t   channels;
    int32_t    samp_sec;
    int32_t    byte_sec;
    int16_t   block_align;
    int16_t   bits_smp;
}       WAV_fmt;

static WAV_fmt wav_hdr; 
static int32_t     dat_pos;

static int32_t
read_long(int fd)
{
    int32_t n;

   (void) _read(fd, &n, sizeof(int32_t));

   return (n);
}

static int16_t
read_short(int fd)
{
    int16_t n;

   (void) _read(fd, &n, sizeof(int16_t));

   return (n);
}

static int
wav_head(const char *fn, VAR *vl, float *fs)
{
    char    id[MINWAV];
    int     fd, status;
    int32_t    fmtsz, datsz, smpsz, frmsz, fpos, size;

    if (vl == NULL) {
	return(-1);
    }
    vl[0].rows = 0;
    vl[0].cols = 0;
    fd = _open(fn, OPMOD);		// open the file
    if (fd < 0) {
	return (-2);
    }
    status = _read(fd, id, 4);
    if (strncmp(id, "RIFF", 4) != 0) {
	_close(fd);
	return (-3);
    }
    size = read_long(fd);
    status = _read(fd, id, 8);
    if (strncmp(id, "WAVEfmt ", 8) != 0) {
	_close(fd);
	return (-4);
    }
    fmtsz = read_long(fd);
    wav_hdr.format = read_short(fd);
    wav_hdr.channels = read_short(fd);
    wav_hdr.samp_sec = read_long(fd);
    wav_hdr.byte_sec = read_long(fd);
    wav_hdr.block_align = read_short(fd);
    wav_hdr.bits_smp = read_short(fd);
    smpsz = (wav_hdr.bits_smp + 7) / 8;
    frmsz = smpsz * wav_hdr.channels;
    datsz = 0;
    fpos = 20 + fmtsz;
    for (;;) {
	_lseek(fd, fpos, 0);
	status = _read(fd, id, 4);
	if (status != 4) {
	    break;
 	}
        size = read_long(fd);
        fpos += 8;
	if (strncmp(id, "data", 4) == 0) {
	    datsz = size;
	    break;
	}
        fpos += size;
    }
    dat_pos = fpos;
    *fs = (float) wav_hdr.samp_sec;
    vl[0].rows = datsz / frmsz;
    vl[0].cols = wav_hdr.channels;
    vl[0].last = 1;

    _close(fd);
    return (0);
}

static void
wav_load(const char *fn, VAR *vl, int *ifr, int *nfr)
{
    char    bv;
    float  *f, sclf;
    int     fd, i;
    int32_t    smpsz, nchn, ifrm, nfrm, nsmp, lv;

    if (vl == NULL) {
	return;
    }
    fd = _open(fn, OPMOD);		// open the file
    if (fd < 0) {
	return;
    }
    ifrm = ifr ? ifr[0] : 0;
    nfrm = vl[0].rows - ifrm;
    if (nfr && nfr[0] < nfrm) {
	nfrm = nfr[0];
    }
    nchn = vl[0].cols;
    nsmp = nfrm  * nchn;
    if (nsmp <= 0) {
	_close(fd);
	return;
    }
    smpsz = (wav_hdr.bits_smp + 7) / 8;
    f = (float *) calloc(nsmp, sizeof(float));
    sclf = (float) (1 / (pow(2.0, wav_hdr.bits_smp - 1.0) - 1));
    _lseek(fd, dat_pos + ifrm * nchn * smpsz, 0);
    for (i = 0; i < nsmp; i++) {
	if (smpsz == 1) {
	    (void) _read(fd, &bv, 1);
	    f[i] = (float) bv * sclf;
	} else if (smpsz == 2) {
	    f[i] = (float) read_short(fd) * sclf;
	} else if (smpsz == 3) {
	    lv = 0;
	    (void) _read(fd, &lv, 3);
	    f[i] = (float) lv * sclf;
	} else if (smpsz == 4) {
	    f[i] = (float) read_long(fd) * sclf;
	} else {
	    break;
	}
    }
    _close(fd);
    vl[0].name = _strdup(fn);
    vl[0].data = f;
    vl[0].dtyp = 8;
    vl[0].rows = nfrm;
}

/**************************************************************/

static void
write_data_1(int fd, VAR *vl)
{
    char *i1;
    double *f8, sc;
    float *f4;
    int j, m, n;
    int32_t *i4;
    int16_t *i2;
    unsigned char *u1;
    uint32_t *u4;
    uint32_t *u2;

    n = vl[0].rows * vl[0].cols;
    n = vl[0].cmpx ? 2 * n : n;
    // convert VAR data to char array
    i1 = (char *) calloc(n, sizeof(char));
    switch (vl[0].dtyp) {
    case 1:		// int*1
	memcpy(i1, vl[0].data, n * sizeof(char));
        break;
    case 2:		// uint*1
        u1 = (unsigned char *) vl[0].data;
        m = 1 << 7;
        for (j = 0; j < n; j++) {
    	    i1[j] = (char) (u1[j] - m);
        }
        break;
    case 3:		// int*2
        i2 = (int16_t *) vl[0].data;
        for (j = 0; j < n; j++) {
    	    i1[j] = (char) (i2[j] >> 8);
        }
        break;
    case 4:		// uint*2
        u2 = (uint32_t *) vl[0].data;
        m = 1 << 15;
        for (j = 0; j < n; j++) {
    	    i1[j] = (char) ((u2[j] - m) >> 8);
        }
        break;
    case 5:		// int*4
        i4 = (int32_t *) vl[0].data;
        for (j = 0; j < n; j++) {
    	    i1[j] = (char) (i4[j] >> 24);
        }
        break;
    case 6:		// uint*4
        u4 = (uint32_t *) vl[0].data;
        m = 1 << 31;
        for (j = 0; j < n; j++) {
    	    i1[j] = (char) ((u4[j] - m) >> 24);
        }
        break;
    case 7:		// float*4
    case 8:		// float*4
        f4 = (float *) vl[0].data;
	sc = pow(2.0, 7.0) - 1;
        for (j = 0; j < n; j++) {
	    i1[j] = (char) (limit(-1, f4[j], 1) * sc);
        }
        break;
    case 9:		// float*8
        f8 = (double *) vl[0].data;
	sc = pow(2.0, 7.0) - 1;
        for (j = 0; j < n; j++) {
	    i1[j] = (char) (limit(-1, f8[j], 1) * sc);
        }
        break;
    }
    // write long array to file
    _write(fd, i1, n * sizeof(char));
    free(i1);
}

static void
write_data_2(int fd, VAR *vl)
{
    char *i1;
    double *f8, sc;
    float *f4;
    int j, m, n;
    int32_t *i4;
    int16_t *i2;
    unsigned char *u1;
    uint32_t *u4;
    uint32_t *u2;

    n = vl[0].rows * vl[0].cols;
    n = vl[0].cmpx ? 2 * n : n;
    // convert VAR data to short array
    i2 = (int16_t *) calloc(n, sizeof(int16_t));
    switch (vl[0].dtyp) {
    case 1:		// int*1
        i1 = (char *) vl[0].data;
        for (j = 0; j < n; j++) {
    	    i2[j] = (int16_t) i1[j] << 8;
        }
        break;
    case 2:		// uint*1
        u1 = (unsigned char *) vl[0].data;
        m = 1 << 7;
        for (j = 0; j < n; j++) {
    	    i2[j] = (int16_t) (u1[j] - m) << 8;
        }
        break;
    case 3:		// int*2
	memcpy(i2, vl[0].data, n * sizeof(int16_t));
        break;
    case 4:		// uint*2
        u2 = (uint32_t *) vl[0].data;
        m = 1 << 15;
        for (j = 0; j < n; j++) {
    	    i2[j] = (int16_t) (u2[j] - m);
        }
        break;
    case 5:		// int*4
        i4 = (int32_t *) vl[0].data;
        for (j = 0; j < n; j++) {
    	    i2[j] = (int16_t) (i4[j] >> 16);
        }
        break;
    case 6:		// uint*4
        u4 = (uint32_t *) vl[0].data;
        m = 1 << 31;
        for (j = 0; j < n; j++) {
    	    i2[j] = (int16_t) ((u4[j] - m) >> 16);
        }
        break;
    case 8:		// float*4
        f4 = (float *) vl[0].data;
	sc = pow(2.0, 15.0) - 1;
        for (j = 0; j < n; j++) {
	    i2[j] = (int16_t) (limit(-1, f4[j], 1) * sc);
        }
        break;
    case 9:		// float*8
        f8 = (double *) vl[0].data;
	sc = pow(2.0, 15.0) - 1;
        for (j = 0; j < n; j++) {
	    i2[j] = (int16_t) (limit(-1, f8[j], 1) * sc);
        }
        break;
    }
    // write short array to file
    _write(fd, i2, n * sizeof(int16_t));
    free(i2);
}

static void
write_data_4(int fd, VAR *vl)
{
    char *i1;
    double *f8, sc;
    float *f4;
    int j, m, n;
    int32_t *i4;
    int16_t *i2;
    unsigned char *u1;
    uint32_t *u4;
    uint32_t *u2;

    n = vl[0].rows * vl[0].cols;
    n = vl[0].cmpx ? 2 * n : n;
    // convert VAR data to long array
    i4 = (int32_t *) calloc(n, sizeof(int32_t));
    switch (vl[0].dtyp) {
    case 1:		// int*1
        i1 = (char *) vl[0].data;
        for (j = 0; j < n; j++) {
    	    i4[j] = (int32_t) i1[j] << 24;
        }
        break;
    case 2:		// uint*1
        u1 = (unsigned char *) vl[0].data;
        m = 1 << 7;
        for (j = 0; j < n; j++) {
    	    i4[j] = (int32_t) (u1[j] - m) << 24;
        }
        break;
    case 3:		// int*2
        i2 = (int16_t *) vl[0].data;
        for (j = 0; j < n; j++) {
    	    i4[j] = (int32_t) i2[j] << 16;
        }
        break;
    case 4:		// uint*2
        u2 = (uint32_t *) vl[0].data;
        m = 1 << 15;
        for (j = 0; j < n; j++) {
    	    i4[j] = (int32_t) (u2[j] - m) << 16;
        }
        break;
    case 5:		// int*4
	memcpy(i4, vl[0].data, n * sizeof(int32_t));
        break;
    case 6:		// uint*4
        u4 = (uint32_t *) vl[0].data;
        m = 1 << 31;
        for (j = 0; j < n; j++) {
    	    i4[j] = (int32_t) (u4[j] - m);
        }
        break;
    case 7:		// float*4
    case 8:		// float*4
        f4 = (float *) vl[0].data;
	sc = pow(2.0, 31.0) - 1;
        for (j = 0; j < n; j++) {
	    i4[j] = (int32_t) (limit(-1, f4[j], 1) * sc);
        }
        break;
    case 9:		// float*8
        f8 = (double *) vl[0].data;
	sc = pow(2.0, 31.0) - 1;
        for (j = 0; j < n; j++) {
	    i4[j] = (int32_t) (limit(-1, f8[j], 1) * sc);
        }
        break;
    }
    // write long array to file
    _write(fd, i4, n * sizeof(int32_t));
    free(i4);
}

static void
write_data(int fd, VAR *vl, int nbits)
{
    if (nbits <= 8) {
	write_data_1(fd, vl);
    } else if (nbits <= 16) {
	write_data_2(fd, vl);
    } else {
	write_data_4(fd, vl);
   }
}

static int
save_wav(const char *fn, VAR *vl, float *fs, int nbits, int nreps)
{
    char    z = 0;
    int     fd, i, nbps, nchn, nsmp;
    int32_t    rifsz, fmtsz, datsz;
    WAV_fmt wfmt;

    fd = _open(fn, OFLAG, PMODE);		// open the file
    if (fd < 0) {
	return (fd);
    }
    nsmp = vl[0].rows;
    nchn = vl[0].cols;
    nbps = (nbits + 7) / 8;
    wfmt.format = 1;		/* PCM */
    wfmt.channels = nchn;
    wfmt.samp_sec = (int32_t) (*fs);
    wfmt.block_align = nbps * nchn;
    wfmt.bits_smp = nbits;
    wfmt.byte_sec = wfmt.samp_sec * nchn * nbps;
    datsz = nsmp * nchn * nbps * nreps;
    fmtsz = sizeof(wfmt);
    rifsz = 5 * 4 + fmtsz + datsz + (datsz % 2);
    _write(fd, "RIFF", 4);
    _write(fd, (char *) &rifsz, 4);
    _write(fd, "WAVE", 4);
    _write(fd, "fmt ", 4);
    _write(fd, (char *) &fmtsz, 4);
    _write(fd, (char *) &wfmt, (unsigned) fmtsz);
    _write(fd, "data", 4);
    _write(fd, (char *) &datsz, 4);
    for (i = 0; i < nreps; i++) {
        write_data(fd, vl, nbits);
    }
    if (datsz % 2) {
        _write(fd, (char *) &z, 1);
    }
    _close(fd);
    return (0);
}

/*..........................................................................*/

FUNC(VAR *) sp_wav_info(
    const char  *fn,	// file name
    float *fs		// sampling rate
)
{
    VAR *vl;

    vl = sp_var_alloc(1);
    if (wav_head(fn, vl, fs) < 0) {
	sp_var_clear(vl);
	return (NULL);
    }

    return (vl);
}

FUNC(VAR *) sp_wav_read(
    const char  *fn,	// file name
    int   *ifr,		// initial sample & channel
    int   *nfr,		// number of samples & channels
    float *fs		// sampling rate
)
{
    VAR *vl;

    vl = sp_var_alloc(1);
    if (wav_head(fn, vl, fs) < 0) {
	sp_var_clear(vl);
	return (NULL);
    }
    wav_load(fn, vl, ifr, nfr);

    return (vl);
}

FUNC(int) sp_wav_write(
    const char  *fn,	// file name
    VAR   *vl,		// variable list
    float *fs,		// sampling rate
    int    nbits        // number of bits
)
{
    return save_wav(fn, vl, fs, nbits, 1);
}

FUNC(int) sp_wav_write_rep(
    const char  *fn,	// file name
    VAR   *vl,		// variable list
    float *fs,		// sampling rate
    int    nbits,       // number of bits
    int    nreps        // number of repetitions
)
{
    return save_wav(fn, vl, fs, nbits, nreps);
}

/*..........................................................................*/
