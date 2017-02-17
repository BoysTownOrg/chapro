// chapro.h - function prototypes for CHA common functions
#ifndef CHAPRO_H
#define CHAPRO_H

#ifdef DLL
#define FUNC(type) __declspec(dllexport) type _stdcall
#else
#define FUNC(type) type
#endif

#define NPTR       64
#define NVAR       16
#define CHA_IVAR   ((int *)cp[_ivar])
#define CHA_DVAR   ((double *)cp[_dvar])

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2           0.693147180559945309417
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2       0.707106781186547524401
#endif
#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880 
#endif
#ifndef WIN32
#define _hypot          hypot
#define __inline        inline
#endif

#define fmin(x,y)       ((x<y)?(x):(y))
#define fmove(x,y,n)    memmove(x,y,(n)*sizeof(float))
#define fcopy(x,y,n)    memcpy(x,y,(n)*sizeof(float))
#define fzero(x,n)      memset(x,0,(n)*sizeof(float))
#define round(x)        ((int)floorf((x)+0.5))
#define log2(x)         (logf(x)/M_LN2)

typedef unsigned long CHA_DATA;
typedef unsigned long *CHA_LPTR;
typedef void **CHA_PTR;

/*****************************************************/

// CHA common functions

FUNC(void *) cha_allocate(CHA_PTR, int, int, int);
FUNC(void)   cha_cleanup(CHA_PTR);
FUNC(int)    cha_data_gen(CHA_PTR, char *);
FUNC(float)  cha_db1(float);
FUNC(float)  cha_db2(float);
FUNC(void)   cha_fft_cr(float *, int);
FUNC(void)   cha_fft_rc(float *, int);
FUNC(void)   cha_fft(float *, int);
FUNC(void)   cha_ifft(float *, int);
FUNC(void)   cha_prepare(CHA_PTR);
FUNC(void)   cha_scale(float *, int, float);
FUNC(float)  cha_undb1(float);
FUNC(float)  cha_undb2(float);
FUNC(char *) cha_version(void);

/*****************************************************/

// reserved pointer indices

#define _size     0
#define _ivar     1
#define _dvar     2
#define _reserve  3

#endif /* CHAPRO_H */
