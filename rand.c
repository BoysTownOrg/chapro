// rand.c

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

// simple generator initializes MT generator
#define _SET(s)  (js=s)
#define _IUNI	(jz=js,js^=(js<<13),js^=(js>>17),js^=(js<<5),jz+js)
#define ulabs(x) ((uint32_t)labs(x))

static uint32_t js = 123456789;
static uint32_t jz = 0;

//-----------------------------------------------------------

// Mersenne Twister - based on public-domain source code [Oct-2009]
// http://www.qbrundage.com/michaelb/pubs/essays/random_number_generation.html

#define IUNI            mt_random()
#define SET(s)		mt_init(s)
#define MT_LEN		624
#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

static int mt_index = -1;
static uint32_t mt_buffer[MT_LEN];

void mt_init(uint32_t s) {
    int i;

    _SET(s);
    for (i = 0; i < MT_LEN; i++) {
        mt_buffer[i] = _IUNI;
    }
    mt_index = 0;
}

uint32_t mt_random() {
    uint32_t * b = mt_buffer;
    uint32_t s;
    int i, idx;
	
    if (mt_index < 0) {
	mt_init(js);
    }
    idx = mt_index;
    if (idx == MT_LEN * sizeof(uint32_t)) {
        idx = 0;
        i = 0;
        for (; i < MT_IB; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
        }
        for (; i < MT_LEN-1; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
        }
        s = TWIST(b, MT_LEN-1, 0);
        b[MT_LEN-1] = b[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }
    mt_index = idx + sizeof(uint32_t);
    return (*(uint32_t *)((unsigned char *)b + idx));
}

//-----------------------------------------------------------

#define UNI		((IUNI + 0.5) * (1.0 / 4294967296.0))
#define RNOR	(hz=IUNI, iz=hz&127, (ulabs(hz)<kn[iz])? hz*wn[iz] : nfix())

static double wn[128],fn[128];
static int32_t hz;
static uint32_t iz, kn[128];

/* zigset - sets the seed and creates the tables */

static void 
zigset()
{  
    const double m1 = 2147483648.0;
    double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
    int i;
    
    /* Set up tables for RNOR */

    q = vn/exp(-0.5 * dn * dn);
    kn[0] = (int) ((dn/q) * m1);
    kn[1] = 0;

    wn[0] = q/m1;
    wn[127] = dn/m1;

    fn[0]=1.;
    fn[127] = exp(-0.5 * dn * dn);

    for(i = 126;i >= 1;i--) {
	dn = (int) sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
	kn[i+1] = (int) ((dn / tn) * m1);
	tn = dn;
	fn[i] = exp(-.5*dn*dn);
	wn[i] = dn/m1;
    }

}

/* nfix - generates variates from the residue when rejection in RNOR occurs. */

static double 
nfix(void)
{
    const double r = 3.442620;     /* The start of the right tail */
    static double x, y;
    static int table_set = 0;

    if (!table_set) {
	zigset();
	table_set++;
    }

    for(;;)  {  
	x = hz * wn[iz];      /* iz==0, handles the base strip */
	if(iz == 0) {
	    do { 
		x=-log(UNI) * 0.2904764;    /* .2904764 is 1/r */
		y=-log(UNI);
	    } while((y + y) < (x * x));
	    return (hz > 0) ? (r + x) : (-r - x);
	}
      
	/* iz>0, handle the wedges of other strips */

	if((fn[iz] + UNI * (fn[iz - 1] - fn[iz])) < exp(-0.5 * x * x) ) 
	    return x;

     /* initiate, try to exit for(;;) for loop */

      hz = IUNI;
      iz = hz & 127;
      if(ulabs(hz)<kn[iz]) 
	  return (hz*wn[iz]);
  }

}

//-----------------------------------------------------------

FUNC(void) sp_rand(float *x, int n) 
{
    int i;

    for (i = 0; i < n; i++)
	x[i] = (float) UNI;
		
}

FUNC(void) sp_randn(float *x, int n) 
{
    int     i;

    for (i = 0; i < n; i++) {
	x[i] = (float) RNOR;
    }
}

FUNC(void) sp_randseed(uint32_t s) 
{
    SET(s);
}
