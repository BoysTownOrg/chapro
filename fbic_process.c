// fbic_process.c - filterbank+instantaneous_compression processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_fb.h"

/***********************************************************/

FUNC(void)
cha_fbic_process(CHA_PTR cp, float *x, float *y)
{
    float *z = (float *) cp[_xx];
    int n = CHA_IVAR[_cs];

    cha_filterbank_analyze(cp, x, z, n);
    cha_compressor_process(cp, z, z, n);
    cha_filterbank_synthesize(cp, z, y, n);
}
