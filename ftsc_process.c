// ftsc_process.c - stft+agc processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

#ifndef ARDUINO
#include "cha_ft.h"
#endif

/***********************************************************/

FUNC(void)
cha_ftsc_process(CHA_PTR cp, float *x, float *y)
{
    float *z = (float *) cp[_xx];
    int n = CHA_IVAR[_cs];

    cha_stft_analyze(cp, x, z, n);
    cha_compressor_process(cp, z, z, n);
    cha_stft_synthesize(cp, z, y, n);
}
