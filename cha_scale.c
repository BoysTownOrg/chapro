// cha_scale.c - apply scale factor to chunk (in place)

#include "chapro.h"

FUNC(void) 
cha_scale(float *x, int cs, float scale)
{
    int k;

    for (k = 0; k < cs; k++) {
        x[k] *= scale;
    }
};
