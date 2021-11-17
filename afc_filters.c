// afc_filters.c - adaptive-feedback-cancelation filter-access support

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

/***********************************************************/

FUNC(int)
cha_afc_filters(CHA_PTR cp, CHA_AFC *afc)
{
    if (afc->pcp != cp) {
        // copy quality-metric info info to CHA_AFC
        afc->qm = (float *) cp[_qm];
        afc->iqm = 0;
        afc->iqmp = (int32_t *) cp[_iqmp];  //updated WEA 11/17/2021 from (int *) to better match the data type (to avoid compiler warning)
        // copy filters info to CHA_AFC
        afc->fbl = CHA_IVAR[_fbl];
        afc->wfl = CHA_IVAR[_wfl];
        afc->pfl = CHA_IVAR[_pfl];
        afc->efbp = (float *) cp[_efbp];
        afc->sfbp = (float *) cp[_sfbp];
        afc->wfrp = (float *) cp[_wfrp];
        afc->ffrp = (float *) cp[_ffrp];
        // initialize afc_process
        CHA_IVAR[_in1] = 0;
        afc->pcp = cp;
   }

    return (0);
}
