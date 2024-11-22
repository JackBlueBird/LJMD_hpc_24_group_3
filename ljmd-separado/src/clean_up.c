#include "mdsys.h"
#include <stdlib.h>

void clean_up(mdsys_t *sys)
{
    free(sys->rx);
    free(sys->ry);
    free(sys->rz);
    free(sys->vx);
    free(sys->vy);
    free(sys->vz);
    free(sys->fx);
    free(sys->fy);
    free(sys->fz);
}
