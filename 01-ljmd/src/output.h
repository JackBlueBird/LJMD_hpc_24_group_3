#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include <string.h>

#include "mdsys.h"

#ifdef __cplusplus
extern "C"
{
#endif
extern void output(mdsys_t *sys, FILE *erg, FILE *traj);
#ifdef __cplusplus
}
#endif

#endif
