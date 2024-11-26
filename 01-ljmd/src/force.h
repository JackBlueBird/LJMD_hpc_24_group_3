#ifndef PBC_FUNC
#define PBC_FUNC
#ifdef __cplusplus
extern "C"
{
#endif
static inline double pbc(double x, const double boxby2);
#ifdef __cplusplus
}
#endif
#endif

#ifndef FORCE_H
#define FORCE_H
#include "mdsys.h"
#include "utilities.h"
#ifdef __cplusplus
extern "C"
{
#endif
extern void force(mdsys_t *sys);
#ifdef __cplusplus
}
#endif
#endif
