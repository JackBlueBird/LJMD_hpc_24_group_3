#ifndef WALLCLOCK_FUNC
#define WALLCLOCK_FUNC
#include <sys/time.h>
#ifdef __cplusplus
extern "C"
{
#endif
extern double wallclock();
#ifdef __cplusplus
}
#endif

#endif

#ifndef AZZERO_FUNC
#define AZZERO_FUNC
#ifdef __cplusplus
extern "C"
{
#endif
extern void azzero(double *d, const int n);
#ifdef __cplusplus
}
#endif
#endif

#ifndef PBC_FUNC
#define PBC_FUNC
#ifdef __cplusplus
extern "C"
{
#endif
extern double pbc(double x, const double boxby2);
#ifdef __cplusplus
}
#endif
#endif
