#ifndef UTILITIES_H
#define UTILITIES_H
#include <sys/time.h>

double wallclock();
void azzero(double *d, const int n);
double pbc(double x, const double boxby2);

#endif
