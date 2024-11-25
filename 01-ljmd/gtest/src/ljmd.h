/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#ifndef LJMD_H
#define LJMD_H

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

/* velocity verlet */
extern void verlet_x_1(mdsys_t *sys);
extern void verlet_x_2(mdsys_t *sys);
extern void verlet_y_1(mdsys_t *sys);
extern void verlet_y_2(mdsys_t *sys);
extern void verlet_z_1(mdsys_t *sys);
extern void verlet_z_2(mdsys_t *sys);

/* helper functions */
extern void azzero(double *d, const int n);
extern double wallclock();
extern void doublesleep(double t);

/* Reader functions */
extern int get_a_line(FILE *fp, char *buf);

/* compute kinetic energy */
extern void ekin(mdsys_t *sys);

/* helper function: apply minimum image convention */
extern double pbc(double x, const double boxby2);
/* compute forces */
extern void force(mdsys_t *sys);

#ifdef __cplusplus
}
#endif

#endif
