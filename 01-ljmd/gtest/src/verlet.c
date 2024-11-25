/*
 * 1-d velocity verlet in reduced units
 */

#include "ljmd.h"

/* a few physical constants */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* velocity verlet x-component */

void verlet_x_1(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
    }
}

void verlet_x_2(mdsys_t *sys)
{
    int i;

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
    }
}

/* velocity verlet y-component */

void verlet_y_1(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->ry[i] += sys->dt*sys->vy[i];
    }
}

void verlet_y_2(mdsys_t *sys)
{
    int i;

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
    }
}

/* velocity verlet z-component */

void verlet_z_1(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rz[i] += sys->dt*sys->vz[i];
    }
}

void verlet_z_2(mdsys_t *sys)
{
    int i;

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}
