#include <math.h>
#include "mdsys.h"
#include "utilities.h"
#include "verlet.h"


void velverlet1(mdsys_t *sys)
{
    int i;
    double mvsq2e_inv = 1.0/mvsq2e;
    double mass_inv = 1.0/(sys->mass);
    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        //sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        //sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        //sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->vx[i] += 0.5*sys->dt * mvsq2e_inv * sys->fx[i] * mass_inv;
        sys->vy[i] += 0.5*sys->dt * mvsq2e_inv * sys->fy[i] * mass_inv;
        sys->vz[i] += 0.5*sys->dt * mvsq2e_inv * sys->fz[i] * mass_inv;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }
}

void velverlet2(mdsys_t *sys)
{
    int i;
    double mvsq2e_inv = 1.0/mvsq2e;
    double mass_inv = 1.0/(sys->mass);
    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        //sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        //sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        //sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->vx[i] += 0.5*sys->dt * mvsq2e_inv * sys->fx[i] * mass_inv;
        sys->vy[i] += 0.5*sys->dt * mvsq2e_inv * sys->fy[i] * mass_inv;
        sys->vz[i] += 0.5*sys->dt * mvsq2e_inv * sys->fz[i] * mass_inv;
    }
}



