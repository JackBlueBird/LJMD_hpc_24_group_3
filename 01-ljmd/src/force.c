#include <math.h>
#include "mdsys.h"
#include "utilities.h"
#include "force.h"
#include <mpi.h>

/* compute forces */
void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    // MPI instruction - set auxiliary data to 0
    double epot = 0.0;
    azzero(sys->cx,sys->natoms);
    azzero(sys->cy,sys->natoms);
    azzero(sys->cz,sys->natoms);

    // MPI instruction - broadcats particles positions from process of rank 0 to other ranks
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int ii;

    for(i=0; i < (sys->natoms); i+=size) {
        ii = i +rank;
        if (ii > sys->natoms) {
            break;
        }
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (ii==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[ii] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[ii] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[ii] - sys->rz[j], 0.5*sys->box);
            // this can be avoid since we are using even powers of r
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {

                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);
                // MPI instruction store results in auxiliary data store
                epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

                sys->cx[ii] += rx/r*ffac;
                sys->cy[ii] += ry/r*ffac;
                sys->cz[ii] += rz/r*ffac;
            }
        }
    }
    // MPI instruction - sum forces contribution from all process to process of rank 0
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // MPI instruction - sum potential energy contribution from all ranks to process of rank 0
    MPI_Reduce(&epot, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/* velocity verlet */
void velverlet(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }

    /* compute forces and potential energy */
    force(sys);

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}
