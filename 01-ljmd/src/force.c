#include <math.h>
#include "mdsys.h"
#include "utilities.h"
#include "force.h"
#if defined(_MPI)
    #include <mpi.h>
#endif

/* compute forces */
void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    # if defined(_MPI)
        int rank,size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    #endif

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    // MPI instruction - set auxiliary data to 0
    #if defined(_MPI)
        double epot = 0.0;
        azzero(sys->cx,sys->natoms);
        azzero(sys->cy,sys->natoms);
        azzero(sys->cz,sys->natoms);
    #endif

    // MPI instruction - broadcats particles positions from process of rank 0 to other ranks
    #if defined(_MPI)
        MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif

    int ii;

    #if defined(_MPI)
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
    #else
        for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);

                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

                sys->fx[i] += rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
            }
        }
    }
    #endif
}



