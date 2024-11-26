#include <math.h>
#include "mdsys.h"
#include "utilities.h"
#include "force.h"
#if defined(_MPI)
    #include <mpi.h>
#endif

/* helper function: apply minimum image convention */
static inline double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute forces */
void force(mdsys_t *sys)
{
    double r,ffac;
    double rsq;
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
        int ii;
        double epot = 0.0;
        azzero(sys->cx,sys->natoms);
        azzero(sys->cy,sys->natoms);
        azzero(sys->cz,sys->natoms);

        // MPI instruction - broadcats particles positions from process of rank 0 to other ranks
        MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // OPT instruction - define constants before the loop
        double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
        double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
        double rcsq = sys->rcut * sys->rcut;

        for(i=0; i < (sys->natoms); i+=size) {
            ii = i +rank;
            if (ii > sys->natoms) {
                break;
            }
            for(j=ii+1; j < (sys->natoms); ++j) {

                /* particles have no interactions with themselves */
                //if (ii==j) continue;

                /* get distance between particle i and j */
                rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
                ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
                rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
                // OPT instruction - work with r2 instead of r
                rsq = rx*rx + ry*ry + rz*rz;

                /* compute force and energy if within cutoff */
                // OPT instruction - comapre r2 with rc2 instead of r with rc
                if (rsq < rcsq) {
                    // OPT instruction - exploit the following mathematical proterties to perform forces computation
                    // nabla 1/r2 = -2*r/r4
                    // nabla 1/r6 = -6*r/r8
                    // babla 1/r12 = -12*r/r14

                    // OPT instruction - internal auxiliary variables
                    double r6, rinv;
                    rinv=1.0/rsq;
                    r6=rinv*rinv*rinv;
                
                    ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                    epot += r6*(c12*r6 - c6);
                    sys->cx[ii] += rx*ffac; sys->cx[j] -= rx*ffac;  
                    sys->cy[ii] += ry*ffac; sys->cy[j] -= ry*ffac;
                    sys->cz[ii] += rz*ffac; sys->cz[j] -= rz*ffac;
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
        // OPT instruction - define constants before the loop
        double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
        double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
        double rcsq = sys->rcut * sys->rcut;

        for(i=0; i < (sys->natoms)-1; ++i) {
            // OPT instruction - set j to start from i+1 to exploit Newton's Third Law
            for(j=i+1; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            // OPT instruction - work with r2 instead of r
            rsq = rx*rx + ry*ry + rz*rz;

            /* compute force and energy if within cutoff */
            // OPT instruction - comapre r2 with rc2 instead of r with rc
            if (rsq < rcsq) {
                // OPT instruction - exploit the following mathematical proterties to perform forces computation
                // nabla 1/r2 = -2*r/r4
                // nabla 1/r6 = -6*r/r8
                // babla 1/r12 = -12*r/r14

                // OPT instruction - internal auxiliary variables
                double r6, rinv;
                rinv=1.0/rsq;
                r6=rinv*rinv*rinv;
                
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                sys->epot += r6*(c12*r6 - c6);
                sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;   
                sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
                sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
            }
        }
    }
    #endif
}



