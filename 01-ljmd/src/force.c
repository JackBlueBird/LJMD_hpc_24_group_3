#include <math.h>
#include "mdsys.h"
#include "utilities.h"
#include "force.h"
#include <omp.h>
#include <mpi.h>

/* helper function: apply minimum image convention */
// OPT - static inline to help compiler otimization
static inline double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute forces */
/*
void force(mdsys_t *sys)
{
    double r,ffac;
    double rsq;
    double rx,ry,rz;
    int i,j;


    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    // OPT instruction - define constants before the loop
    double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;

    for(i=0; i < (sys->natoms)-1; ++i) {
        // OPT instruction - set j to start from i+1 to exploit Newton's Third Law
        for(j=i+1; j < (sys->natoms); ++j) {
            // particles have no interactions with themselves 

            // get distance between particle i and j 
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            // OPT instruction - work with r2 instead of r
            rsq = rx*rx + ry*ry + rz*rz;

            // compute force and energy if within cutoff 
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
}
*/


/* compute forces */
void force(mdsys_t *sys) {

    // OMP set up
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_num_threads();
    #endif

    // MPI set up
    //# if defined(_MPI)
        int rank,size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    //#endif

    double sigma, sigma6;
    sigma=sys->sigma;
    sigma6=sigma*sigma*sigma*sigma*sigma*sigma;

    double c12 = 4.0 * sys->epsilon * sigma6*sigma6;//pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * sigma6;//pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    // OMP - thread potential energy
    double epot_local = 0.0;

    double r,ffac;
    double rsq;
    double rx,ry,rz;
    double r6, rinv; //OPT
    int i,j;

    // Set zero energy and forces
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);


    //# if defined(_MPI)
        // OMP and MPI - Arrays of support cx -> helper for fx
        // gz: put a flag in the future
        int ii;
        // MPI - rank potential energy
        double epot = 0.0;
        azzero(sys->cx, nthreads * sys->natoms);
        azzero(sys->cy, nthreads * sys->natoms);
        azzero(sys->cz, nthreads * sys->natoms);

        // MPI instruction - broadcats particles positions from process of rank 0 to other ranks
        MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        #ifdef _OPENMP
        #pragma omp parallel num_threads(sys->nthreads)
        #endif
        {   
            // OMP - set number of threads
            int tid=0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif

            // OMP - define pointers and shift them
            double *cx = sys->cx + tid * sys->natoms;
            double *cy = sys->cy + tid * sys->natoms;
            double *cz = sys->cz + tid * sys->natoms;

            // OMP - divide the work between threads
            int grids = nthreads*size;

            // OPT - Newton thir's law
            for (int i = 0; i < sys->natoms - 1; i += grids) {
                int ii = i + tid + sys->nthreads*rank;
                if (ii >= sys->natoms - 1) {
                    break;
                }
                // perhaps ask later
                double rx1 = sys->rx[ii];
                double ry1 = sys->ry[ii];
                double rz1 = sys->rz[ii];

                for (int j = ii + 1; j < sys->natoms; ++j) {
                    // distances bettween particles
                    rx = pbc(rx1 - sys->rx[j], 0.5 * sys->box);
                    ry = pbc(ry1 - sys->ry[j], 0.5 * sys->box);
                    rz = pbc(rz1 - sys->rz[j], 0.5 * sys->box);
                    rsq = rx * rx + ry * ry + rz * rz;

                    // If is inside of the cutoff then compute forces and energies
                    if (rsq < rcsq) {
                        // OPT - populate auxiliary variables
                        double rsqinv = 1.0 / rsq;
                        double r6 = rsqinv * rsqinv * rsqinv;
                        ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rsqinv;

                        // OMP - protect writing of this variable from race conditions
                        #ifdef _OPENMP
                        #pragma omp atomic
                        #endif
                        epot_local += r6 * (c12 * r6 - c6);

                        // OPT - Newton thir's law
                        cx[ii] += rx * ffac;
                        cy[ii] += ry * ffac;
                        cz[ii] += rz * ffac;

                        cx[j] -= rx * ffac;
                        cy[j] -= ry * ffac;
                        cz[j] -= rz * ffac;
                    }
                }
            } // END FOR LOOP

            // OMP - sincronize 
            #ifdef _OPENMP
            #pragma omp barrier
            #endif

            int chunk_size = (sys->natoms + sys->nthreads - 1) / sys->nthreads;
            int start = tid * chunk_size;
            int end = (start + chunk_size > sys->natoms) ? sys->natoms : start + chunk_size;

            // OMP - combine the local  forces into global forces 
            for (int t = 1; t < sys->nthreads; ++t) {
                int offset = t * sys->natoms;
                for (int i = start; i < end; ++i) {
                    sys->cx[i] += sys->cx[offset + i];
                    sys->cy[i] += sys->cy[offset + i];
                    sys->cz[i] += sys->cz[offset + i];
                }
            }
            // OMP + MPI - potential energy is stored in a distributed variable for each process
            epot = epot_local;
        }
        // end OPENMPI


        
        // MPI instruction - sum forces contribution from all process to process of rank 0
        MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        // MPI instruction - sum potential energy contribution from all ranks to process of rank 0
        MPI_Reduce(&epot, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //# else
    //# endif


    /*
    #ifdef _OPENMP
    #pragma omp parallel num_threads(sys->nthreads)
    #endif
    {
        int tid=0;
        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #endif

        // allocated as natoms*nthreads
        double *cx = sys->cx + tid * sys->natoms;
        double *cy = sys->cy + tid * sys->natoms;
        double *cz = sys->cz + tid * sys->natoms;

        // OPT instruction - declare auxiliary variables
        double rx, ry, rz, rsq, ffac;

        // OMP instruction - divide the work between threads
        int grids = sys->nthreads;

        for (int i = 0; i < sys->natoms - 1; i += grids) {
            int ii = i + tid;
            if (ii >= sys->natoms - 1) break;
                double rx1 = sys->rx[ii];
                double ry1 = sys->ry[ii];
                double rz1 = sys->rz[ii];
            for (int j = ii + 1; j < sys->natoms; ++j) {
                // distances bettween particles
                rx = pbc(rx1 - sys->rx[j], 0.5 * sys->box);
                ry = pbc(ry1 - sys->ry[j], 0.5 * sys->box);
                rz = pbc(rz1 - sys->rz[j], 0.5 * sys->box);
                rsq = rx * rx + ry * ry + rz * rz;

                // If is inside of the cutoff then compute forces and energies
                if (rsq < rcsq) {
                    // OPT instruction - exploit the following mathematical proterties to perform forces computation
                    // nabla 1/r2 = -2*r/r4
                    // nabla 1/r6 = -6*r/r8
                    // nabla 1/r12 = -12*r/r14

                    // OPT instruction - internal auxiliary variables
                    double rsqinv = 1.0 / rsq;
                    double r6 = rsqinv * rsqinv * rsqinv;
                    ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rsqinv;

                    // OMP instruction - protect sum from multiple access
                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    epot_local += r6 * (c12 * r6 - c6);

                    cx[ii] += rx * ffac;
                    cy[ii] += ry * ffac;
                    cz[ii] += rz * ffac;

                    cx[j] -= rx * ffac;
                    cy[j] -= ry * ffac;
                    cz[j] -= rz * ffac;
                }
            }
        }

        // OMP instruction - reduce forces over threads 
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        int chunk_size = (sys->natoms + sys->nthreads - 1) / sys->nthreads;
        int start = tid * chunk_size;
        int end = (start + chunk_size > sys->natoms) ? sys->natoms : start + chunk_size;
        for (int t = 1; t < sys->nthreads; ++t) {
            int offset = t * sys->natoms;
            for (int i = start; i < end; ++i) {
                sys->cx[i] += sys->cx[offset + i];
                sys->cy[i] += sys->cy[offset + i];
                sys->cz[i] += sys->cz[offset + i];
            }
        }
    } // OMP end parallel region 

    // OMP instruction - update global forces and epot
    for (int i = 0; i < sys->natoms; ++i) {
        sys->fx[i] = sys->cx[i];
        sys->fy[i] = sys->cy[i];
        sys->fz[i] = sys->cz[i];
    }
    sys->epot = epot_local;
    */
}



