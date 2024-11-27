/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>

/* Include version number */
#include "LJMDConfig.h"

/* #define BLEN 200 */
#include "mdsys.h"
#include "read_input.h"
#include "utilities.h"
#include "ekin.h"
#include "force.h"
#include "verlet.h"
#include "output.h"


/* main */
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0)
        {printf("rank 0 is stepping here, file reading \n");}
        
    if (rank == 1)
        {printf("rank 1 is stepping here, file reading \n");}

    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    /*Initialice ntheads*/
    #if defined(_OPENMP)
    #pragma omp parallel
        sys.nthreads = omp_get_num_threads();
    #else
        sys.nthreads = 1;
    #endif

    int nthreads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    #endif

    printf("LJMD version %3.1f\n", LJMD_VERSION);
    printf("Number of threads equal to %d \n",nthreads);
    t_start = wallclock();

    /* read input file */
    if (rank == 0) {
        if(get_a_line(stdin,line)) return 1;
        sys.natoms=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.mass=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.epsilon=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.sigma=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.rcut=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.box=atof(line);
        if(get_a_line(stdin,restfile)) return 1;
        if(get_a_line(stdin,trajfile)) return 1;
        if(get_a_line(stdin,ergfile)) return 1;
        if(get_a_line(stdin,line)) return 1;
        sys.nsteps=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.dt=atof(line);
        if(get_a_line(stdin,line)) return 1;
        nprint=atoi(line);
    }

    // MPI - Broadcast to other process
    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nprint, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    // OMP instruction - allocate support array for forces
	sys.cx = (double*)malloc(nthreads * sys.natoms * sizeof(double));
	sys.cy = (double*)malloc(nthreads * sys.natoms * sizeof(double));
	sys.cz = (double*)malloc(nthreads * sys.natoms * sizeof(double));

    /* read restart */
    if (rank == 0) {
        fp=fopen(restfile,"r");
        if(fp) 
        {
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
            }
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
            }
            fclose(fp);
            azzero(sys.fx, sys.natoms);
            azzero(sys.fy, sys.natoms);
            azzero(sys.fz, sys.natoms);
        } else {
            perror("cannot read restart file");
            return 3;
        }
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);

    printf("Startup time, for rank %d. : %10.3fs\n",rank ,wallclock()-t_start);
        if (rank == 0) {
            erg=fopen(ergfile,"w");
            traj=fopen(trajfile,"w");        
            printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
            printf("     NFI            TEMP       instruction     EKIN                 EPOT              ETOT\n");
            output(&sys, erg, traj);
        }    

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if (rank == 0) {
            if ((sys.nfi % nprint) == 0) {
                output(&sys, erg, traj);
            }
        }
        velverlet1(&sys);
        force(&sys);
        velverlet2(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done for rank %d. Run time: %10.3fs\n",rank,wallclock()-t_start);
    if (rank == 0) {    
        fclose(erg);
        fclose(traj);
    }

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);

    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    // OMP instruction - de allocate auxiliary storages
    free(sys.cx);
    free(sys.cy);
    free(sys.cz);

    MPI_Finalize();
    return 0;
}