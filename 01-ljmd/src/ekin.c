#include "ekin.h"

// compute kinetic energy 
/*
void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
*/

void ekin(mdsys_t* sys) {
	int i, ii;
	double local_ekin = 0.0;
	sys->ekin = 0.0;

	//MPI_Bcast(sys->vx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(sys->vy, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(sys->vz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(i, ii) reduction(+ : local_ekin)
#endif
	for (i = 0; i < sys->natoms; i += sys->nsize) {
		ii = i; //+ sys->mpirank;
		if (ii < (sys->natoms)) {
			local_ekin +=
				0.5 * mvsq2e * sys->mass *
				(sys->vx[ii] * sys->vx[ii] + sys->vy[ii] * sys->vy[ii] + sys->vz[ii] * sys->vz[ii]);
		}
	}

    sys->ekin = local_ekin;
    //MPI_Allreduce(&local_ekin, &(sys->ekin), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // remember to put MPI_Barrier in main after restar, before sys.nfi=0, 193 in main 
	sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}