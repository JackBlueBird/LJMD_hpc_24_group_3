// compute kinetic energy

void ekin(mdsys_t* sys) {
	int i, ii;
	double partial_ekin = 0.0;
	sys->ekin = 0.0;

	MPI_Bcast(sys->vx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	MPI_Bcast(sys->vy, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	MPI_Bcast(sys->vz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(i, ii) reduction(+ : partial_ekin)
#endif
	for (i = 0; i < sys->natoms; i += sys->nsize) {
		ii = i + sys->mpirank;
		if (ii < (sys->natoms)) {
			partial_ekin +=
				0.5 * mvsq2e * sys->mass *
				(sys->vx[ii] * sys->vx[ii] + sys->vy[ii] * sys->vy[ii] + sys->vz[ii] * sys->vz[ii]);
		}
	}

	MPI_Allreduce(&partial_ekin, &(sys->ekin), 1, MPI_DOUBLE, MPI_SUM, sys->mpicomm);
	sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}