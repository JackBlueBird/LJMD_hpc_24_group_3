/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few pcysical constants */
const double kboltz = 0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */

struct _mdsys {
  MPI_Comm mpicomm;
  int natoms, nfi, nsteps;
  double dt, mass, epsilon, sigma, box, rcut;
  double ekin, epot, temp;
  double *rx, *ry, *rz;
  double *vx, *vy, *vz;
  double *fx, *fy, *fz;
  double *cx, *cy, *cz;
  int nthreads, nsize, mpirank;

};
typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE *fp, char *buf) {
  char tmp[BLEN], *ptr;

  /* read a line and cut of comments and blanks */
  if (fgets(tmp, BLEN, fp)) {
    int i;

    ptr = strchr(tmp, '#');
    if (ptr)
      *ptr = '\0';
    i = strlen(tmp);
    --i;
    while (isspace(tmp[i])) {
      tmp[i] = '\0';
      --i;
    }
    ptr = tmp;
    while (isspace(*ptr)) {
      ++ptr;
    }
    i = strlen(ptr);
    strcpy(buf, tmp);
    return 0;
  } else {
    perror("problem reading input");
    return -1;
  }
  return 0;
}

/* helper function: get current time in seconds since epoch */


static double wallclock() {
  struct timeval t;
  gettimeofday(&t, 0);
  return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

/* helper function: zero out an array */
// gz: setta un array a 0
static void azzero(double *d, const int n) {
  int i;
  for (i = 0; i < n; ++i) {
    d[i] = 0.0;
  }
}

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2) {
  while (x > boxby2)
    x -= 2.0 * boxby2;
  while (x < -boxby2)
    x += 2.0 * boxby2;
  return x;
}

// compute kinetic energy 
/*
static void ekin(mdsys_t *sys) {
  int i;

  sys->ekin = 0.0;
  for (i = 0; i < sys->natoms; ++i) {
    sys->ekin += 0.5 * mvsq2e * sys->mass *
                 (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] +
                  sys->vz[i] * sys->vz[i]);
  }
  sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}
*/
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

/* compute forces */
static void force(mdsys_t *sys) {

    double sigma, sigma6;
    sigma=sys->sigma;
    sigma6=sigma*sigma*sigma*sigma*sigma*sigma;

    double c12 = 4.0 * sys->epsilon * sigma6*sigma6;//pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * sigma6;//pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    double epot = 0.0;
    int tid;

    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	  MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	  MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

    // Arrys of support cx -> helper x
    azzero(sys->cx, sys->nthreads * sys->natoms);
    azzero(sys->cy, sys->nthreads * sys->natoms);
    azzero(sys->cz, sys->nthreads * sys->natoms);

    #ifdef _OPENMP
    #pragma omp parallel private(tid)
    #endif
    {
        int tid=0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        //#else
        //tid = 0;
        #endif
        //double nsize = 1;// this will become sys->nsize
        double rx, ry, rz, rsq, ffac;
        double *cx = sys->cx + tid * sys->natoms;
        double *cy = sys->cy + tid * sys->natoms;
        double *cz = sys->cz + tid * sys->natoms;

        // divide the work between threads
        int grids = sys->nthreads * sys->nsize;// nsize is the number of ranks
        //double mpirank =0;//sys->mpirank;
        
        for (int i = 0; i < sys->natoms - 1; i += grids) {
            int ii = i + tid + sys->nthreads * sys->mpirank; // mpirank will become sys->mpirank
            if (ii >= sys->natoms - 1) break;
                double rx1 = sys->rx[ii];
                double ry1 = sys->ry[ii];
                double rz1 = sys->rz[ii];
            for (int j = ii + 1; j < (sys->natoms); ++j) {
                // distances bettween particles
                rx = pbc(rx1 - sys->rx[j], 0.5 * sys->box);
                ry = pbc(ry1 - sys->ry[j], 0.5 * sys->box);
                rz = pbc(rz1 - sys->rz[j], 0.5 * sys->box);
                rsq = rx * rx + ry * ry + rz * rz;

                // If is inside of the cutoff then compute forces and energies
                if (rsq < rcsq) {
                    double rsqinv = 1.0 / rsq;
                    double r6 = rsqinv * rsqinv * rsqinv;
                    ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rsqinv;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    epot += r6 * (c12 * r6 - c6);

                    cx[ii] += rx * ffac;
                    cy[ii] += ry * ffac;
                    cz[ii] += rz * ffac;

                    cx[j] -= rx * ffac;
                    cy[j] -= ry * ffac;
                    cz[j] -= rz * ffac;
                }
            }
        }

    // sincronization and combine the local  forces into global forces 
    #ifdef _OPENMP
	  #pragma omp barrier	 // sync everything
    #endif
		int i = 1 + (sys->natoms / sys->nthreads);
		int fromidx = tid * i;
		int toidx = fromidx + i;
		if (toidx > sys->natoms)
			toidx = sys->natoms;
		for (i = 1; i < sys->nthreads; ++i) {
			int offs = i * sys->natoms;
			for (int j = fromidx; j < toidx; ++j) {
				sys->cx[j] += sys->cx[offs + j];
				sys->cy[j] += sys->cy[offs + j];
				sys->cz[j] += sys->cz[offs + j];
			}
		}
	}

	MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}


/* velocity verlet */
static void velverlet(mdsys_t *sys) {
  int i;
  double inv_mass = 1 / sys->mass;
	double inv_mvsq2e = 1 / mvsq2e;

  /* first part: propagate velocities by half and positions by full step */
  // vx = 0.5 * dt * f / ( mass * other constant)
  for (i = 0; i < sys->natoms; ++i) {
    sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
    sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
    sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    sys->rx[i] += sys->dt * sys->vx[i];
    sys->ry[i] += sys->dt * sys->vy[i];
    sys->rz[i] += sys->dt * sys->vz[i];
  }

  /* compute forces and potential energy */
  force(sys);

  /* second part: propagate velocities by another half step */
  for (i = 0; i < sys->natoms; ++i) {
    sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
    sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
    sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
  }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj) {
  int i;

  printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
         sys->ekin, sys->epot, sys->ekin + sys->epot);
  fprintf(erg, "% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
          sys->ekin, sys->epot, sys->ekin + sys->epot);
  fprintf(traj, "%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi,
          sys->ekin + sys->epot);
  for (i = 0; i < sys->natoms; ++i) {
    fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i],
            sys->rz[i]);
  }
}

/* main */
int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  int nprint, i;
  char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
  FILE *fp, *traj, *erg;
  mdsys_t sys;
  double t_start;

  // MPI parameters
  sys.mpirank = rank;
  sys.nsize = size;
  sys.mpicomm = MPI_COMM_WORLD;

  /*Initialice ntheads*/
  #if defined(_OPENMP)
    #pragma omp parallel
    sys.nthreads = omp_get_num_threads();
  #else
    sys.nthreads = 1;
  #endif

  if(sys.mpirank==0){
  printf("num theads= %d\n",sys.nthreads);
  printf("LJMD version %3.1f\n", LJMD_VERSION);

  t_start = wallclock();
  
  /* read input file */
  if (get_a_line(stdin, line))
    return 1;
  sys.natoms = atoi(line);
  if (get_a_line(stdin, line))
    return 1;
  sys.mass = atof(line);
  if (get_a_line(stdin, line))
    return 1;
  sys.epsilon = atof(line);
  if (get_a_line(stdin, line))
    return 1;
  sys.sigma = atof(line);
  if (get_a_line(stdin, line))
    return 1;
  sys.rcut = atof(line);
  if (get_a_line(stdin, line))
    return 1;
  sys.box = atof(line);
  if (get_a_line(stdin, restfile))
    return 1;
  if (get_a_line(stdin, trajfile))
    return 1;
  if (get_a_line(stdin, ergfile))
    return 1;
  if (get_a_line(stdin, line))
    return 1;
  sys.nsteps = atoi(line);
  if (get_a_line(stdin, line))
    return 1;
  sys.dt = atof(line);
  if (get_a_line(stdin, line))
    return 1;
  nprint = atoi(line);
  }// end reading data 

  // Send data to all process 
  /*Sending data to all MPI processors*/
	MPI_Bcast(&(sys.natoms), 1, MPI_INT, 0, sys.mpicomm);
	MPI_Bcast(&(sys.nsteps), 1, MPI_INT, 0, sys.mpicomm);
	MPI_Bcast(&(sys.mass), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.epsilon), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.sigma), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.box), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.rcut), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.dt), 1, MPI_DOUBLE, 0, sys.mpicomm);

  if (sys.mpirank == 0) {
		printf("Communication time: %10.3fs\n", wallclock() - t_start);
	}

  /* allocate memory */
  sys.rx = (double *)malloc(sys.natoms * sizeof(double));
  sys.ry = (double *)malloc(sys.natoms * sizeof(double));
  sys.rz = (double *)malloc(sys.natoms * sizeof(double));
  sys.vx = (double *)malloc(sys.natoms * sizeof(double));
  sys.vy = (double *)malloc(sys.natoms * sizeof(double));
  sys.vz = (double *)malloc(sys.natoms * sizeof(double));
  sys.fx = (double *)malloc(sys.natoms * sizeof(double));
  sys.fy = (double *)malloc(sys.natoms * sizeof(double));
  sys.fz = (double *)malloc(sys.natoms * sizeof(double));
  // allocate support array for forces
	sys.cx = (double*)malloc(sys.nthreads * sys.natoms * sizeof(double));
	sys.cy = (double*)malloc(sys.nthreads * sys.natoms * sizeof(double));
	sys.cz = (double*)malloc(sys.nthreads * sys.natoms * sizeof(double));


  /* read restart */
  if(sys.mpirank==0){
  fp = fopen(restfile, "r");
  if (fp) {
    for (i = 0; i < sys.natoms; ++i) {
      fscanf(fp, "%lf%lf%lf", sys.rx + i, sys.ry + i, sys.rz + i);
    }
    for (i = 0; i < sys.natoms; ++i) {
      fscanf(fp, "%lf%lf%lf", sys.vx + i, sys.vy + i, sys.vz + i);
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
  MPI_Barrier(sys.mpicomm);
  sys.nfi = 0;

  force(&sys);
  ekin(&sys);

if(sys.mpirank==0){
  erg = fopen(ergfile, "w");
  traj = fopen(trajfile, "w");

  printf("Startup time: %10.3fs\n", wallclock() - t_start);
  printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms,
         sys.nsteps);
  printf("     NFI            TEMP            EKIN                 EPOT        "
         "      ETOT\n");
  output(&sys, erg, traj);
}
  /* reset timer */
  t_start = wallclock();

  /**************************************************/
  /* main MD loop */
  for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {

    /* write output, if requested */
    if(sys.mpirank==0){
    if ((sys.nfi % nprint) == 0)
      output(&sys, erg, traj);
    }
    /* propagate system and recompute energies */
    velverlet(&sys);
    ekin(&sys);
  }
  /**************************************************/

  /* clean up: close files, free memory */
  if(sys.mpirank==0){
  printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
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
  free(sys.cx);
  free(sys.cy);
  free(sys.cz);

  MPI_Finalize();
  return 0;
}