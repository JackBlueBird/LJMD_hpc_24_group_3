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

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz = 0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */

struct _mdsys {
  int natoms, nfi, nsteps;
  double dt, mass, epsilon, sigma, box, rcut;
  double ekin, epot, temp;
  double *rx, *ry, *rz;
  double *vx, *vy, *vz;
  double *fx, *fy, *fz;
  double *hx, *hy, *hz;
  int nthreads;
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

/* compute kinetic energy */
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

/* compute forces */
static void force(mdsys_t *sys) {
    double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    double epot_local = 0.0;

    // Arrys of support hx -> helper x
    azzero(sys->hx, sys->nthreads * sys->natoms);
    azzero(sys->hy, sys->nthreads * sys->natoms);
    azzero(sys->hz, sys->nthreads * sys->natoms);

    #ifdef _OPENMP
    #pragma omp parallel num_threads(sys->nthreads)
    #endif
    {
        int tid;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #else
        tid = 0;
        #endif

        double rx, ry, rz, rsq, ffac;
        double *hx = sys->hx + tid * sys->natoms;
        double *hy = sys->hy + tid * sys->natoms;
        double *hz = sys->hz + tid * sys->natoms;

        // divide the work between threads
        int grids = sys->nthreads;
        for (int i = 0; i < sys->natoms - 1; i += grids) {
            int ii = i + tid; // new index 
            if (ii >= sys->natoms - 1) break;

            for (int j = ii + 1; j < sys->natoms; ++j) {
                // distances bettween particles
                rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
                ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
                rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
                rsq = rx * rx + ry * ry + rz * rz;

                // If is inside of the cutoff then compute forces and energies
                if (rsq < rcsq) {
                    double rsqinv = 1.0 / rsq;
                    double r6 = rsqinv * rsqinv * rsqinv;
                    ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rsqinv;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    epot_local += r6 * (c12 * r6 - c6);

                    hx[ii] += rx * ffac;
                    hy[ii] += ry * ffac;
                    hz[ii] += rz * ffac;

                    hx[j] -= rx * ffac;
                    hy[j] -= ry * ffac;
                    hz[j] -= rz * ffac;
                }
            }
        }

        // sincronization and combine the local  forces into global forces 
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        int chunk_size = (sys->natoms + sys->nthreads - 1) / sys->nthreads;
        int start = tid * chunk_size;
        int end = (start + chunk_size > sys->natoms) ? sys->natoms : start + chunk_size;

        for (int t = 1; t < sys->nthreads; ++t) {
            int offset = t * sys->natoms;
            for (int i = start; i < end; ++i) {
                sys->hx[i] += sys->hx[offset + i];
                sys->hy[i] += sys->hy[offset + i];
                sys->hz[i] += sys->hz[offset + i];
            }
        }
    }

    // Update global forces and epot
    for (int i = 0; i < sys->natoms; ++i) {
        sys->fx[i] = sys->hx[i];
        sys->fy[i] = sys->hy[i];
        sys->fz[i] = sys->hz[i];
    }

    sys->epot = epot_local;
}


/* velocity verlet */
static void velverlet(mdsys_t *sys) {
  int i;

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
  int nprint, i;
  char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
  FILE *fp, *traj, *erg;
  mdsys_t sys;
  double t_start;

  /*Initialice ntheads*/
  #if defined(_OPENMP)
    #pragma omp parallel
    sys.nthreads = omp_get_num_threads();
  #else
    sys.nthreads = 1;
  #endif
  
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
	sys.hx = (double*)malloc(sys.nthreads * sys.natoms * sizeof(double));
	sys.hy = (double*)malloc(sys.nthreads * sys.natoms * sizeof(double));
	sys.hz = (double*)malloc(sys.nthreads * sys.natoms * sizeof(double));


  /* read restart */
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

  /* initialize forces and energies.*/
  sys.nfi = 0;
  force(&sys);
  ekin(&sys);

  erg = fopen(ergfile, "w");
  traj = fopen(trajfile, "w");

  printf("Startup time: %10.3fs\n", wallclock() - t_start);
  printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms,
         sys.nsteps);
  printf("     NFI            TEMP            EKIN                 EPOT        "
         "      ETOT\n");
  output(&sys, erg, traj);

  /* reset timer */
  t_start = wallclock();

  /**************************************************/
  /* main MD loop */
  for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {

    /* write output, if requested */
    if ((sys.nfi % nprint) == 0)
      output(&sys, erg, traj);

    /* propagate system and recompute energies */
    velverlet(&sys);
    ekin(&sys);
  }
  /**************************************************/

  /* clean up: close files, free memory */
  printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
  fclose(erg);
  fclose(traj);

  free(sys.rx);
  free(sys.ry);
  free(sys.rz);
  free(sys.vx);
  free(sys.vy);
  free(sys.vz);
  free(sys.fx);
  free(sys.fy);
  free(sys.fz);
  free(sys.hx);
  free(sys.hy);
  free(sys.hz);

  return 0;
}