/* generic file- or pathname buffer length */


#ifndef BLEN
#define BLEN 200
#endif

/* a few physical constants */
#ifndef PHYS_CONSTS
#define PHYS_CONSTS
extern const double kboltz;     /* boltzman constant in kcal/mol/K */
extern const double mvsq2e;     /* m*v^2 in kcal/mol */
#endif

/* structure to hold the complete information
 * about the MD system */
#ifndef MDSYS_H
#define MDSYS_H
#ifdef __cplusplus
extern "C"
{
#endif

struct _mdsys {
    //  WATCH OUT FOR THE COMMUNICATOR
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    double *cx, *cy, *cz;
    int natoms, nfi, nsteps;
    int nthreads, nsize, mpirank;
};
/*
struct _mdsys {
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    int natoms,nfi,nsteps;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    # if defined(_MPI)
        double *cx, *cy, *cz;
    #endif
};
*/
typedef struct _mdsys mdsys_t;
#ifdef __cplusplus
}
#endif
#endif
