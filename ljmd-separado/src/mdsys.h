#ifndef BLEN 
#define BLEN 200
#endif 

#ifndef PHYS_CONST 
#define PHYS_CONST
extern const double kboltz;
extern const double mvsq2e;
#endif


#ifndef MDSYS_H
#define MDSYS_H

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

#endif
