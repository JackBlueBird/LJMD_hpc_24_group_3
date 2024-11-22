#ifndef MDSYS_H
#define MDSYS_H

typedef struct {
    int natoms, nsteps, nfi;
    double dt, mass, epsilon, sigma, box, rcut;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    double epot, ekin, temp;
} mdsys_t;

extern const double kboltz;
extern const double mvsq2e;

#endif
