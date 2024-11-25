#include <math.h>
#include "mdsys.h"
#include "utilities.h"
#include "force.h"

/* compute forces */
void force(mdsys_t *sys)
{
    double r,ffac;
    double rsq;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
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
}



