#include <stdio.h>

#include "gtest/gtest.h"

#include "mdsys.h"
#include "read_input.h"
#include "utilities.h"
#include "ekin.h"
#include "force.h"
#include "verlet.h"
#include "output.h"

/* generic file- or pathname buffer length */
#define BLEN 200

TEST(ReadInput, params)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;

    /* read input file */
    if(get_a_line(stdin,line)) return;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return;
    sys.mass=atof(line);
    if(get_a_line(stdin,line)) return;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return;
    if(get_a_line(stdin,trajfile)) return;
    if(get_a_line(stdin,ergfile)) return;
    if(get_a_line(stdin,line)) return;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return;
    nprint=atoi(line);

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

    /* Assert input data */
    ASSERT_EQ(sys.natoms,3);
    ASSERT_DOUBLE_EQ(sys.mass,39.948);
    ASSERT_DOUBLE_EQ(sys.epsilon,0.2379);
    ASSERT_DOUBLE_EQ(sys.sigma,3.405);
    ASSERT_DOUBLE_EQ(sys.rcut,8.5);
    ASSERT_DOUBLE_EQ(sys.box,17.1580);
    ASSERT_STREQ(restfile,"argon_3.rest");
    ASSERT_STREQ(trajfile,"argon_3.xyz");
    ASSERT_STREQ(ergfile,"argon_3.dat");
    ASSERT_EQ(sys.nsteps,10);
    ASSERT_DOUBLE_EQ(sys.dt,5.0);
    ASSERT_EQ(nprint,1);

    /* read restart */
    fp=fopen(restfile,"r");

    ASSERT_FALSE(fp == NULL);

    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
    } else {
        perror("cannot read restart file");
        return;
    }

    /* Assert initial status */
    ASSERT_DOUBLE_EQ(sys.rx[0],6.67);
    ASSERT_DOUBLE_EQ(sys.ry[0],-10.61);
    ASSERT_DOUBLE_EQ(sys.rz[0],12.63);
    ASSERT_DOUBLE_EQ(sys.vx[0],-1.56e-03);
    ASSERT_DOUBLE_EQ(sys.vy[0],4.84e-04);
    ASSERT_DOUBLE_EQ(sys.vz[0],-4.33e-04);

    ASSERT_DOUBLE_EQ(sys.rx[1],1.06);
    ASSERT_DOUBLE_EQ(sys.ry[1],-3.33);
    ASSERT_DOUBLE_EQ(sys.rz[1],-2.59);
    ASSERT_DOUBLE_EQ(sys.vx[1],4.16e-04);
    ASSERT_DOUBLE_EQ(sys.vy[1],2.28e-05);
    ASSERT_DOUBLE_EQ(sys.vz[1],-6.19e-04);

    ASSERT_DOUBLE_EQ(sys.rx[2],-1.78);
    ASSERT_DOUBLE_EQ(sys.ry[2],-16.52);
    ASSERT_DOUBLE_EQ(sys.rz[2],4.61);
    ASSERT_DOUBLE_EQ(sys.vx[2],-7.56e-04);
    ASSERT_DOUBLE_EQ(sys.vy[2],4.07e-04);
    ASSERT_DOUBLE_EQ(sys.vz[2],-4.65e-04);

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
}
