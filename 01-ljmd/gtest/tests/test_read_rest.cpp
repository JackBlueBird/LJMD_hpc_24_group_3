// unit test example with test fixture
#include <iostream>
#include <fstream>
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

class InputRestTest: public ::testing::Test {

protected:

    mdsys_t *sys;
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;

    void SetUp()
    {
      sys = new mdsys_t;

      /* read input file */
      if(get_a_line(stdin,line)) return;
      sys->natoms=atoi(line);
      if(get_a_line(stdin,line)) return;
      sys->mass=atof(line);
      if(get_a_line(stdin,line)) return;
      sys->epsilon=atof(line);
      if(get_a_line(stdin,line)) return;
      sys->sigma=atof(line);
      if(get_a_line(stdin,line)) return;
      sys->rcut=atof(line);
      if(get_a_line(stdin,line)) return;
      sys->box=atof(line);
      if(get_a_line(stdin,restfile)) return;
      if(get_a_line(stdin,trajfile)) return;
      if(get_a_line(stdin,ergfile)) return;
      if(get_a_line(stdin,line)) return;
      sys->nsteps=atoi(line);
      if(get_a_line(stdin,line)) return;
      sys->dt=atof(line);
      if(get_a_line(stdin,line)) return;
      nprint=atoi(line);

      sys->rx = new double[sys->natoms];
      sys->ry = new double[sys->natoms];
      sys->rz = new double[sys->natoms];
      sys->vx = new double[sys->natoms];
      sys->vy = new double[sys->natoms];
      sys->vz = new double[sys->natoms];
    }

    void TearDown()
        {
          delete[] sys->rx;
          delete[] sys->ry;
          delete[] sys->rz;
          delete[] sys->vx;
          delete[] sys->vy;
          delete[] sys->vz;

          delete sys;
        }
};

TEST_F(InputRestTest, step2)
{
    ASSERT_NE(sys,nullptr);

    /* read restart */
    fp=fopen(restfile,"r");
    ASSERT_FALSE(fp == NULL);

    SUCCEED() << "Read .rest data file.";

    if(fp) {
        for (i=0; i<sys->natoms; ++i) {
          fscanf(fp,"%lf%lf%lf",sys->rx+i, sys->ry+i, sys->rz+i);
        }
        for (i=0; i<sys->natoms; ++i) {
          fscanf(fp,"%lf%lf%lf",sys->vx+i, sys->vy+i, sys->vz+i);
        }
        fclose(fp);
    } else {
        perror("cannot read restart file");
        return;
    }

    /* Assert initial status */
    ASSERT_DOUBLE_EQ(sys->rx[0],6.67);
    ASSERT_DOUBLE_EQ(sys->ry[0],-10.61);
    ASSERT_DOUBLE_EQ(sys->rz[0],12.63);
    ASSERT_DOUBLE_EQ(sys->vx[0],-1.56e-03);
    ASSERT_DOUBLE_EQ(sys->vy[0],4.84e-04);
    ASSERT_DOUBLE_EQ(sys->vz[0],-4.33e-04);

    ASSERT_DOUBLE_EQ(sys->rx[1],1.06);
    ASSERT_DOUBLE_EQ(sys->ry[1],-3.33);
    ASSERT_DOUBLE_EQ(sys->rz[1],-2.59);
    ASSERT_DOUBLE_EQ(sys->vx[1],4.16e-04);
    ASSERT_DOUBLE_EQ(sys->vy[1],2.28e-05);
    ASSERT_DOUBLE_EQ(sys->vz[1],-6.19e-04);

    ASSERT_DOUBLE_EQ(sys->rx[2],-1.78);
    ASSERT_DOUBLE_EQ(sys->ry[2],-16.52);
    ASSERT_DOUBLE_EQ(sys->rz[2],4.61);
    ASSERT_DOUBLE_EQ(sys->vx[2],-7.56e-04);
    ASSERT_DOUBLE_EQ(sys->vy[2],4.07e-04);
    ASSERT_DOUBLE_EQ(sys->vz[2],-4.65e-04);

}

