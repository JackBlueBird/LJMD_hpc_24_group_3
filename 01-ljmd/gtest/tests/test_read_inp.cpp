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

class InputDataTest: public ::testing::Test {

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

TEST_F(InputDataTest, step1)
{
    ASSERT_NE(sys,nullptr);
    /* Assert input data */
    ASSERT_EQ(sys->natoms,3);
    ASSERT_DOUBLE_EQ(sys->mass,39.948);
    ASSERT_DOUBLE_EQ(sys->epsilon,0.2379);
    ASSERT_DOUBLE_EQ(sys->sigma,3.405);
    ASSERT_DOUBLE_EQ(sys->rcut,8.5);
    ASSERT_DOUBLE_EQ(sys->box,17.1580);
    ASSERT_STREQ(restfile,"argon_3.rest");
    ASSERT_STREQ(trajfile,"argon_3.xyz");
    ASSERT_STREQ(ergfile,"argon_3.dat");
    ASSERT_EQ(sys->nsteps,10);
    ASSERT_DOUBLE_EQ(sys->dt,5.0);
    ASSERT_EQ(nprint,1);

}
