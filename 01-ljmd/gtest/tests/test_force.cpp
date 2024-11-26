// unit test example with test fixture
#include "gtest/gtest.h"

#include "mdsys.h"
#include "read_input.h"
#include "utilities.h"
#include "ekin.h"
#include "force.h"
#include "verlet.h"
#include "output.h"

const double abs_err=1e-6;

class ForceTest: public ::testing::Test {

  protected:

    mdsys_t *sys;

    void SetUp()
    {
      sys = new mdsys_t;
      sys->natoms = 4;
      sys->mass = 39.948;
      sys->epsilon = 0.2379;
      sys->sigma = 3.405;
      sys->rcut = 2.6;
      sys->box = 8.0;
      sys->dt = 5.0;

      sys->rx = new double[4];
      sys->fx = new double[4];

      sys->ry = new double[4];
      sys->fy = new double[4];

      sys->rz = new double[4];
      sys->fz = new double[4];

      // Zero energy and forces
      sys->epot = 0.0;

      sys->fx[0] = 0.0;
      sys->fx[1] = 0.0;
      sys->fx[2] = 0.0;
      sys->fx[3] = 0.0;
      sys->fy[0] = 0.0;
      sys->fy[1] = 0.0;
      sys->fy[2] = 0.0;
      sys->fy[3] = 0.0;
      sys->fz[0] = 0.0;
      sys->fz[1] = 0.0;
      sys->fz[2] = 0.0;
      sys->fz[3] = 0.0;

      // Position
      sys->rx[0] = 0.5; // |2     |
      sys->rx[1] = 2.0; // |    3 |
      sys->rx[2] = 0.5; // | 1    |
      sys->rx[3] = 4.0; // |0     |

      sys->ry[0] = 0.5;
      sys->ry[1] = 2.0;
      sys->ry[2] = 0.5;
      sys->ry[3] = 4.0;

      sys->rz[0] = 0.5;
      sys->rz[1] = 2.0;
      sys->rz[2] = 7.5;
      sys->rz[3] = 4.0;
    }

    void TearDown()
    {
      delete[] sys->rx;
      delete[] sys->fx;

      delete[] sys->ry;
      delete[] sys->fy;

      delete[] sys->rz;
      delete[] sys->fz;

      delete sys;
    }
};

// Unit test - force - inside cutoff - directly (p0,p1)
TEST_F(ForceTest, inside_directly)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->fx[0],0.0,abs_err);
  ASSERT_NEAR(sys->fy[0],0.0,abs_err);
  ASSERT_NEAR(sys->fz[0],0.0,abs_err);
  ASSERT_NEAR(sys->fx[1],0.0,abs_err);
  ASSERT_NEAR(sys->fy[1],0.0,abs_err);
  ASSERT_NEAR(sys->fz[1],0.0,abs_err);
  force(sys);
  ASSERT_NEAR(sys->fx[1],58.73411942344227,abs_err);
  ASSERT_NEAR(sys->fy[1],58.73411942344227,abs_err);
  ASSERT_NEAR(sys->fz[1],58.73411942344227,abs_err);
}

// Unit test - force - inside cutoff - via PBC (p0,p2)
TEST_F(ForceTest, inside_viaPBC)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->fx[0],0.0,abs_err);
  ASSERT_NEAR(sys->fy[0],0.0,abs_err);
  ASSERT_NEAR(sys->fz[0],0.0,abs_err);
  ASSERT_NEAR(sys->fx[2],0.0,abs_err);
  ASSERT_NEAR(sys->fy[2],0.0,abs_err);
  ASSERT_NEAR(sys->fz[2],0.0,abs_err);
  force(sys);
  ASSERT_NEAR(sys->fx[2],0.0,abs_err);
  ASSERT_NEAR(sys->fy[2],0.0,abs_err);
  ASSERT_NEAR(sys->fz[2],-27726925.7743613,abs_err); // negative by Newton 3
}

// Unit test - force - outside cutoff (p0,p3)
TEST_F(ForceTest, outside_directly)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->fx[0],0.0,abs_err);
  ASSERT_NEAR(sys->fy[0],0.0,abs_err);
  ASSERT_NEAR(sys->fz[0],0.0,abs_err);
  ASSERT_NEAR(sys->fx[3],0.0,abs_err);
  ASSERT_NEAR(sys->fy[3],0.0,abs_err);
  ASSERT_NEAR(sys->fz[3],0.0,abs_err);
  force(sys);
  ASSERT_NEAR(sys->fx[3],0.0,abs_err);
  ASSERT_NEAR(sys->fy[3],0.0,abs_err);
  ASSERT_NEAR(sys->fz[3],0.0,abs_err);
}
