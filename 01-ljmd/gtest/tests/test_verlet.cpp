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

class VerletTest: public ::testing::Test {

  protected:

    mdsys_t *sys;

    void SetUp()
    {
      sys = new mdsys_t;
      sys->natoms = 3;
      sys->mass = 39.948;
      sys->dt = 5.0;

      sys->rx = new double[3];
      sys->vx = new double[3];
      sys->fx = new double[3];

      sys->ry = new double[3];
      sys->vy = new double[3];
      sys->fy = new double[3];

      sys->rz = new double[3];
      sys->vz = new double[3];
      sys->fz = new double[3];

      sys->rx[0] = -1.0;
      sys->rx[1] = 1.0;
      sys->rx[2] = 0.0;
      sys->vx[0] = 0.0;
      sys->vx[1] = 0.0;
      sys->vx[2] = 0.0;
      sys->fx[0] = 1.0;
      sys->fx[1] = 0.2;
      sys->fx[2] = 0.0;

      sys->ry[0] = -1.0;
      sys->ry[1] = 1.0;
      sys->ry[2] = -3.5;
      sys->vy[0] = 0.0;
      sys->vy[1] = 0.0;
      sys->vy[2] = -0.5;
      sys->fy[0] = 1.0;
      sys->fy[1] = 0.2;
      sys->fy[2] = 0.5;

      sys->rz[0] = -1.0;
      sys->rz[1] = 8.0;
      sys->rz[2] = 0.0;
      sys->vz[0] = 0.0;
      sys->vz[1] = 2.5;
      sys->vz[2] = 0.0;
      sys->fz[0] = 1.0;
      sys->fz[1] = -10.0;
      sys->fz[2] = 0.0;
    }

    void TearDown()
    {
      delete[] sys->rx;
      delete[] sys->vx;
      delete[] sys->fx;

      delete[] sys->ry;
      delete[] sys->vy;
      delete[] sys->fy;

      delete[] sys->rz;
      delete[] sys->vz;
      delete[] sys->fz;

      delete sys;
    }
};

// Unit test - time integration - x component
TEST_F(VerletTest, step1_x)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->rx[0],-1.0,abs_err);
  ASSERT_NEAR(sys->vx[0],0.0,abs_err);
  velverlet1(sys);
  ASSERT_NEAR(sys->rx[0],-0.9998690798037535,abs_err);
  ASSERT_NEAR(sys->vx[0],0.00002618403924930834,abs_err);
}

TEST_F(VerletTest, step2_x)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->rx[0],-1.0,abs_err);
  ASSERT_NEAR(sys->vx[0],0.0,abs_err);
  velverlet2(sys);
  ASSERT_NEAR(sys->rx[0],-1.0,abs_err);
  ASSERT_NEAR(sys->vx[0],0.00002618403924930834,abs_err);
}

// Unit test - time integration - y component
TEST_F(VerletTest, step1_y)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->ry[2],-3.5,abs_err);
  ASSERT_NEAR(sys->vy[2],-0.5,abs_err);
  velverlet1(sys);
  ASSERT_NEAR(sys->ry[2],-5.999934539901877,abs_err);
  ASSERT_NEAR(sys->vy[2],-0.4999869079803754,abs_err);
}

TEST_F(VerletTest, step2_y)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->ry[2],-3.5,abs_err);
  ASSERT_NEAR(sys->vy[2],-0.5,abs_err);
  velverlet2(sys);
  ASSERT_NEAR(sys->ry[2],-3.5,abs_err);
  ASSERT_NEAR(sys->vy[2],-0.4999869079803754,abs_err);
}

// Unit test - time integration - z component
TEST_F(VerletTest, step1_z)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->rz[1],8.0,abs_err);
  ASSERT_NEAR(sys->vz[1],2.5,abs_err);
  velverlet1(sys);
  ASSERT_NEAR(sys->rz[1],20.49869079803754,abs_err);
  ASSERT_NEAR(sys->vz[1],2.499738159607507,abs_err);
}

TEST_F(VerletTest, step2_z)
{
  ASSERT_NE(sys,nullptr);
  ASSERT_NEAR(sys->rz[1],8.0,abs_err);
  ASSERT_NEAR(sys->vz[1],2.5,abs_err);
  velverlet2(sys);
  ASSERT_NEAR(sys->rz[1],8.0,abs_err);
  ASSERT_NEAR(sys->vz[1],2.499738159607507,abs_err);
}
