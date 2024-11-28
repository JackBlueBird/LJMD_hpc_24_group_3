// unit test example with test fixture
#include "gtest/gtest.h"

#include "mdsys.h"
#include "read_input.h"
#include "utilities.h"
#include "ekin.h"
#include "force.h"
#include "verlet.h"
#include "output.h"

#if defined(_MPI)
#include <mpi.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

const double abs_err=1e-6;

class ForceTest: public ::testing::Test {

  protected:

    mdsys_t *sys;

    void SetUp()
    {
#if defined(_MPI)
      char** argv;
      int argc = 0;
      int mpiError = MPI_Init(&argc, &argv);
      ASSERT_FALSE(mpiError);
#endif

      sys = new mdsys_t;
      sys->natoms = 4;
      sys->mass = 39.948;
      sys->epsilon = 0.2379;
      sys->sigma = 3.405;
      sys->rcut = 2.6;
      sys->box = 8.0;
      sys->dt = 5.0;

      // MPI - parameters
# if defined(_MPI)
      int rank,size;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      sys->mpirank = rank;
      sys->nsize = size;
# else
      sys->mpirank = 0;
      sys->nsize = 1;
# endif

#if defined(_OPENMP)
#pragma omp parallel
      {
      sys->nthreads = omp_get_num_threads();
      }
#else
      sys->nthreads = 1;
#endif

// #if defined(_MPI) || defined(_OPENMP)
      sys->cx = new double[sys->nthreads*4];
      sys->cy = new double[sys->nthreads*4];
      sys->cz = new double[sys->nthreads*4];
// #endif

      sys->rx = new double[sys->nthreads*4];
      sys->fx = new double[sys->nthreads*4];

      sys->ry = new double[sys->nthreads*4];
      sys->fy = new double[sys->nthreads*4];

      sys->rz = new double[sys->nthreads*4];
      sys->fz = new double[sys->nthreads*4];

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

#if defined(_MPI) || defined(_OPENMP)
      delete[] sys->cx;
      delete[] sys->cy;
      delete[] sys->cz;
#endif
      delete sys;

#if defined(_MPI)
      int mpiError = MPI_Finalize();
      ASSERT_FALSE(mpiError);
#endif
    }
};

TEST_F(ForceTest, interaction)
{
  ASSERT_NE(sys,nullptr);

  ASSERT_NEAR(sys->fx[0],0.0,abs_err);
  ASSERT_NEAR(sys->fy[0],0.0,abs_err);
  ASSERT_NEAR(sys->fz[0],0.0,abs_err);
  ASSERT_NEAR(sys->fx[1],0.0,abs_err);
  ASSERT_NEAR(sys->fy[1],0.0,abs_err);
  ASSERT_NEAR(sys->fz[1],0.0,abs_err);
  ASSERT_NEAR(sys->fx[2],0.0,abs_err);
  ASSERT_NEAR(sys->fy[2],0.0,abs_err);
  ASSERT_NEAR(sys->fz[2],0.0,abs_err);
  ASSERT_NEAR(sys->fx[3],0.0,abs_err);
  ASSERT_NEAR(sys->fy[3],0.0,abs_err);
  ASSERT_NEAR(sys->fz[3],0.0,abs_err);

  force(sys);

  // Unit test - force - inside cutoff - directly (p0,p1)
  ASSERT_NEAR(sys->fx[1],58.73411942344227,abs_err);
  ASSERT_NEAR(sys->fy[1],58.73411942344227,abs_err);
  ASSERT_NEAR(sys->fz[1],58.73411942344227,abs_err);
  // Unit test - force - inside cutoff - via PBC (p0,p2)
  ASSERT_NEAR(sys->fx[2],0.0,abs_err);
  ASSERT_NEAR(sys->fy[2],0.0,abs_err);
  ASSERT_NEAR(sys->fz[2],-27726925.7743613,abs_err); // negative by Newton 3

  // Unit test - force - outside cutoff (p0,p3)
  ASSERT_NEAR(sys->fx[3],0.0,abs_err);
  ASSERT_NEAR(sys->fy[3],0.0,abs_err);
  ASSERT_NEAR(sys->fz[3],0.0,abs_err);
}
