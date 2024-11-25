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

const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

class EkinTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 3;
        sys->mass = 1.0;
        sys->ekin = 0.0;
        // sys->temp = 0.0;
        sys->vx = new double[3];
        sys->vy = new double[3];
        sys->vz = new double[3];
        sys->vx[0] = 0.300;
        sys->vx[1] = 0.150;
        sys->vx[2] = 0.070;
        sys->vy[0] = 0.200;
        sys->vy[1] = 0.100;
        sys->vy[2] = 0.050;
        sys->vz[0] = 0.100;
        sys->vz[1] = 0.050;
        sys->vz[2] = 0.025;
    }

    void TearDown()
        {
            delete[] sys->vx;
            delete[] sys->vy;
            delete[] sys->vz;

            delete sys;
        }
};

TEST_F(EkinTest, step1)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_NEAR(sys->ekin,0.0,abs_err);
    ekin(sys);
    ASSERT_NEAR(sys->ekin,218.7201242973335,abs_err);
}

TEST_F(EkinTest, step2)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_NEAR(sys->ekin,0.0,abs_err);
    ekin(sys);
    ASSERT_NEAR(sys->temp,36688.03456586131,abs_err);
}
