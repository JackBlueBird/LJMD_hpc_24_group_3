
#include "gtest/gtest.h"

#include "mdsys.h"
#include "read_input.h"
#include "utilities.h"
#include "ekin.h"
#include "force.h"
#include "verlet.h"
#include "output.h"


// TEST(TestWallClock, one)
// {
//     double t_test=-wallclock();
//     doublesleep(0.1);
//     t_test += wallclock();
//     ASSERT_TRUE((t_test>=0.1));
//     ASSERT_TRUE((t_test<1.1));
// }

TEST(TestAzzero, doubles)
{
    double *buf = new double[10];
    for (int i=0; i<10; ++i) buf[i]=i+1;
    ASSERT_DOUBLE_EQ(buf[1],2.0);
    ASSERT_DOUBLE_EQ(buf[5],6.0);
    ASSERT_DOUBLE_EQ(buf[9],10.0);

    azzero(buf,10);
    ASSERT_DOUBLE_EQ(buf[1],0.0);
    ASSERT_DOUBLE_EQ(buf[5],0.0);
    ASSERT_DOUBLE_EQ(buf[9],0.0);
}

TEST(TestAzzero, half)
{
    double *buf = new double[10];
    for (int i=0; i<10; ++i) buf[i]=i+1;
    ASSERT_DOUBLE_EQ(buf[1],2.0);
    ASSERT_DOUBLE_EQ(buf[5],6.0);
    ASSERT_DOUBLE_EQ(buf[9],10.0);

    azzero(buf,5);
    ASSERT_DOUBLE_EQ(buf[1],0.0);
    ASSERT_DOUBLE_EQ(buf[4],0.0);
    ASSERT_DOUBLE_EQ(buf[5],6.0);
    ASSERT_DOUBLE_EQ(buf[9],10.0);
}
