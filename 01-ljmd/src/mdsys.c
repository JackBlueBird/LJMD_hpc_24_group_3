/* a few physical constants */
#include "mdsys.h"
#if defined(_MPI)
    #include <mpi.h>
#endif

const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */
