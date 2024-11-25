/*
 */

#include "ljmd.h"

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <stdint.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

/* helper function: get current time in seconds since epoch */

double wallclock()
{
#if defined(_WIN32)
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTimeAsFileTime(&file_time);
    time = ((uint64_t) file_time.dwHighDateTime) << 32;
    time |= (uint64_t) file_time.dwLowDateTime;
    return ((double) time)*1.0e-7;
#else
    struct timeval t;
    gettimeofday(&t,0);
    return ((double) t.tv_sec) + 1.0e-6*((double) t.tv_usec);
#endif
}

/* helper function: sleep for some time with support for fractions of a second */

void doublesleep(double t)
{
#if defined(_WIN32)
    DWORD msec = (DWORD) (t*1000.0);
    Sleep(msec);
#else
    struct timespec ts;
    ts.tv_sec = (time_t) t;
    ts.tv_nsec = (long)((t-(double)ts.tv_sec)*1000000000.0);
    nanosleep(&ts, NULL);
#endif
}

/* helper function: zero out an array of doubles */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}
