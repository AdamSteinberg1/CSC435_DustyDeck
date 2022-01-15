// Demonstrates how to use time.h and gettimeofday function to return wallclock time

#include <stddef.h> /* defines NULL */
#include <sys/time.h>

double factor = 1.0e-6;

extern "C" {
    double walltime_();
    }

double walltime_() 
{
    struct timeval tp;
    int rtn;
    double seconds;

    rtn=gettimeofday(&tp, NULL);

    seconds = tp.tv_sec + factor * tp.tv_usec;

    return  seconds ; 
}
