// Demonstrates how to use time.h to return cputime

#include <time.h>

extern "C" {
    double cputime_();
    }


double cputime_() 
{
    return  (double) clock() / (double) CLOCKS_PER_SEC;
}
