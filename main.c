#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "cudd.h"

#include "HMM.h"
#include "learn.h"



int main (int argc, char *argv[])
{
    const int N = 2;
    const int M = 2;


    // The observation O
    const int T = 2;
    int O[T];
    O[0]= 0;
    O[1]= 1;


    struct HMM L = learn(N, M, T, O);

    return 0;
}