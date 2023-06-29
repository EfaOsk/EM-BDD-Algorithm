#ifndef LEARN
#define LEARN

#include <stdio.h>
#include "cudd.h"
#include "HMM.h"


struct HMM learn(const int N, const int M, int T, int O[T]);


#endif