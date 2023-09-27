#ifndef LEARN
#define LEARN

#include <stdio.h>
#include "cudd.h"
#include "HMM.h"


HMM* learn(HMM *orignal_hmm, int T, int O[T]);


#endif