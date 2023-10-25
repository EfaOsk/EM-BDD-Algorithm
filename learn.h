#ifndef LEARN
#define LEARN

#include <stdio.h>
#include "cudd.h"
#include "HMM.h"


HMM* learn(HMM *orignal_hmm, int T, int O[T]);

typedef struct {
    DdNode *node;
    double forward0;
    double forward1;
} NodeWithForward;


#endif