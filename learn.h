#ifndef LEARN
#define LEARN

#include <stdio.h>
#include "cudd.h"
#include "HMM.h"


HMM* learn(HMM *orignal_hmm, int T, int O[T]);

typedef struct NodeDataNode {
    DdNode* node;
    double forward[2];  // Forward probabilities for binary values 0 and 1
    double backward[2]; // Backward probabilities for binary values 0 and 1
    struct NodeDataNode* next;
} NodeDataNode;

typedef struct NodeDataList {
    NodeDataNode* head;
    NodeDataNode* tail;
} NodeDataList;


#endif