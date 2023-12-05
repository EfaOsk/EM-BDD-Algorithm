#ifndef LEARN
#define LEARN

#include <stdio.h>
#include "cudd.h"
#include "HMM.h"


HMM* learn(HMM *hypothesis_hmm, int T, int NO, int **O, double epsilon, const char *logs_folder, const char *result_file);

typedef struct NodeDataNode {
    DdNode* node;
    double forward[2];  // Forward probabilities for binary values 0 and 1
    double backward; // Backward probabilities for binary values 0 and 1
    struct NodeDataNode* next;
} NodeDataNode;

typedef struct NodeDataList {
    NodeDataNode* head;
    NodeDataNode* tail;
} NodeDataList;


#endif