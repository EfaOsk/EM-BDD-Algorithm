#ifndef BDD
#define BDD

#include <stdio.h>
#include "cudd.h"
#include "HMM.h"

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


// Build BDDs
DdNode *build_F_single_seq_O(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M], int observation[T]);
void encode_variables(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M], int **lookup_table_variables);
DdNode **build_F_seq(DdManager *manager, int N, int M, int num_sequences, int T, int **observations, int **lookup_table_variables);


void free_lookup_table_variables(int numVars, int **lookup_table_variables);


// Learn with BDDs
double Backward(DdManager* manager, DdNode* node, const HMM *hmm, int** lookup_table_variables);
void CalculateForward(DdManager* manager, DdNode** F_seq, const HMM *hmm, int T, int num_sequences, int** lookup_table_variables);
void computeConditionalExpectations(DdManager *manager, const HMM *hmm, int T, double ***eta,  double ***gamma,  double *D, int **lookup_table_variables);
HMM* BDD_update(HMM *hmm, double ***eta);
HMM* EMBDD_learn(HMM *hypothesis_hmm, int num_sequences, int **observations, int T, double epsilon);



#endif