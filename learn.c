#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cudd.h"
#include "HMM.h"
#include "learn.h"
#include "helpers.h"



/**
 * Builds the F_O BDD (Binary Decision Diagram) for a single sequence.
 * 
 *      F_O = OR_{i=0}^{n_s-1} ( "S_0=i" AND "O_0=o[0]" AND F_O(1,i)" )
 *      
 *      F_O(t, i) = OR_{j=0}^{n_s-1} ( "S_t^(i)=j" AND "O_t^(j)=o[t]" AND F_O(t+1,j)" ) 
 * 
 *      F_O(T, i) = true
 * 
 *
 * @param manager The Cudd Manager.
 * @param N       The number of variables.
 * @param M       The number of variables in AO.
 * @param T       The length of the sequence.
 * @param AS1     Array of BDDs representing AS1.
 * @param AS      Array of BDDs representing AS.
 * @param AO      Array of BDDs representing AO.
 * @param O       Array representing O.
 * @return        The resulting F_O BDD.
 */

DdNode *build_F_single_seq_O(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M], int O[T]) {
    // Create FO array
    DdNode *FO[T+1][N];
    // Base case: F_O^(L, i) = true for all i
    for (int i = 0; i < N; i++) {
        FO[T][i] = Cudd_Not(Cudd_ReadLogicZero(manager));
        Cudd_Ref(FO[T][i]);
    }

    // Recursive case: F_O^(t, i) = (AS[i][t][j] & AO[i][t][0] & F_O^(t+1, j)) for 2 <= t <= L-1
    for (int t = T - 1; t >= 1; t--) {
        for (int i = 0; i < N; i++) {
            FO[t][i] = Cudd_ReadLogicZero(manager);
            Cudd_Ref(FO[t][i]);

            for (int j = 0; j < N; j++) {
                DdNode *temp = Cudd_bddAnd(manager, AS[i][t-1][j], AO[j][t][O[t]]);
                Cudd_Ref(temp);

                DdNode *recursive = Cudd_bddAnd(manager, temp, FO[t+1][j]);
                Cudd_Ref(recursive);

                DdNode *disjunction = Cudd_bddOr(manager, FO[t][i], recursive);
                Cudd_Ref(disjunction);

                Cudd_RecursiveDeref(manager, FO[t][i]);
                Cudd_RecursiveDeref(manager, temp);
                Cudd_RecursiveDeref(manager, recursive);

                FO[t][i] = disjunction;
            }
        }
    }

    // F_O = (AS1[i] & AO[i][1][0] & F_O^(2, i)) for i = 0 to N-1
    for (int i = 0; i < N; i++) {
        DdNode *temp = Cudd_bddAnd(manager, AS1[i], AO[i][0][O[0]]);
        Cudd_Ref(temp);

        DdNode *recursive = Cudd_bddAnd(manager, temp, FO[1][i]);
        Cudd_Ref(recursive);

        FO[0][i] = recursive;

        Cudd_RecursiveDeref(manager, temp);
    }

    DdNode *FO_ = Cudd_ReadLogicZero(manager);
    Cudd_Ref(FO_);
    for (int i = 0; i < N; i++) {
        DdNode *temp = Cudd_bddOr(manager, FO_, FO[0][i]);
        Cudd_Ref(temp);
        Cudd_RecursiveDeref(manager, FO_);
        FO_ = temp;
    }




    for (int t = 0; t < T+1; t++) {
        for (int i = 0; i < N; i++) {
            Cudd_RecursiveDeref(manager, FO[t][i]);
        }
    }

    // Cudd_RecursiveDeref(manager, FO_);
    return FO_;

}


/**
 * Encodes the variables in the given BDD arrays.
 * 
 *
 * @param manager The Cudd Manager.
 * @param N       The number of variables.
 * @param M       The number of variables in AO.
 * @param T       The length of the sequence.
 * @param AS1     Array of BDDs representing AS1.
 * @param AS      Array of BDDs representing AS.
 * @param AO      Array of BDDs representing AO.
 */
void encode_variables(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M], int **lookup_table_variables) {
    // DdNode *AO_enc[N][T][M-1];
    // DdNode *AS1_enc[N - 1];
    // DdNode *AS_enc[N][T - 1][N-1];

    int id = 0;


    // Allocate memory for AS1_enc
    DdNode **AS1_enc = (DdNode **)malloc((N) * sizeof(DdNode *));

    // Encode AS1
    for (int u = 0; u < N - 1; u++) {
        // printf("S_0 = %d\n", u);
        lookup_table_variables[id][0]= 1;
        lookup_table_variables[id][1]= 0;
        lookup_table_variables[id][2]= 0;
        lookup_table_variables[id][3]= u;
        AS1_enc[u] = Cudd_bddIthVar(manager, id++);
    }

    for (int u = 0; u < N; u++) {
        if (u == 0) {
            AS1[0] = AS1_enc[0];
            Cudd_Ref(AS1[0]);
        } 
        else if (u < N-1) {
            AS1[u] = Cudd_Not(Cudd_ReadLogicZero(manager));
            Cudd_Ref(AS1[u]);
            for (int j = 0; j < u; j++){
                DdNode *temp = Cudd_bddAnd(manager, AS1[u], Cudd_Not(AS1_enc[j]));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, AS1[u]);
                AS1[u] = temp;
            }
            DdNode *temp = Cudd_bddAnd(manager, AS1[u], AS1_enc[u]);
            Cudd_Ref(temp);
            Cudd_RecursiveDeref(manager, AS1[u]);
            AS1[u] = temp;
            
        } else {
            AS1[u] = Cudd_Not(Cudd_ReadLogicZero(manager));
            Cudd_Ref(AS1[u]);
            for (int j = 0; j < u; j++){
                DdNode *temp = Cudd_bddAnd(manager, AS1[u], Cudd_Not(AS1_enc[j]));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, AS1[u]);
                AS1[u] = temp;
            }
            
        }
    }
    // Free AS1_enc
    free(AS1_enc);

    // Allocate memory for AS_enc
    DdNode ****AS_enc = (DdNode ****)malloc(N * sizeof(DdNode ***));
    for (int u = 0; u < N; u++) {
        AS_enc[u] = (DdNode ***)malloc(T * sizeof(DdNode **));
        for (int t = 0; t < T - 1; t++) {
            AS_enc[u][t] = (DdNode **)malloc((N) * sizeof(DdNode *));
        }
    }

    // Encode AS
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N-1; v++) {
                // printf("S^%d_%d = %d\n", u, t+1, v);
                lookup_table_variables[id][0]= 2;
                lookup_table_variables[id][1]= u;
                lookup_table_variables[id][2]= t;
                lookup_table_variables[id][3]= v;
                AS_enc[u][t][v] = Cudd_bddIthVar(manager, id++);
            }
        }
    }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N; v++) {
                if (v == 0) {
                    AS[u][t][0] = AS_enc[u][t][0];
                    Cudd_Ref(AS[u][t][0]);
                } 
                else if (v == N-1)
                {
                    AS[u][t][v] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AS[u][t][v]);
                    for (int v0 = 0; v0 < v; v0++){
                        DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], Cudd_Not(AS_enc[u][t][v0]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AS[u][t][v]);
                        AS[u][t][v] = temp;
                    }
                } else {
                    AS[u][t][v] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AS[u][t][v]);
                    for (int v0 = 0; v0 < v; v0++){
                        DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], Cudd_Not(AS_enc[u][t][v0]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AS[u][t][v]);
                        AS[u][t][v] = temp;
                    }

                    DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], AS_enc[u][t][v]);
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(manager, AS[u][t][v]);
                    AS[u][t][v] = temp;
                    
                }
            }
        }
    }
    
    // Free AS_enc
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T - 1; t++) {
            free(AS_enc[u][t]);
        }
        free(AS_enc[u]);
    }
    free(AS_enc);

    
    // Allocate memory for AO_enc
    DdNode ****AO_enc = (DdNode ****)malloc(N * sizeof(DdNode ***));
    for (int u = 0; u < N; u++) {
        AO_enc[u] = (DdNode ***)malloc(T * sizeof(DdNode **));
        for (int t = 0; t < T; t++) {
            AO_enc[u][t] = (DdNode **)malloc((M) * sizeof(DdNode *));
        }
    }

    // Encode AO
    for (int t = 0; t < T; t++) {
        for (int u = 0; u < N; u++) {
            for (int o = 0; o < M-1; o++) {
                // printf("O^%d_%d = %d\n", u, t, o);
                lookup_table_variables[id][0]= 0;
                lookup_table_variables[id][1]= u;
                lookup_table_variables[id][2]= t;
                lookup_table_variables[id][3]= o;
                AO_enc[u][t][o] = Cudd_bddIthVar(manager, id++);
            }
        }
    }
    for (int v = 0; v < N; v++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M; o++) {
                if (o == 0) {
                    AO[v][t][0] = AO_enc[v][t][0];
                    Cudd_Ref(AO[v][t][0]);
                } 
                else if ( o == M-1 )
                {
                    AO[v][t][o] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AO[v][t][o]);
                    for (int o0 = 0; o0 <= M-2; o0++){
                        DdNode *temp = Cudd_bddAnd(manager, AO[v][t][o], Cudd_Not(AO_enc[v][t][o0]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AO[v][t][o]);
                        AO[v][t][o] = temp;
                    }
                    
                } else {
                    AO[v][t][o] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AO[v][t][o]);
                    for (int o0 = 0; o0 <= o-1; o0++){
                        DdNode *temp = Cudd_bddAnd(manager, AO[v][t][o], Cudd_Not(AO_enc[v][t][o0]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AO[v][t][o]);
                        AO[v][t][o] = temp;
                    }
                    DdNode *temp = Cudd_bddAnd(manager, AO[v][t][o], AO_enc[v][t][o]);
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(manager, AO[v][t][o]);
                    AO[v][t][o] = temp;
                    
                }
            }
        }
    }


    // Free AO_enc
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            free(AO_enc[u][t]);
        }
        free(AO_enc[u]);
    }
    free(AO_enc);

}

DdNode *build_C_A(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M]) {


    // L_A = TRUE
    DdNode* L_A = Cudd_Not(Cudd_ReadLogicZero(manager));
    Cudd_Ref(L_A);
    
    for (int t = 0; t < T; t++){
        DdNode* disjunction = Cudd_ReadLogicZero(manager);
        Cudd_Ref(disjunction);
        for (int u = 0; u < N; u++){
            for (int o = 0; o < M; o++){
                DdNode* temp = Cudd_bddOr(manager, AO[u][t][o], disjunction);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, disjunction);
                disjunction = temp;
            }
        }
        DdNode* temp = Cudd_bddAnd(manager, disjunction, L_A);
        Cudd_Ref(temp);
        Cudd_RecursiveDeref(manager, disjunction);
        Cudd_RecursiveDeref(manager, L_A);
        L_A= temp;
    }
    
    DdNode* disjunction = Cudd_ReadLogicZero(manager);
    Cudd_Ref(disjunction);
    for (int u = 0; u < N; u++){
        DdNode* temp = Cudd_bddOr(manager, AS1[u], disjunction);
        Cudd_Ref(temp);
        Cudd_RecursiveDeref(manager, disjunction);
        disjunction = temp;
    }
    DdNode* temp = Cudd_bddAnd(manager, disjunction, L_A);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(manager, disjunction);
    Cudd_RecursiveDeref(manager, L_A);
    L_A= temp;
    
    for (int t = 0; t < T - 1; t++){
        DdNode* disjunction = Cudd_ReadLogicZero(manager);
        Cudd_Ref(disjunction);
        for (int u = 0; u < N; u++){
            for (int v = 0; v < N; v++){
                DdNode* temp = Cudd_bddOr(manager, AS[u][t][v], disjunction);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, disjunction);
                disjunction = temp;
            }
        }
        DdNode* temp = Cudd_bddAnd(manager, disjunction, L_A);
        Cudd_Ref(temp);
        Cudd_RecursiveDeref(manager, disjunction);
        Cudd_RecursiveDeref(manager, L_A);
        L_A= temp;
    }

    DdNode* M_A = Cudd_Not(Cudd_ReadLogicZero(manager));
    Cudd_Ref(M_A);
    for (int t = 0; t < T; t++){
        for (int u = 0; u < N; u++){
            for (int o = 0; o < M; o++){
                for (int u0 = 0; u0 < N; u0++){
                    for (int o0 = 0; o0 < M; o0++){
                        if ( !( (o0 == o) && (u0 == u) ) ) {
                            DdNode* temp0 = Cudd_bddAnd(manager, AO[u][t][o], AO[u0][t][o0]);
                            Cudd_Ref(temp0);
                            DdNode* temp1 = Cudd_Not(temp0);
                            Cudd_Ref(temp1);
                            Cudd_RecursiveDeref(manager, temp0);
                            DdNode* temp2 = Cudd_bddAnd(manager, temp1, M_A);
                            Cudd_Ref(temp2);
                            Cudd_RecursiveDeref(manager, M_A);
                            Cudd_RecursiveDeref(manager, temp1);
                            M_A = temp2;
                        }
                    }
                }
            }
        }
    }


    for (int u = 0; u < N; u++){
        for (int j = 0; j < N; j++){
            if (j != u) {
                DdNode* temp0 = Cudd_bddAnd(manager, AS1[u], AS1[j]);
                Cudd_Ref(temp0);
                DdNode* temp1 = Cudd_Not(temp0);
                Cudd_Ref(temp1);
                Cudd_RecursiveDeref(manager, temp0);
                DdNode* temp2 = Cudd_bddAnd(manager, temp1, M_A);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(manager, M_A);
                Cudd_RecursiveDeref(manager, temp1);
                M_A = temp2;
            }
        }
    }


    for (int t = 0; t < T-1; t++){
        for (int u = 0; u < N; u++){
            for (int v = 0; v < N; v++){
                for (int u0 = 0; u0 < N; u0++){
                    for (int v0 = 0; v0 < N; v0++){
                        if ( !( (u == u0) && (v == v0) ) ) {
                            DdNode* temp0 = Cudd_bddAnd(manager, AS[u][t][v], AS[u0][t][v0]);
                            Cudd_Ref(temp0);
                            DdNode* temp1 = Cudd_Not(temp0);
                            Cudd_Ref(temp1);
                            Cudd_RecursiveDeref(manager, temp0);
                            DdNode* temp2 = Cudd_bddAnd(manager, temp1, M_A);
                            Cudd_Ref(temp2);
                            Cudd_RecursiveDeref(manager, M_A);
                            Cudd_RecursiveDeref(manager, temp1);
                            M_A = temp2;
                        }
                    }
                }
            }
        }
    }

    
    DdNode* C_A = Cudd_bddAnd(manager, M_A, L_A);
    Cudd_Ref(C_A);
    Cudd_RecursiveDeref(manager, M_A);
    Cudd_RecursiveDeref(manager, L_A);

    return C_A;

}


void free_lookup_table_variables(int numVars, int **lookup_table_variables) {
    for (int id = 0; id < numVars; id++) {
        free(lookup_table_variables[id]);
    }
    free(lookup_table_variables);
}

/**
 * Builds the FO nodes for a single sequences O.
 *
 * @param manager The Cudd Manager.
 * @param N       The number of states.
 * @param M       The number of letters in alphabet.
 * @param T       The length of the sequence.
 * @param O       The observation sequence to represent
 *
 * @return An array of FO nodes for all possible sequences.
 */
DdNode **build_F_seq(DdManager *manager, int N, int M, int NO, int T, int **O, int **lookup_table_variables) {
    
    // Define the variables
    DdNode *AS1[N];             // AS1[u] := "S1 = u"
    DdNode *AS[N][T-1][N];      // AS[u][t][v] := "S^u_(t+1) = v"
    DdNode *AO[N][T][M];        // AO[u][t][o] := "O^u_t = o"


    encode_variables(manager, N, M, T, AS1, AS, AO, lookup_table_variables);

    // If not encode, set the initial varables
        // for (int u = 0; u < N; u++){
        //     for (int t = 0; t < T; t++){
        //         for (int o = 0; o < M; o++){
        //             AO[u][t][o] = Cudd_bddNewVar(manager);
        //             printf("O^%d_%d = %d\n", u, t, o);
        //         }
        //     }
        // }

        // for (int u = 0; u < N; u++){
        //     AS1[u] = Cudd_bddNewVar(manager);
        //     printf("S_0 = %d\n", u);
        // }

        // for (int u = 0; u < N; u++){
        //     for (int t = 0; t < T-1; t++){
        //         for (int v = 0; v < N; v++){
        //             AS[u][t][v] = Cudd_bddNewVar(manager);
        //             printf("S^%d_%d = %d\n", u, t+1, v);
        //         }
        //     }
        // }
    
    // If not encode, set the initial varables
    // DdNode* C_A= build_C_A(manager, N, M, T, AS1, AS, AO);
    
    DdNode** F_seq = malloc((NO) * sizeof(DdNode*));
    for (int obs_i = 0; obs_i < NO; obs_i++) {
        DdNode* temp = build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, O[obs_i]);
        F_seq[obs_i] = temp;
    }



    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N; v++) {
                Cudd_RecursiveDeref(manager, AS[u][t][v]);
            }
        }
    }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M; o++) {
                Cudd_RecursiveDeref(manager, AO[u][t][o]);
            }
        }
    }
    for (int u = 0; u < N; u++) {
        Cudd_RecursiveDeref(manager, AS1[u]);
    }

    return F_seq;
}


double get_theta(const HMM *hmm, int x, int i, int j) { 
    if (x == 0) {
        return hmm->B[i][j];
    } else if (x == 1) {
        return hmm->C[j];
    } else if (x == 2) {
        return hmm->A[i][j];
    } else {
        perror("Error: Invalid value of x");
        return -1.0; 
    }
}


double get_sigma(const HMM *hmm, int x, int i, int j) { 
    double sum = 0.0;
    if (x == 0) {
        for (int j0 = j; j0 < hmm->M; j0++) {
            sum += hmm->B[i][j0];
        }
    } else if (x == 1) {
        for (int j0 = j; j0 < hmm->N; j0++){
            sum += hmm->C[j0];
        }
    } else if (x == 2) {
        for (int j0 = j; j0 < hmm->N; j0++) {
            sum += hmm->A[i][j0];
        }
    } else {
        perror("Error: Invalid value of x");
        return -1.0; 
    }
    return sum;
}


double get_prob_encoded(DdManager* manager, const HMM *hmm, DdNode *n, int b, int **lookup_table_variables) { 
    int id = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(n));
    int x = lookup_table_variables[id][0];
    int i = lookup_table_variables[id][1];
    // int t = lookup_table_variables[id][2];
    int j = lookup_table_variables[id][3];

    if (b == 0) {
        return get_sigma(hmm, x, i, j+1) / get_sigma(hmm, x, i, j);
    } else if (b == 1) {
        return get_theta(hmm, x, i, j) / get_sigma(hmm, x, i, j);
    } else {
        perror("Error: Invalid value of b");
        return -1.0; 
    }
}




NodeDataList** nodeData;


NodeDataNode* FindTargetNodeAtLevel(DdManager* manager, int targetLevel, DdNode* targetNode) {
    if (targetLevel == -1) {
        int numVars = Cudd_ReadSize(manager);
        targetLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(targetNode));
        if (targetLevel > numVars) {
            targetLevel = numVars;
        }
    }
    NodeDataNode* current = nodeData[targetLevel]->head;
    while (current != NULL) {
        if (current->node == targetNode) {
            return current; // Return the NodeDataNode containing the target node
        }
        current = current->next;
    }
    printf("Notfound");
    return NULL; // Node not found at the specified level
}



/**
 * @brief Backward procedure where
 * 
 *      /Beta(x, n) = probability that paths logically reach the terminal node x from node n
 * 
 *  for level in vars: // from lowest level to highest
 *      for node n in level:
 *          B_1[n] = 0.5* B_1[Child(True)] + 0.5* B_1[Child(False)]
 * 
 * @param manager 
 * @param node 
 * @param M 
 * @return double
 */
double Backward(DdManager* manager, DdNode* node, const HMM *hmm, int** lookup_table_variables) {
    NodeDataNode *targetNodeData = FindTargetNodeAtLevel(manager, -1, node);
    
    if (targetNodeData->backward >= 0) {
        return targetNodeData->backward; 
    }


    DdNode* high = Cudd_T(node);
    DdNode* low = Cudd_E(node);
    
    double prob_high = get_prob_encoded(manager, hmm, node, 1, lookup_table_variables) * Backward(manager, high, hmm, lookup_table_variables);
    double prob_low;

    if (!Cudd_IsComplement(low)) {
        prob_low = get_prob_encoded(manager, hmm, node, 0, lookup_table_variables) * Backward(manager, low, hmm, lookup_table_variables);
    } else {
        // Adjusting for negative (complemented) edge
        prob_low = get_prob_encoded(manager, hmm, node, 0, lookup_table_variables) * (1.0 - Backward(manager, Cudd_Regular(low), hmm, lookup_table_variables));
    }
    
    // Calculate the probability for the current node
    double prob = prob_low + prob_high;
    targetNodeData->backward = prob; // Store in the lookup table

    return prob;
}


/**
 * Calculate forward values for nodes in a Binary Decision Diagram (BDD).
 *
 * @param manager   The CUDD manager for the BDD operations.
 * @param F_all     An array of BDD nodes representing the roots of the BDDs.
 * @param hmm       A pointer to a Hidden Markov Model (HMM) or a similar probabilistic model.
 * @param num_roots The number of BDD roots in the 'F_all' array.
 */
void CalculateForward(DdManager* manager, DdNode** F_seq, const HMM *hmm, int T, int NO, int** lookup_table_variables) {
    int numVars = Cudd_ReadSize(manager); // Cudd_ReadNodeCount(manager)+(hmm->N)*T+2;
 
    for (int r = 0; r < NO; r++) {
        DdNode *targetNode = F_seq[r];
        int isNegated = Cudd_IsComplement(targetNode);
        if (isNegated) {
            targetNode = Cudd_Regular(targetNode);
        }
        int level = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(targetNode));
        // Find target node
        NodeDataNode* targetNodeData = FindTargetNodeAtLevel(manager, level, targetNode);
        if (targetNodeData == NULL) {
            // raise error, nod not found
            printf("ERROR!");
            return;
        }
        double backwardVal = Backward(manager, targetNode, hmm, lookup_table_variables);
        if (!isNegated) {
            // even number of comple edges
            targetNodeData->forward[1] += 1 / (1 - backwardVal);
        } else {
            // odd number of comple edges
            targetNodeData->forward[0] += 1 / (1 - backwardVal);
        }
    }


    for (int level = 0 ; level < numVars; level++) {
        NodeDataNode* targetNode = nodeData[level]->head;
        while (targetNode != NULL) {
            // Find children:
            DdNode *lowChild = Cudd_E(targetNode->node);
            DdNode *highChild = Cudd_T(targetNode->node);
            
            int lowLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(lowChild));
            if (lowLevel > numVars) {
                lowLevel = numVars;
            }
            int highLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(highChild));
            if (highLevel > numVars) {
                highLevel = numVars;
            }
            NodeDataNode* lowNode = FindTargetNodeAtLevel(manager, lowLevel, Cudd_Regular(lowChild));
            NodeDataNode* highNode = FindTargetNodeAtLevel(manager, highLevel, highChild);

            double ProbLowEdge = get_prob_encoded(manager, hmm, targetNode->node, 0, lookup_table_variables);
            double ProbHighEdge = get_prob_encoded(manager, hmm, targetNode->node, 1, lookup_table_variables);

            highNode->forward[0] += targetNode->forward[0]* ProbHighEdge; 
            highNode->forward[1] += targetNode->forward[1]* ProbHighEdge; 

            if (!Cudd_IsComplement(lowChild)) {
                lowNode->forward[0] += targetNode->forward[0]* ProbLowEdge; 
                lowNode->forward[1] += targetNode->forward[1]* ProbLowEdge; 
            } else {
                lowNode->forward[0] += targetNode->forward[1]* ProbLowEdge; 
                lowNode->forward[1] += targetNode->forward[0]* ProbLowEdge; 
            }

            // next node
            targetNode = targetNode->next;
        }
    }
}


void InitNodeData(DdManager* manager, DdNode** F_seq, int T, int NO) {
    int numVars = Cudd_ReadSize(manager); 
    nodeData = (NodeDataList **)malloc((numVars + 1) * sizeof(NodeDataList*));

    if (nodeData == NULL) {
        // Handle memory allocation failure
        perror("Memory allocation failed");
    }
    for (int i = 0; i < numVars + 1; i++) {
        nodeData[i] = (NodeDataList *)malloc(sizeof(NodeDataList));
        nodeData[i]->head = NULL;
        nodeData[i]->tail = NULL;
    }

    DdNode *node;
    DdGen *gen;


    for (int r = 0; r < NO; r++){
        Cudd_ForeachNode(manager, F_seq[r], gen, node) {
            if (node == Cudd_ReadLogicZero(manager) || node == Cudd_Not(Cudd_ReadLogicZero(manager))) {
                // Terminal node
            } else {
                // Non-terminal nod
                // Get level
                int level = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(node));

                // Check if the node is already in nodeData[level]
                NodeDataNode* current = nodeData[level]->head;
                int nodeExists = 0;

                while (current != NULL) {
                    if (current->node == node) {
                        nodeExists = 1;
                        break;
                    }
                    current = current->next;
                }

                if (!nodeExists) {
                    NodeDataNode* newNode = (NodeDataNode*)malloc(sizeof(NodeDataNode));
                    newNode->node = node;
                    newNode->forward[0] = 0.0; // Inizilize forward with 0
                    newNode->forward[1] = 0.0; // Inizilize forward with 0
                    newNode->backward = -1.0; // Inizilize backward with 0
                    
                    
                    
                    newNode->next = NULL;

                    // Append the new node to the linked list for the level
                    if (nodeData[level]->tail == NULL) {
                        nodeData[level]->head = newNode;
                        nodeData[level]->tail = newNode;
                    } else {
                        nodeData[level]->tail->next = newNode;
                        nodeData[level]->tail = newNode;
                    }
                }
            }
        }
    }

    NodeDataNode* trueTerminal = (NodeDataNode*)malloc(sizeof(NodeDataNode));
    trueTerminal->node = Cudd_Not(Cudd_ReadLogicZero(manager));
    trueTerminal->forward[0] = 0.0;
    trueTerminal->forward[1] = 0.0;
    trueTerminal->backward = 1.0;
    trueTerminal->next = NULL;
    nodeData[numVars]->head = trueTerminal;
    nodeData[numVars]->tail = trueTerminal;
    
}


void CleanNodeData(int numVars) {
    for (int level = 0; level <= numVars; level++) {
        NodeDataNode* node = nodeData[level]->head;
        while (node != NULL) {
            node->forward[0] = 0.0; // Inizilize forward with 0
            node->forward[1] = 0.0; // Inizilize forward with 0
            node->backward = -1.0; // Inizilize backward with 0
            // next node
            node = node->next;
        }
    }
    nodeData[numVars]->head->backward = 1.0;
}

void FreeNodeData(int numVars) {
    // Clean up the allocated memory
    for (int i = 0; i < numVars + 1; i++) {
        NodeDataNode* current = nodeData[i]->head;
        while (current != NULL) {
            NodeDataNode* temp = current;
            current = current->next;
            free(temp);
        }
        free(nodeData[i]);
    }
    free(nodeData);
}


void computeConditionalExpectations(DdManager *manager, const HMM *hmm, int T, double ***eta,  double ***gamma,  double *D, int **lookup_table_variables) {

    int N = hmm->N;
    int NX = (hmm->N > hmm->M) ? hmm->N : hmm->M;
    int numVars = Cudd_ReadSize(manager);

    double e0, e1, temp;

    // Reset gamma and eta
    for (int x = 0; x < 3; x++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < NX; j++) {
                gamma[x][i][j] = 0.0;
                eta[x][i][j] = 0.0;
            }
        }
    }

    // Reset D
    for (int i = 0; i < numVars+1; i++) {
        D[i] = 0.0;
    }

    for (int level = 0; level < numVars; level++) {
        NodeDataNode* node = nodeData[level]->head;
        int x = lookup_table_variables[level][0];
        int i = lookup_table_variables[level][1];
        int t = lookup_table_variables[level][2];
        int j = lookup_table_variables[level][3];
        while (node != NULL) {
            DdNode *lowChild = Cudd_E(node->node);
            DdNode *highChild = Cudd_T(node->node);
                
            int lowLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(Cudd_Regular(lowChild)));
            if (lowLevel > numVars) {
                lowLevel = numVars;
            }
            int highLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(highChild));
            if (highLevel > numVars) {
                highLevel = numVars;
            }
            NodeDataNode *lowChildData = FindTargetNodeAtLevel(manager, lowLevel, Cudd_Regular(lowChild));
            NodeDataNode *highChildData = FindTargetNodeAtLevel(manager, highLevel, highChild);
            double PrLow = get_prob_encoded(manager, hmm, node->node, 0, lookup_table_variables);
            double PrHigh = get_prob_encoded(manager, hmm, node->node, 1, lookup_table_variables);

            e1 = node->forward[0]*PrHigh*(1-highChildData->backward) + node->forward[1]*PrHigh*highChildData->backward;
            
            if (!Cudd_IsComplement(lowChild)) {
                e0 = node->forward[0]*PrLow*(1-lowChildData->backward) + node->forward[1]*PrLow*lowChildData->backward;
            } else {
                e0 = node->forward[0]*PrLow*(lowChildData->backward) + node->forward[1]*PrLow*(1-lowChildData->backward);
            }
            eta[x][i][j] += e1;

            gamma[x][i][j+1] += e0;
            gamma[x][i][j] -= e0 + e1;

            D[level+1] += e0 + e1;
            D[highLevel] -= e1;
            D[lowLevel] -= e0;

            // next node
            node = node->next;
        }
    }

    temp = 0;
    for (int level = 0; level < numVars; level++) {
        int x = lookup_table_variables[level][0];
        int i = lookup_table_variables[level][1];
        int t = lookup_table_variables[level][2];
        int j = lookup_table_variables[level][3];
        temp += D[level];
        if (j == 0) {
            gamma[x][i][0] = temp;
        }
    }

    for (int x = 0; x < 3; x++) { 
        int tempN = N;
        if (x == 1) { tempN = 1; }
        for (int i = 0; i < tempN; i++) {
            temp = 0;
            if (x > 0) { NX = N; } 
            else { NX = hmm->M; }
            for (int j = 0; j < NX; j++) {
                temp += gamma[x][i][j] / get_sigma(hmm, x, i, j); // TODO: / sigma[mu(i)][j];
                eta[x][i][j] += temp * get_theta(hmm, x, i, j);  // TODO: * thetda[mu(i)][j];
            }
        }
    }

}


HMM* update(HMM *hmm, double ***eta) {
    HMM *new_hmm = HMM_create(hmm->N, hmm->M, "Updated model");
    double sum;
    double min_p_f = 0.00001;
    // update b
    for (int u = 0; u < hmm->N; u++) {
        sum = min_p_f*hmm->M;
        // ToDo: remove this quick fix, handles when obs is in alphabet, but not data to learn
        for (int o = 0; o < hmm->M; o++) {
            if (eta[0][u][o]<=0){
                eta[0][u][o] = 0.0;
            }
            sum += eta[0][u][o];
        }
        for (int o = 0; o < hmm->M; o++) {
            new_hmm->B[u][o] = (eta[0][u][o]+min_p_f) / sum;
        }
    }
    // update pi
    sum = hmm->N*min_p_f;
    for (int v = 0; v < hmm->N; v++) {
        sum += eta[1][0][v];
    }
    for (int v = 0; v < hmm->N; v++) {
        new_hmm->C[v] = (eta[1][0][v]+min_p_f) / sum;
    }
    // update a
    for (int u = 0; u < hmm->N; u++) {
        sum = hmm->N*min_p_f;
        for (int v = 0; v < hmm->N; v++) {
            sum += eta[2][u][v];
        }
        for (int v = 0; v < hmm->N; v++) {
            new_hmm->A[u][v] = (eta[2][u][v]+min_p_f) / sum;
        }
    }
    return new_hmm;
}


/**
 * @brief This function builds an HMM that "best" fits a given observation sequence O with the given number of states and observations.
 * 
 *      Uses the cour structure of the Baum-Welch algorithm, and BDDs, in order to improve the time and memory complexity
 *       
 * @param hypothesis_hmm  Pointer to an HMM structure, which is the initial hypothesis for the HMM.
 * @param T               The length of the observation sequence.
 * @param O               Observation sequence, an array of integers of length T.
 * @param epsilon         The convergence threshold; when the change in log likelihood is less than or equal to this value, learning stops.
 * @param logs_folder     Path to the folder where logs and model files will be saved.
 * @return HMM*           Pointer to the learned HMM structure.
 */
HMM* learn(HMM *hypothesis_hmm, int T, int NO, int **O, double epsilon, const char *logs_folder, const char *result_file)
{
    /*

    EM on BDD
        (1) Build (S)BDD
        (1.5) Ordered encoding
        
        (2) initilize M (= some random HMM) 

        (repeat steps 3-5 until converged)

        (3) E-step
            (a) Backward
            (b) Forward
            (c) Conditional expectations
        
        (4) M-step
            (a) update M
        
        (5) Calculate the log-likelyhood of M

        (6) retrun M

    */
    
    int N = hypothesis_hmm->N;
    int M = hypothesis_hmm->M;
    HMM *model = HMM_create(N, M, "model");

    HMM_copy(model, hypothesis_hmm);

    // for loging
    int iteration = 0;
    char log_filename[256];
    char model_filename[256];

    // Construct the log filename
    // sprintf(log_filename, "%s/log.txt", logs_folder);
    sprintf(log_filename, "log.txt");

    // Open the log file
    FILE *log_file = fopen(log_filename, "a");
    if (log_file == NULL) {
        perror("Error opening log file");
        return NULL;
    }

    sprintf(model_filename, "%s/models", logs_folder);
    mkdir(model_filename, 0777);
    

    // Step 1 build (S)BDD
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0);  
    // Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);
    
    int **lookup_table_variables;
    lookup_table_variables = (int **)malloc((N*T*(M-1)+(N-1)+(N*(T-1)*(N-1))) * sizeof(int *));
    for (int id = 0; id <(N*T*(M-1)+(N-1)+(N*(T-1)*(N-1))); id++) {
        (lookup_table_variables)[id] = (int *)malloc(4 * sizeof(int));
    }

    DdNode **F_obs = build_F_seq(manager, N, M, NO, T, O, lookup_table_variables);

    int numVars = Cudd_ReadSize(manager);

    // Step 2: initilize M (= some random HMM) 
        // ToDo currently input


    double prob_priv, prob_original, prob_new;
    prob_original = log_likelihood_forward_multiple(model, O, NO, T);
    prob_priv = prob_original;
    int converged = 0;

    int tmp = (N > M) ? N : M;
    double ***eta = allocate_3D_matrix(3, N, tmp, 0.0);
    double ***gamma = allocate_3D_matrix(3, N, tmp, 0.0);
    double *D = malloc((numVars+1) * sizeof(double));

    InitNodeData(manager, F_obs, T, NO);

    while (!converged)
    {
        clock_t start_time = clock();

        // Step 3: E-step
        
        CleanNodeData(numVars);
    
        // Step 3 (a) : Backward

        // Step 3 (b) : Forward
    
        CalculateForward(manager, F_obs, model, T, NO, lookup_table_variables);

        // Step 3 (c) : Conditional Expectations

        computeConditionalExpectations(manager, model, T, eta, gamma, D, lookup_table_variables);
    
        // Step 4: M-step
        // Step 4 (a) : update M

        const HMM *new_hmm = update(model, eta);


        prob_new = log_likelihood_forward_multiple(new_hmm, O, NO, T);
        
        
        HMM_copy(model, new_hmm); 
        HMM_destroy(new_hmm);

        // HMM_print(model);
        HMM_validate(model);

        clock_t end_time = clock();
        double iteration_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        printf("\timprovement: %f\n", prob_new-prob_priv);
        fprintf(log_file, "Iteration: %d, Log Likelihood: %f, Improvement: %f, Time: %f\n",
                iteration, prob_new, prob_new-prob_priv, iteration_time);
        fflush(log_file); 

        sprintf(model_filename, "%s/models/model_%d", logs_folder, iteration);
        HMM_save(model, model_filename); 
        if (prob_new <= prob_priv+epsilon) {
            converged = 1;
        }
        prob_priv = prob_new;
        iteration++;
    }

    fclose(log_file);

    // Open the result file in append mode
    FILE *result_fp = fopen(result_file, "a");
    if (result_fp == NULL) {
        perror("Error opening result file");
        return NULL;
    }

    // Append the results to the result file
    fprintf(result_fp, "%d, %f, %f\n", iteration, prob_new, prob_new - prob_original);


    // ToDo: remove (For Debuging):
    fprintf(result_fp, "N : %d | ", N );
    fprintf(result_fp, "M : %d | ", M );
    fprintf(result_fp, "T : %d | ", T );
    fprintf(result_fp, "Encode: TRUE | ");
    fprintf(result_fp, "DdManager vars: %d | ", numVars ); // Returns the number of BDD variables in existence
    fprintf(result_fp, "DdManager nodes: %ld | ", Cudd_ReadNodeCount(manager) ); // Reports the number of live nodes in BDDs and ADDs
    fprintf(result_fp, "DdManager reorderings: %d | ", Cudd_ReadReorderings(manager) ); // Returns the number of times reordering has occurred
    fprintf(result_fp, "DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager) ); // Returns the memory in use by the manager measured in bytes
    
    // Close the result file
    fclose(result_fp);

    char BDDfilename[526];
    sprintf(BDDfilename, "%s/BDD.dot", logs_folder); // Write .dot filename to a string
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(BDDfilename,"w");
    Cudd_DumpDot(manager, (1), F_obs, NULL, NULL, outfile);
    fclose(outfile);

    // Clean up
    for (int x = 0; x < 3; x++) {
        for (int i = 0; i < N; i++) {
            free(gamma[x][i]);
            free(eta[x][i]);
        }
        free(gamma[x]);
        free(eta[x]);
    }
    free(gamma);
    free(eta);

    free(D);


    free(F_obs);
    FreeNodeData(numVars);
    free_lookup_table_variables(numVars, lookup_table_variables);
    Cudd_Quit(manager);

    // Step 6: Return the learned model

    return model;
}
