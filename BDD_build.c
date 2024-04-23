#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "cudd.h"
#include "HMM.h"
#include "BDD.h"
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
 * @param observation   Array representing O.
 * @return        The resulting F_O BDD.
 */
DdNode *build_F_single_seq_O(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][N], DdNode *AO[N][M], int observations[T]) {
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
                DdNode *temp = Cudd_bddAnd(manager, AS[i][j], AO[j][observations[t]]);
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
        DdNode *temp = Cudd_bddAnd(manager, AS1[i], AO[i][observations[0]]);
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
void encode_variables(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][N], DdNode *AO[N][M], int **lookup_table_variables) {
    // DdNode *AO_enc[N][M-1];
    // DdNode *AS1_enc[N - 1];
    // DdNode *AS_enc[N][N-1];

    int id = 0;


    // Allocate memory for AS1_enc
    DdNode **AS1_enc = (DdNode **)malloc((N) * sizeof(DdNode *));

    // Encode AS1
    for (int u = 0; u < N - 1; u++) {
        // printf("S_0 = %d\n", u);
        lookup_table_variables[id][0]= 1;
        lookup_table_variables[id][1]= 0;
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
    DdNode ***AS_enc = (DdNode ***)malloc(N * sizeof(DdNode **));
    for (int u = 0; u < N; u++) {
        AS_enc[u] = (DdNode **)malloc((N) * sizeof(DdNode *));
    }

    // Encode AS
    for (int u = 0; u < N; u++) {
        for (int v = 0; v < N-1; v++) {
            // printf("S^%d_%d = %d\n", u, v);
            lookup_table_variables[id][0]= 2;
            lookup_table_variables[id][1]= u;
            lookup_table_variables[id][3]= v;
            AS_enc[u][v] = Cudd_bddIthVar(manager, id++);
        }
    }
    for (int u = 0; u < N; u++) {
        for (int v = 0; v < N; v++) {
            if (v == 0) {
                AS[u][0] = AS_enc[u][0];
                Cudd_Ref(AS[u][0]);
            } 
            else if (v == N-1)
            {
                AS[u][v] = Cudd_Not(Cudd_ReadLogicZero(manager));
                Cudd_Ref(AS[u][v]);
                for (int v0 = 0; v0 < v; v0++){
                    DdNode *temp = Cudd_bddAnd(manager, AS[u][v], Cudd_Not(AS_enc[u][v0]));
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(manager, AS[u][v]);
                    AS[u][v] = temp;
                }
            } else {
                AS[u][v] = Cudd_Not(Cudd_ReadLogicZero(manager));
                Cudd_Ref(AS[u][v]);
                for (int v0 = 0; v0 < v; v0++){
                    DdNode *temp = Cudd_bddAnd(manager, AS[u][v], Cudd_Not(AS_enc[u][v0]));
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(manager, AS[u][v]);
                    AS[u][v] = temp;
                }

                DdNode *temp = Cudd_bddAnd(manager, AS[u][v], AS_enc[u][v]);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, AS[u][v]);
                AS[u][v] = temp;
                
            }
        }
    }

    // Free AS_enc
    for (int u = 0; u < N; u++) {
        free(AS_enc[u]);
    }
    free(AS_enc);

    
    // Allocate memory for AO_enc
    DdNode ***AO_enc = (DdNode ***)malloc(N * sizeof(DdNode **));
    for (int u = 0; u < N; u++) {
        AO_enc[u] = (DdNode **)malloc(M * sizeof(DdNode *));
    }

    // Encode AO
    for (int u = 0; u < N; u++) {
        for (int o = 0; o < M-1; o++) {
            // printf("O^%d_%d = %d\n", u, t, o);
            lookup_table_variables[id][0]= 0;
            lookup_table_variables[id][1]= u;
            lookup_table_variables[id][3]= o;
            AO_enc[u][o] = Cudd_bddIthVar(manager, id++);
        }
    }
    for (int v = 0; v < N; v++) {
        for (int o = 0; o < M; o++) {
            if (o == 0) {
                AO[v][0] = AO_enc[v][0];
                Cudd_Ref(AO[v][0]);
            } 
            else if ( o == M-1 )
            {
                AO[v][o] = Cudd_Not(Cudd_ReadLogicZero(manager));
                Cudd_Ref(AO[v][o]);
                for (int o0 = 0; o0 <= M-2; o0++){
                    DdNode *temp = Cudd_bddAnd(manager, AO[v][o], Cudd_Not(AO_enc[v][o0]));
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(manager, AO[v][o]);
                    AO[v][o] = temp;
                }
                
            } else {
                AO[v][o] = Cudd_Not(Cudd_ReadLogicZero(manager));
                Cudd_Ref(AO[v][o]);
                for (int o0 = 0; o0 <= o-1; o0++){
                    DdNode *temp = Cudd_bddAnd(manager, AO[v][o], Cudd_Not(AO_enc[v][o0]));
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(manager, AO[v][o]);
                    AO[v][o] = temp;
                }
                DdNode *temp = Cudd_bddAnd(manager, AO[v][o], AO_enc[v][o]);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, AO[v][o]);
                AO[v][o] = temp;
                
            }
        }
    }


    // Free AO_enc
    for (int u = 0; u < N; u++) {
        free(AO_enc[u]);
    }
    free(AO_enc);

}


/**
 * Builds the FO nodes for a single sequences O.
 *
 * @param manager The Cudd Manager.
 * @param N       The number of states.
 * @param M       The number of letters in alphabet.
 * @param T       The length of the sequence.
 * @param observations       The observation sequence to represent
 *
 * @return An array of FO nodes for all possible sequences.
 */
DdNode **build_F_seq(DdManager *manager, int N, int M, int NO, int T, int **observations, int **lookup_table_variables) {
    // Define the variables
    DdNode *AS1[N];             // AS1[u] := "S1 = u"
    DdNode *AS[N][N];      // AS[u][t][v] := "S^u_(t+1) = v"
    DdNode *AO[N][M];        // AO[u][t][o] := "O^u_t = o"

    encode_variables(manager, N, M, T, AS1, AS, AO, lookup_table_variables);

    DdNode** F_seq = malloc((NO) * sizeof(DdNode*));
    for (int obs_i = 0; obs_i < NO; obs_i++) {
        DdNode* temp = build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, observations[obs_i]);
        F_seq[obs_i] = temp;
    }




    for (int u = 0; u < N; u++) {
        for (int v = 0; v < N; v++) {
            Cudd_RecursiveDeref(manager, AS[u][v]);
        }
    }
    for (int u = 0; u < N; u++) {
        for (int o = 0; o < M; o++) {
            Cudd_RecursiveDeref(manager, AO[u][o]);
        }
    }
    for (int u = 0; u < N; u++) {
        Cudd_RecursiveDeref(manager, AS1[u]);
    }
    return F_seq;
}


void free_lookup_table_variables(int numVars, int **lookup_table_variables) {
    for (int id = 0; id < numVars; id++) {
        free(lookup_table_variables[id]);
    }
    free(lookup_table_variables);
}