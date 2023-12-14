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
    
    for (int t = 0; t < T; t++) {
        for (int u = 0; u < N; u++){
            DdNode* disjunction = Cudd_ReadLogicZero(manager);
            Cudd_Ref(disjunction);
            for (int o = 0; o < M; o++){
                DdNode* temp = Cudd_bddOr(manager, AO[u][t][o], disjunction);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, disjunction);
                disjunction = temp;
            }
            DdNode* temp = Cudd_bddAnd(manager, disjunction, L_A);
            Cudd_Ref(temp);
            Cudd_RecursiveDeref(manager, disjunction);
            Cudd_RecursiveDeref(manager, L_A);
            L_A= temp;
        }
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
        for (int u = 0; u < N; u++){
            DdNode* disjunction = Cudd_ReadLogicZero(manager);
            Cudd_Ref(disjunction);
            for (int v = 0; v < N; v++){
                DdNode* temp = Cudd_bddOr(manager, AS[u][t][v], disjunction);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, disjunction);
                disjunction = temp;
            }
            DdNode* temp = Cudd_bddAnd(manager, disjunction, L_A);
            Cudd_Ref(temp);
            Cudd_RecursiveDeref(manager, disjunction);
            Cudd_RecursiveDeref(manager, L_A);
            L_A= temp;
        }
    }

    DdNode* M_A = Cudd_Not(Cudd_ReadLogicZero(manager));
    Cudd_Ref(M_A);
    for (int t = 0; t < T; t++){
        for (int u = 0; u < N; u++){
            for (int o = 0; o < M; o++){
                for (int o0 = 0; o0 < M; o0++){
                    if ( !(o0 == o)) {
                        DdNode* temp0 = Cudd_bddAnd(manager, AO[u][t][o], AO[u][t][o0]);
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
                for (int v0 = 0; v0 < N; v0++){
                    if ( !(v == v0)) {
                        DdNode* temp0 = Cudd_bddAnd(manager, AS[u][t][v], AS[u][t][v0]);
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

    
    DdNode* C_A = Cudd_bddAnd(manager, M_A, L_A);
    Cudd_Ref(C_A);
    Cudd_RecursiveDeref(manager, M_A);
    Cudd_RecursiveDeref(manager, L_A);

    return C_A;

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
        //             // printf("O^%d_%d = %d\n", u, t, o);
        //         }
        //     }
        // }

        // for (int u = 0; u < N; u++){
        //     AS1[u] = Cudd_bddNewVar(manager);
        //     // printf("S_0 = %d\n", u);
        // }

        // for (int u = 0; u < N; u++){
        //     for (int t = 0; t < T-1; t++){
        //         for (int v = 0; v < N; v++){
        //             AS[u][t][v] = Cudd_bddNewVar(manager);
        //             // printf("S^%d_%d = %d\n", u, t+1, v);
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


void free_lookup_table_variables(int numVars, int **lookup_table_variables) {
    for (int id = 0; id < numVars; id++) {
        free(lookup_table_variables[id]);
    }
    free(lookup_table_variables);
}