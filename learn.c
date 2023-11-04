#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "cudd.h"
#include "HMM.h"
#include "learn.h"
#include "helpers.h"



int **lookup_table_variables;

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

DdNode **build_F_single_seq_O(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M], int O[T]) {
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




    DdNode **flattenedArray = (DdNode **)malloc((T) * N * sizeof(DdNode *));

    // Flatten the 2D array into the 1D array.
    int count = 0;
    for (int t = 0; t < T; t++) {
        for (int n = 0; n < N; n++) {
            flattenedArray[count] = FO[t][n]; // Cudd_BddToAdd(manager, FO[t][n]);
            count++;
        }
    }

    // for (int t = 0; t < T+1; t++) {
    //     for (int i = 0; i < N; i++) {
    //         Cudd_RecursiveDeref(manager, FO[t][i]);
    //     }
    // }
    Cudd_RecursiveDeref(manager, FO_);
    return flattenedArray;

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
void encode_variables(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M]) {
    // DdNode *AO_enc[N][T][M-1];
    // DdNode *AS1_enc[N - 1];
    // DdNode *AS_enc[N][T - 1][N-1];

    // Allocate memory for AO_enc
    DdNode ****AO_enc = (DdNode ****)malloc(N * sizeof(DdNode ***));
    for (int u = 0; u < N; u++) {
        AO_enc[u] = (DdNode ***)malloc(T * sizeof(DdNode **));
        for (int t = 0; t < T; t++) {
            AO_enc[u][t] = (DdNode **)malloc((M) * sizeof(DdNode *));
        }
    }

    // Allocate memory for AS1_enc
    DdNode **AS1_enc = (DdNode **)malloc((N) * sizeof(DdNode *));

    // Allocate memory for AS_enc
    DdNode ****AS_enc = (DdNode ****)malloc(N * sizeof(DdNode ***));
    for (int u = 0; u < N; u++) {
        AS_enc[u] = (DdNode ***)malloc(T * sizeof(DdNode **));
        for (int t = 0; t < T - 1; t++) {
            AS_enc[u][t] = (DdNode **)malloc((N) * sizeof(DdNode *));
        }
    }

    int id = 0;

    // Encode AO
    for (int t = 0; t < T; t++) {
        for (int u = 0; u < N; u++) {
            for (int o = 0; o < M; o++) {
                if ((o == M-1) && (u == N-1)){

                } else {
                    // printf("O^%d_%d = %d\n", u, t, o);
                    lookup_table_variables[id][0]= 0;
                    lookup_table_variables[id][1]= u;
                    lookup_table_variables[id][2]= t;
                    lookup_table_variables[id][3]= o;
                    AO_enc[u][t][o] = Cudd_bddIthVar(manager, id++);
                }
            }
        }
    }
    for (int v = 0; v < N; v++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M; o++) {
                if ((v == 0) && ((o == 0))) {
                    AO[0][t][0] = AO_enc[0][t][0];
                    Cudd_Ref(AO[0][t][0]);
                } 
                else if ( (v == N-1) && (o == M-1) )
                {
                    AO[v][t][o] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AO[v][t][o]);
                    for (int v0 = 0; v0 <= N-2; v0++){
                        for (int o0 = 0; o0 <= M-1; o0++){
                            DdNode *temp = Cudd_bddAnd(manager, AO[v][t][o], Cudd_Not(AO_enc[v0][t][o0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AO[v][t][o]);
                            AO[v][t][o] = temp;
                        }
                    }
                    for (int o0 = 0; o0 <= M-2; o0++){
                        DdNode *temp = Cudd_bddAnd(manager, AO[v][t][o], Cudd_Not(AO_enc[v][t][o0]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AO[v][t][o]);
                        AO[v][t][o] = temp;
                    }
                    
                } else {
                    

                    AO[v][t][o] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AO[v][t][o]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int o0 = 0; o0 <= M-1; o0++){
                            DdNode *temp = Cudd_bddAnd(manager, AO[v][t][o], Cudd_Not(AO_enc[v0][t][o0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AO[v][t][o]);
                            AO[v][t][o] = temp;
                        }
                    }
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


    // Encode AS
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N; v++) {
                if ((v == N-1) && (u == N-1)){

                } else {
                    // printf("S^%d_%d = %d\n", u, t+1, v);
                    lookup_table_variables[id][0]= 2;
                    lookup_table_variables[id][1]= u;
                    lookup_table_variables[id][2]= t;
                    lookup_table_variables[id][3]= v;
                    AS_enc[u][t][v] = Cudd_bddIthVar(manager, id++);
                }
            }
        }
    }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N; v++) {
                if ((v == 0) && ((u == 0))) {
                    AS[0][t][0] = AS_enc[0][t][0];
                    Cudd_Ref(AS[0][t][0]);
                } 
                else if ( (v == N-1) && (u == N-1) )
                {
                    AS[u][t][v] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AS[u][t][v]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int u0 = 0; u0 <= N-1; u0++){
                            DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AS[u][t][v]);
                            AS[u][t][v] = temp;
                        }
                    }
                    for (int u0 = 0; u0 < u; u0++){
                        DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AS[u][t][v]);
                        AS[u][t][v] = temp;
                    }
                    
                } else {

                    AS[u][t][v] = Cudd_Not(Cudd_ReadLogicZero(manager));
                    Cudd_Ref(AS[u][t][v]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int u0 = 0; u0 <= N-1; u0++){
                            DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AS[u][t][v]);
                            AS[u][t][v] = temp;
                        }
                    }
                    for (int u0 = 0; u0 < u; u0++){
                        DdNode *temp = Cudd_bddAnd(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v]));
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

    // Free AO_enc
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            free(AO_enc[u][t]);
        }
        free(AO_enc[u]);
    }
    free(AO_enc);

    // Free AS_enc
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T - 1; t++) {
            free(AS_enc[u][t]);
        }
        free(AS_enc[u]);
    }
    free(AS_enc);

    // Free AS1_enc
    free(AS1_enc);

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

void free_lookup_table_variables(int numVars) {
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
DdNode **build_F_seq(DdManager *manager, int N, int M, int T, int O[T]) {
    
    // Define the variables
    DdNode *AS1[N];             // AS1[u] := "S1 = u"
    DdNode *AS[N][T-1][N];      // AS[u][t][v] := "S^u_(t+1) = v"
    DdNode *AO[N][T][M];        // AO[u][t][o] := "O^u_t = o"


    encode_variables(manager, N, M, T, AS1, AS, AO);

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
    

    DdNode** F_seq = malloc(N*(T+1) * sizeof(DdNode*));

    DdNode** temp = build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, O);
    // F_seq =  Cudd_BddToAdd(manager, Cudd_bddAnd(manager, C_A, temp)); // If not encode
    F_seq =  temp;
    // Cudd_Ref(F_seq);
    // Cudd_RecursiveDeref(manager, temp);


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


/**
 * @brief Calculate the edge probability in a SBDD reprentation of a HMM.
 *
 * This function calculates the probability of transitioning to an encoded state in a SBDD 
 * reprentation of a HMM based on the provided parameters.
 *
 * @param hmm     A pointer to the HMM structure.
 * @param i       The current state index.
 * @param b       A flag indicating whether it's a true edge (b=0) or false edge (b=1).
 *
 * @return        The probability of transitioning to the encoded state.
 *
 * @note          The function uses the following formulas:
 *                - For a true edge (b=0):
 *                  \Prob{\langle c_i^{enc}, true\rangle} = \frac{\pi(i)}{\sum_{i'=i}^{N-1}\pi(i')}
 *                - For a false edge (b=1):
 *                  \Prob{\langle c_i^{enc}, false\rangle} = 
 *                      \frac{
 *                          \sum_{i'=i+1}^{N-1}\pi(i')
 *                      }{
 *                          \sum_{i'=i}^{N-1}\pi(i')
 *                      }
 */
double get_prob_AS1_encoded(const HMM *hmm, int i, int b){

    
    if (b==0){ // false edge
        double sum = 0.0;

        for (int i0 = i; i0 < hmm->N; i0++){
            sum += hmm->C[i0];
        }
        double sum_ = 0.0;

        for (int i0 = i+1; i0 < hmm->N; i0++){
            sum_ += hmm->C[i0];
        }

        // Sum_{i'=i+1}^{i=N-1} pi(i'/ Sum_{i'=i}^{i=N-1} pi(i') 
        return sum_ / sum ;
    } else { // true edge
        double sum = 0.0;

        for (int i0 = i; i0 < hmm->N; i0++){
            sum += hmm->C[i0];
        }

        //pi(i)/ Sum_{i'=i}^{i=N-1} pi(i') 
        return (hmm->C[i]) / sum ;
    }

}


/**
 * @brief Calculate the edge probability in a SBDD reprentation of a HMM.
 *
 * This function calculates the probability of transitioning to an encoded state in a SBDD 
 * reprentation of a HMM based on the provided parameters.
 *
 * @param hmm     A pointer to the HMM structure.
 * @param i       The previus state index.
 * @param j       The current state index.
 * @param b       A flag indicating whether it's a true edge (b=0) or false edge (b=1).
 *
 * @return        The probability of transitioning to the encoded state.
 *
 * @note          The function uses the following formulas:
 *                - For a true edge (b=0):
 *                  \Prob{\langle d_{i,j}^{enc}, true\rangle} = 
 *                      \frac{
 *                          a(i)(j)
 *                      }{
 *                          (\sum_{i'=i+1}^{N-1}\sum_{j'=0}^{N-1}   a(i')(j'))+(\sum_{j'=j}^{N-1} a(i)(j'))
 *                      }
 *                - For a false edge (b=1):
 *                  \Prob{\langle d_i^{enc}, false\rangle} = 
 *                      \frac{
 *                          (\sum_{i'=i+1}^{N-1}\sum_{j'=0}^{N-1} a(i')(j'))+ (\sum_{j'=j+1}^{N-1} a(i)(j')
 *                      }{
 *                          \sum_{i'=i+1}^{N-1}\sum_{j'=0}^{N-1} a(i')(j'))+(\sum_{j'=j}^{N-1} a(i)(j'))
 *                      }
 */
double get_prob_AS_encoded(const HMM *hmm, int i, int j, int b){

    if (b==0){ // false edge
        double sum_num = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++){
            for (int j0 = 0; j0<= hmm->N -1; j0++){
                sum_num += hmm->A[i0][j0];
            }
        }
        for (int j0 = j+1; j0<= hmm->N -1; j0++){
            sum_num += hmm->A[i][j0];
        }

        double sum_den = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++){
            for (int j0 = 0; j0<= hmm->N -1; j0++){
                sum_den += hmm->A[i0][j0];
            }
        }
        for (int j0 = j; j0<= hmm->N -1; j0++){
            sum_den += hmm->A[i][j0];
        }

        // 
        return sum_num / sum_den ;
    } else { // true edge
        double sum_den = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++){
            for (int j0 = 0; j0<= hmm->N -1; j0++){
                sum_den += hmm->A[i0][j0];
            }
        }
        for (int j0 = j; j0<= hmm->N -1; j0++){
            sum_den += hmm->A[i][j0];
        }

        // a(i)(j) / 
        return (hmm->A[i][j]) / sum_den ;
    }

}


/**
 * @brief Calculate the edge probability in a SBDD reprentation of a HMM.
 *
 * This function calculates the probability of transitioning to an encoded state in a SBDD 
 * reprentation of a HMM based on the provided parameters.
 *
 * @param hmm     A pointer to the HMM structure.
 * @param i       The current state index.
 * @param j       The current observation index.
 * @param b       A flag indicating whether it's a true edge (b=0) or false edge (b=1).
 *
 * @return        The probability of transitioning to the encoded state.
 *
 * @note          The function uses the following formulas:
 *                - For a true edge (b=0):
 *                  \Prob{\langle e_{i,j}^{enc}, true\rangle} = 
 *                      \frac{
 * \                        b(i)(j)
 *                      }{
 *                          (\sum_{i'=i+1}^{N-1}\sum_{j'=0}^{M-1} b(i')(j'))+(\sum_{j'=j}^{M-1} b(i)(j'))
 *                      }
 *                - For a false edge (b=1):
 *                  \Prob{\langle e_i^{enc}, false\rangle} = 
 *                      \frac{
 *                          (\sum_{i'=i+1}^{N-1}\sum_{j'=0}^{M-1} b(i')(j'))+ (\sum_{j'=j+1}^{M-1} b(i)(j')
 *                      }{
 *                          \sum_{i'=i+1}^{N-1}\sum_{j'=0}^{M-1} b(i')(j'))+(\sum_{j'=j}^{M-1} b(i)(j'))
 *                      }
 */
double get_prob_AO_encoded(const HMM *hmm, int i, int j, int edgeType) {

    // Calculate sum_denominator which is common for both edgeType
    double sum_denominator = 0.0;

    for (int i0 = i+1; i0 <= hmm->N - 1; i0++) {
        for (int j0 = 0; j0 <= hmm->M - 1; j0++) {
            sum_denominator += hmm->B[i0][j0];
        }
    }
    for (int j0 = j; j0 <= hmm->M - 1; j0++) {
        sum_denominator += hmm->B[i][j0];
    }

    if (edgeType == 0) { // false edge
        double sum_numerator = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++) {
            for (int j0 = 0; j0 <= hmm->M - 1; j0++) {
                sum_numerator += hmm->B[i0][j0];
            }
        }
        for (int j0 = j+1; j0 <= hmm->M - 1; j0++) {
            sum_numerator += hmm->B[i][j0];
        }

        return sum_numerator / sum_denominator;

    } else { // true edge
        return hmm->B[i][j] / sum_denominator;
    }
}

double get_prob_encoded(const HMM *hmm, DdNode *n, int b) { 
    int id = Cudd_NodeReadIndex(n);
    int x = lookup_table_variables[id][0];
    int i = lookup_table_variables[id][1];
    // int t = lookup_table_variables[id][2];
    int j = lookup_table_variables[id][3];

    if (x == 0) {
        return get_prob_AO_encoded(hmm, i, j, b);
    } else if (x == 1) {
        return get_prob_AS1_encoded(hmm, j, b);
    } else if (x == 2) {
        return get_prob_AS_encoded(hmm, i, j, b);
    } else {
        perror("Error: Invalid value of x");
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
double Backward(DdManager* manager, DdNode* node, const HMM *hmm) {
    NodeDataNode *targetNodeData = FindTargetNodeAtLevel(manager, -1, node);
    
    if (targetNodeData->backward >= 0) {
        return targetNodeData->backward; 
    }


    DdNode* high = Cudd_T(node);
    DdNode* low = Cudd_E(node);
    
    double prob_high = get_prob_encoded(hmm, node, 1) * Backward(manager, high, hmm);
    double prob_low;

    if (!Cudd_IsComplement(low)) {
        prob_low = get_prob_encoded(hmm, node, 0) * Backward(manager, low, hmm);
    } else {
        // Adjusting for negative (complemented) edge
        prob_low = get_prob_encoded(hmm, node, 0) * (1.0 - Backward(manager, Cudd_Regular(low), hmm));
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
void CalculateForward(DdManager* manager, DdNode** F_seq, const HMM *hmm, int T) {
    int numVars = Cudd_ReadSize(manager); // Cudd_ReadNodeCount(manager)+(hmm->N)*T+2;
 
    for (int r = 0; r < hmm->N*T; r++) {
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
        double backwardVal = Backward(manager, targetNode, hmm);
        if (!isNegated) {
            // even number of comple edges
            targetNodeData->forward[1] += 1 / backwardVal;
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

            double ProbLowEdge = get_prob_encoded(hmm, targetNode->node, 0);
            double ProbHighEdge = get_prob_encoded(hmm, targetNode->node, 1);

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


void InitNodeData(DdManager* manager, DdNode** F_seq, const HMM *hmm, int T) {
    int numVars = Cudd_ReadSize(manager); // Cudd_ReadNodeCount(manager)+(hmm->N)*T+2;
    nodeData = (NodeDataList **)malloc((numVars + 1) * sizeof(NodeDataList*));

    if (nodeData == NULL) {
        // Handle memory allocation failure
        perror("Memory allocation failed");
        Cudd_Quit(manager);
    }
    for (int i = 0; i < numVars + 1; i++) {
        nodeData[i] = (NodeDataList *)malloc(sizeof(NodeDataList));
        nodeData[i]->head = NULL;
        nodeData[i]->tail = NULL;
    }

    DdNode *node;
    DdGen *gen;


    for (int r = 0; r < (hmm->N)*T; r++){
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
        for (int i0 = i+1; i0 < hmm->N; i0++){
            for (int j0 = 0; j0< hmm->M; j0++){
                sum += hmm->B[i0][j0];
            }
        }
        for (int j0 = j; j0 < hmm->M; j0++){
            sum += hmm->B[i][j0];
        }
        return sum;
    } else if (x == 1) {
        for (int j0 = j; j0 < hmm->N; j0++){
            sum += hmm->C[j0];
        }
        return sum;
    } else if (x == 2) {
        for (int i0 = i+1; i0 < hmm->N; i0++){
            for (int j0 = 0; j0 < hmm->N; j0++){
                sum += hmm->A[i0][j0];
            }
        }
        for (int j0 = j; j0 < hmm->N; j0++){
            sum += hmm->A[i][j0];
        }
        return sum;
    } else {
        perror("Error: Invalid value of x");
        return -1.0; 
    }
}


void computeConditionalExpectations(DdManager *manager, const HMM *hmm, int T, double ***eta) {

    int N = hmm->N;
    int NX = hmm->M; //(hmm->N > hmm->M) ? hmm->N : hmm->M;
    int numVars = Cudd_ReadSize(manager);
    // int T = 3;
    double e0, e1, temp;
    // *eta = (double***)malloc(3 * sizeof(double**));
    double ****etaX = (double****)malloc(3 * sizeof(double***)); // x is the leading dimension

    // Check memory allocation for etaX
    if (!etaX) { perror("Failed to allocate memory for etaX"); exit(EXIT_FAILURE); }

    double ****gamma = (double****)malloc(3 * sizeof(double***)); // Making gamma 3D with x as leading dimension

    // Check memory allocation for gamma
    if (!gamma) { perror("Failed to allocate memory for gamma"); exit(EXIT_FAILURE); }

    double *D = malloc((numVars+1) * sizeof(double));

    // Check memory allocation for D
    if (!D) { perror("Failed to allocate memory for D"); exit(EXIT_FAILURE); }

    for (int x = 0; x < 3; x++) {
        etaX[x] = (double***)malloc(N * sizeof(double**));
        if (!etaX[x]) { perror("Failed to allocate memory for etaX[x]"); exit(EXIT_FAILURE); }

        gamma[x] = (double***)malloc(N * sizeof(double**));
        if (!gamma[x]) { perror("Failed to allocate memory for gamma[x]"); exit(EXIT_FAILURE); }

        for (int i = 0; i < N; i++) {
            etaX[x][i] = (double**)malloc(T * sizeof(double*));
            if (!etaX[x][i]) { perror("Failed to allocate memory for etaX[x][i]"); exit(EXIT_FAILURE); }


            gamma[x][i] = (double**)malloc(T * sizeof(double*));
            if (!gamma[x][i]) { perror("Failed to allocate memory for gamma[x][i]"); exit(EXIT_FAILURE); }

            for (int t = 0; t < T; t++) {
                etaX[x][i][t] = (double*)malloc(NX * sizeof(double));
                if (!etaX[x][i][t]) { perror("Failed to allocate memory for etaX[x][i][t]"); exit(EXIT_FAILURE); }

                gamma[x][i][t] = (double*)malloc((NX+1) * sizeof(double));
                if (!gamma[x][i][t]) { perror("Failed to allocate memory for etaX[x][i][t]"); exit(EXIT_FAILURE); }
                

                for (int j = 0; j < NX; j++) {
                    etaX[x][i][t][j] = 0.0;
                    gamma[x][i][t][j] = 0.0;
                }
                gamma[x][i][t][NX] = 0.0;
            }
        }
    }

    for (int i = 0; i < numVars+1; i++) {
        D[i] = 0;
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
            double PrLow = get_prob_encoded(hmm, node->node, 0);
            double PrHigh = get_prob_encoded(hmm, node->node, 1);

            e1 = node->forward[0]*PrHigh*(1-highChildData->backward) + node->forward[1]*PrHigh*highChildData->backward;
            
            if (!Cudd_IsComplement(lowChild)) {
                e0 = node->forward[0]*PrLow*(1-lowChildData->backward) + node->forward[1]*PrLow*lowChildData->backward;
            } else {
                e0 = node->forward[0]*PrLow*(lowChildData->backward) + node->forward[1]*PrLow*(1-lowChildData->backward);
            }
            etaX[x][i][t][j] += e1;

            gamma[x][i][t][j+1] += e0;
            gamma[x][i][t][j] -= e0 + e1;

            D[level+1] += e0 + e1;
            D[highLevel] -= e1;
            D[lowLevel] -= e0;

            // next node
            node = node->next;
        }
    }

    temp = 0;
    for (int level = 1; level < numVars; level++) {
        int x = lookup_table_variables[level][0];
        int i = lookup_table_variables[level][1];
        int t = lookup_table_variables[level][2];
        int j = lookup_table_variables[level][3];
        temp += D[level];
        if (j == 0) {
            gamma[x][i][t][0] = temp;
        }
    }

    for (int x = 0; x < 3; x++) { 
        int tempN = N;
        if (x == 1) { tempN = 1; }
        for (int i = 0; i < tempN; i++) {
            for (int t = 0; t < T; t++) {
                temp = 0;
                if (x > 0) { NX = N; } 
                else { NX = hmm->M; }
                for (int j = 0; j < NX; j++) {
                    temp += gamma[x][i][t][j] / get_sigma(hmm, x, i, j); // TODO: / sigma[mu(i)][j];
                    etaX[x][i][t][j] += temp * get_theta(hmm, x, i, j);  // TODO: * thetda[mu(i)][j];
                }
            }
        }
    }

    for (int x = 0; x < 3; x++) { 
        int tempN = N;
        if (x == 1) { tempN = 1; }
        for (int i = 0; i < tempN; i++) {
            if (x > 0) { NX = N; } 
            else { NX = hmm->M; }
            for (int j = 0; j < NX; j++) {
                for (int t = 0; t < T; t++) {
                    eta[x][i][j]+= etaX[x][i][t][j];
                }
            }
        }
    }

    // Clean up
    for (int x = 0; x < 3; x++) {
        for (int i = 0; i < N; i++) {
            for (int t = 0; t < T; t++) {
                free(etaX[x][i][t]);
            }
            free(etaX[x][i]);
        }
        free(etaX[x]);
    }
    free(etaX);

    for (int x = 0; x < 3; x++) {
        for (int i = 0; i < N; i++) {
            for (int t = 0; t < T; t++) {
                free(gamma[x][i][t]);
            }
            free(gamma[x][i]);
        }
        free(gamma[x]);
    }
    free(gamma);

    free(D);
}

HMM* update(HMM *hmm, double ***eta) {
    HMM *new_hmm = HMM_create(hmm->N, hmm->M, "Updated model");
    double sum;
    // update b
    for (int u = 0; u < hmm->N; u++) {
        sum = 0.0;
        // ToDo: remove this quick fix, handles when obs is in alphabet, but not data to learn
        for (int o = 0; o < hmm->M; o++) {
            if (eta[0][u][o]>=0){
                sum += eta[0][u][o];
            }
        }
        for (int o = 0; o < hmm->M; o++) {
            // ToDo: remove this quick fix, handles when obs is in alphabet, but not data to learn
            if (eta[0][u][o]>=0){
                new_hmm->B[u][o] = eta[0][u][o] / sum;
            } else {
                new_hmm->B[u][o] = 0.0;
            }
        }
    }
    // update pi
    sum = 0.0;
    for (int v = 0; v < hmm->N; v++) {
        sum += eta[1][0][v];
    }
    for (int v = 0; v < hmm->N; v++) {
        new_hmm->C[v] = eta[1][0][v] / sum;
    }
    // update a
    for (int u = 0; u < hmm->N; u++) {
        sum = 0.0;
        for (int v = 0; v < hmm->N; v++) {
            sum += eta[2][u][v];
        }
        for (int v = 0; v < hmm->N; v++) {
            new_hmm->A[u][v] = eta[2][u][v] / sum;

        }
    }
    return new_hmm;
}


/**
 * @brief This function builds an HMM that "best" fits a given observation sequence O with the given number of states and observations.
 * 
 *      Uses the cour structure of the Baum-Welch algorithm, and BDDs, in order to improve the time and memory complexity
 *       
 * @param N       The number of variables.
 * @param M       The number of variables in AO.
 * @param T       The length of the sequence.
 * @param O       Observation sequence
 * @return struct HMM 
 */
HMM* learn(HMM *hypothesis_hmm, int T, int O[T])
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


    

    // Step 1 build (S)BDD
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0);  
    // Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);
    
    lookup_table_variables = (int **)malloc((N*T*M-T+N-1+(N*N*(T-1))) * sizeof(int *));
    for (int id = 0; id <(N*T*M-T+N-1+(N*N*(T-1))); id++) {
        (lookup_table_variables)[id] = (int *)malloc(4 * sizeof(int));
    }

    DdNode **F_obs = build_F_seq(manager, N, M, T, O);

    // Step 2: initilize M (= some random HMM) 
        // ToDo currently input


    double prob_priv, prob_new;
    int converged = 0;
    while (!converged)
    {
        prob_priv = log_likelihood_forward(model, O, T);

        // Step 3: E-step

        InitNodeData(manager, F_obs, model, T);
    
        // Step 3 (a) : Backward

        // Step 3 (b) : Forward
    
        CalculateForward(manager, F_obs, model, T);

        // Look at forward and backward values 
        // for (int level = Cudd_ReadSize(manager) ; level >= 0; level--) {
        //     NodeDataNode* current = nodeData[level]->head;
        //     while (current != NULL) {
        //         printf("Level %d: Node Address: %p \t %f\t %f\t %f\t %f\n", level, (void*)current->node,current->forward[0], current->forward[1], 1-current->backward, current->backward);
        //         current = current->next;
        //     }
        // }

        // Step 3 (c) : Conditional Expectations
        double ***eta;
        eta = (double***)malloc(3 * sizeof(double**));
        for (int i = 0; i < 3; ++i) {
            eta[i] = (double**)malloc(N * sizeof(double*));
            for (int j = 0; j < N; ++j) {
                eta[i][j] = (double*)malloc(M * sizeof(double));
                for (int k = 0; k < M; ++k) {
                    eta[i][j][k] = 0.0;
                }
            }
        }
        computeConditionalExpectations(manager, model, T, eta);
    
        // Step 4: M-step
        // Step 4 (a) : update M

        const HMM *new_hmm = update(model, eta);

        for (int x = 0; x < 3; x++) {
            for (int i = 0; i < N; i++) {
                free(eta[x][i]); // Free the innermost arrays
            }
            free(eta[x]); // Free the second level pointers
        }
        free(eta); // Finally, free the top level pointer

        double prob_new = log_likelihood_forward(new_hmm, O, T);
        
        
        HMM_copy(model, new_hmm); 
        HMM_destroy(new_hmm);

        printf("improvemment %f : %f, %f\n", prob_new-prob_priv, prob_new, prob_priv);
        // HMM_print(model);
        validate_hmm(model);
        if (prob_new <= prob_priv+0.05) {
            converged = 1;
        }

    }
    


    // ToDo: remove (For Debuging):
    // printf("N : %d | ", N );
    // printf("M : %d | ", M );
    // printf("T : %d | ", T );
    // printf("Encode: TRUE | ");
    // printf("DdManager vars: %d | ", Cudd_ReadSize(manager) ); /*Returns the number of BDD variables in existence*/
    // printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(manager) ); // countUniqueNodes(manager, pow(M,T), F_all) );/*Reports the number of live nodes in BDDs and ADDs*/
    // printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(manager) ); /*Returns the number of times reordering has occurred*/
    // printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager) ); /*Returns the memory in use by the manager measured in bytes*/


    // Step 6: Return the learned model

    return hypothesis_hmm;
}
