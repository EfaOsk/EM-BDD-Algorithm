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




    // DdNode **flattenedArray = (DdNode **)malloc((T + 1) * N * sizeof(DdNode *));

    // // Flatten the 2D array into the 1D array.
    // int count = 0;
    // for (int t = 0; t <= T; t++) {
    //     for (int n = 0; n < N; n++) {
    //         flattenedArray[count] = Cudd_BddToAdd(manager, FO[t][n]);
    //         count++;
    //     }
    // }
    // char filename[30];
    // sprintf(filename, "./graphs/graph_.dot");
    // FILE *outfile; 
    // outfile = fopen(filename,"w");
    // Cudd_DumpDot(manager, T*N, flattenedArray, NULL, NULL, outfile);
    // fclose(outfile);

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
void encode_variables(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M]) {
    // DdNode *AO_enc[N][T][M-1];
    // DdNode *AS1_enc[N - 1];
    // DdNode *AS_enc[N][T - 1][N-1];

    // Allocate memory for AO_enc
    DdNode ****AO_enc = (DdNode ****)malloc(N * sizeof(DdNode ***));
    for (int u = 0; u < N; u++) {
        AO_enc[u] = (DdNode ***)malloc(T * sizeof(DdNode **));
        for (int t = 0; t < T; t++) {
            AO_enc[u][t] = (DdNode **)malloc((M - 1) * sizeof(DdNode *));
        }
    }

    // Allocate memory for AS1_enc
    DdNode **AS1_enc = (DdNode **)malloc((N - 1) * sizeof(DdNode *));

    // Allocate memory for AS_enc
    DdNode ****AS_enc = (DdNode ****)malloc(N * sizeof(DdNode ***));
    for (int u = 0; u < N; u++) {
        AS_enc[u] = (DdNode ***)malloc(T * sizeof(DdNode **));
        for (int t = 0; t < T - 1; t++) {
            AS_enc[u][t] = (DdNode **)malloc((N - 1) * sizeof(DdNode *));
        }
    }

    int id = 0;

    // Encode AO
    for (int t = 0; t < T; t++) {
        for (int u = 0; u < N; u++) {
            for (int o = 0; o < M; o++) {
                if ((o == M-1) && (u == N-1)){

                } else {
                    printf("O^%d_%d = %d\n", u, t, o);
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
        printf("S_0 = %d\n", u);
        lookup_table_variables[id][0]= 1;
        lookup_table_variables[id][1]= u;
        lookup_table_variables[id][2]= -1;
        lookup_table_variables[id][3]= -1;
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
                    printf("S^%d_%d = %d\n", u, t+1, v);
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

    // // Free AO_enc
    // for (int u = 0; u < N; u++) {
    //     for (int t = 0; t < T; t++) {
    //         free(AO_enc[u][t]);
    //     }
    //     free(AO_enc[u]);
    // }
    // free(AO_enc);

    // // Free AS_enc
    // for (int u = 0; u < N; u++) {
    //     for (int t = 0; t < T - 1; t++) {
    //         for (int v = 0; v < N - 1; v++) {
    //             free(AS_enc[u][t][v]);
    //         }
    //         free(AS_enc[u][t]);
    //     }
    //     free(AS_enc[u]);
    // }
    // free(AS_enc);

    // // Free AS1_enc
    // for (int u = 0; u < N - 1; u++) {
    //     free(AS1_enc[u]);
    // }
    // free(AS1_enc);

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

void free_lookup_table_variables( int N, int M, int T) {
    for (int id = 0; id < (N*T*M-T+N-1+(N*N*(T-1))); id++) {
        free(lookup_table_variables[id]);
    }
    free(lookup_table_variables);
}

/**
 * Builds the FO nodes for all possible sequences.
 *
 * @param manager The Cudd Manager.
 * @param N       The number of variables.
 * @param M       The number of variables in AO.
 * @param T       The length of the sequence.
 *
 * @return An array of FO nodes for all possible sequences.
 */
DdNode **build_F_all_seq(DdManager *manager, int N, int M, int T) {
    
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
    


    // Calculate the total number of possible sequences
    int totalSequences = pow(M, T);

    // Create the array to hold the FO nodes for all possible sequences
    // DdNode* F_all[totalSequences];
    DdNode** F_all = malloc(totalSequences * sizeof(DdNode*));

    // Generate all possible sequences and build FO for each sequence
    for (int i = 0; i < totalSequences; i++) {
        int sequence[T];

        // Generate the sequence based on the current index 'i'
        for (int j = 0; j < T; j++) {
            int symbolIndex = (i / (int)pow(M, j)) % M;
            sequence[j] = symbolIndex;
        }

        // Build FO for the current sequence
        // TODO: Remove BddToAdd
        // F_all[i] = build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, sequence);


        DdNode* temp = build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, sequence);
        // F_all[i] =  Cudd_BddToAdd(manager, Cudd_bddAnd(manager, C_A, temp)); 
        F_all[i] =  Cudd_BddToAdd(manager, temp); 
        Cudd_Ref(F_all[i]);
        Cudd_RecursiveDeref(manager, temp);

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

    // Cudd_RecursiveDeref(manager, C_A);
    return F_all;
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
double get_prob_AO_encoded(const HMM *hmm, int i, int j, int b){

    if (b==0){ // false edge
        double sum_num = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++){
            for (int j0 = 0; j0<= hmm->M -1; j0++){
                sum_num += hmm->B[i0][j0];
            }
        }
        for (int j0 = j+1; j0<= hmm->M -1; j0++){
            sum_num += hmm->B[i][j0];
        }

        double sum_den = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++){
            for (int j0 = 0; j0<= hmm->M -1; j0++){
                sum_den += hmm->B[i0][j0];
            }
        }
        for (int j0 = j; j0<= hmm->M -1; j0++){
            sum_den += hmm->B[i][j0];
        }

        // 
        return sum_num / sum_den ;
    } else { // true edge
        double sum_den = 0.0;

        for (int i0 = i+1; i0 <= hmm->N - 1; i0++){
            for (int j0 = 0; j0<= hmm->M -1; j0++){
                sum_den += hmm->B[i0][j0];
            }
        }
        for (int j0 = j; j0<= hmm->M -1; j0++){
            sum_den += hmm->B[i][j0];
        }

        // a(i)(j) / 
        return (hmm->B[i][j]) / sum_den ;
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
        return get_prob_AS1_encoded(hmm, i, b);
    } else if (x == 2) {
        return get_prob_AS_encoded(hmm, i, j, b);
    } else {
        // error
    }
}


double** lookupTableBackward;

/**
 * @brief Backward prosedure where
 * 
 *      /Beta(x, n) = probability that paths logicallyreach the terminal node x from node n
 * 
 *  for level in vars: // from lowest level to highest
 *      for node n in level:
 *          B_1[n] = 0.5* B_1[Child(True)] + 0.5* B_1[Child(False)]
 * 
 * @param manager 
 * @param bdd 
 * @param M 
 * @return double
 */
double Backward(DdManager* manager, DdNode* bdd, const HMM *hmm, int b) {

    if (bdd == Cudd_ReadLogicZero(manager)) {
        // Terminal node: 0
        if (b==0){
            return 1.0;
        } else {
            return 0.0;
        }
    }
    if (bdd == Cudd_Not(Cudd_ReadLogicZero(manager))) {
        // Terminal node: 1
        if (b==0){
            return 0.0;
        } else {
            return 1.0;
        }
    }

    int id = Cudd_NodeReadIndex(bdd);

    // Check if the probability is already computed for this node
    if (lookupTableBackward[id][b] >= 0.0) {
        return lookupTableBackward[id][b];
    }



    // Compute probabilities for the true and false children
    DdNode* high = Cudd_E(bdd);
    DdNode* low = Cudd_T(bdd);

    // Calculate the probability for the current node
    double prob_high, prob_low, prob;
    prob_high =  get_prob_encoded(hmm, bdd, 1)*Backward(manager, high, hmm, b);
    prob_low = get_prob_encoded(hmm, bdd, 0)*Backward(manager, low, hmm, b);
    prob = prob_low + prob_high;
    lookupTableBackward[id][b] = prob; // Store in the lookup table


    return prob;
}

// double** lookupTableForward;

NodeWithForward *nodeData;

int findNodeIndex(DdNode *searchNode, int numNodes) {
    for (int i = 0; i < numNodes; i++) {
        if (nodeData[i].node == searchNode) {
            // Found the node in the list, return its index
            return i;
        }
    }
    // Node not found in the list
    return -1; // You can choose a suitable sentinel value for "not found"
}

/**
 * Calculate forward values for nodes in a Binary Decision Diagram (BDD).
 *
 * @param manager   The CUDD manager for the BDD operations.
 * @param F_all     An array of BDD nodes representing the roots of the BDDs.
 * @param hmm       A pointer to a Hidden Markov Model (HMM) or a similar probabilistic model.
 * @param num_roots The number of BDD roots in the 'F_all' array.
 */
void CalculateForward(DdManager* manager, DdNode** F_all, const HMM *hmm, int num_roots) {
    int numNodes = Cudd_ReadNodeCount(manager)+num_roots+2;
    nodeData = (NodeWithForward *)malloc(numNodes* sizeof(NodeWithForward));

    if (nodeData == NULL) {
        // Handle memory allocation failure
        perror("Memory allocation failed");
        Cudd_Quit(manager);
    }

    DdNode *node;
    DdGen *gen;
    int index = 0;

    int r= 1;
    for (int r = 0; r<num_roots; r++){
        Cudd_ForeachNode(manager, F_all[r], gen, node) {
            int i;
            for (i = 0; i < index; i++) {
                if (node == nodeData[i].node)
                    break;
            }
            // If the node is not in the set, add it
            if (i == index) {
                nodeData[index].node = node;
                nodeData[index].forward0 = 0.0; 
                nodeData[index].forward1 = 0.0; 
                // printf("\t%d", index);
                index++;
            }
        }
        nodeData[index].node = F_all[r];
        nodeData[index].forward1 = 1; //1/ Backward(manager, F_all[r], hmm, 0);
        nodeData[index].forward0 = 0;
        index++;

    }
    
    nodeData[index].node = Cudd_Not(Cudd_ReadLogicZero(manager));
    nodeData[index].forward0 = 0.0;
    nodeData[index].forward1 = 0.0;
    index++;
    nodeData[index].node = Cudd_ReadLogicZero(manager);
    nodeData[index].forward0 = 0.0;
    nodeData[index].forward1 = 0.0;
    index++;

    
    for (int i = 0; i < index; i++){
        node = nodeData[i].node;
        // cheack for terminal nodes
        if ((node !=  Cudd_Not(Cudd_ReadLogicZero(manager))) && (node !=  Cudd_ReadLogicZero(manager))){
            int id_high = findNodeIndex(Cudd_E(node), index+1);
            int id_low = findNodeIndex(Cudd_T(node), index+1);
            nodeData[id_low].forward0 += nodeData[i].forward0*get_prob_encoded(hmm, node, 0);
            nodeData[id_high].forward0 += nodeData[i].forward0*get_prob_encoded(hmm, node, 1);
            nodeData[id_low].forward1 += nodeData[i].forward1*get_prob_encoded(hmm, node, 0);
            nodeData[id_high].forward1 += nodeData[i].forward1*get_prob_encoded(hmm, node, 1);
        }
    }

    for (int i = 0; i < index; i++){
        printf("\t%f\t\t%f\n", nodeData[i].forward0, nodeData[i].forward1);
    }
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

    // TODO make the  HMM N and M an input ?
    int N = hypothesis_hmm->N;
    int M = hypothesis_hmm->M;

    lookup_table_variables = (int **)malloc((N*T*M-T+N-1+(N*N*(T-1))) * sizeof(int *));
    for (int id = 0; id <(N*T*M-T+N-1+(N*N*(T-1))); id++) {
        (lookup_table_variables)[id] = (int *)malloc(4 * sizeof(int));
    }

    // Step 1 build (S)BDD
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0);  
    // Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);
    

    DdNode **F_all = build_F_all_seq(manager, N, M, T);
    
        // test probabilty: ToDo: Remove :P
        // for (int u = 0; u < N - 1; u++) {
        //     printf("\t%f \n", get_prob_AS1_encoded(hypothesis_hmm, u, 1));
        //     printf("\t%f \n\n", get_prob_AS1_encoded(hypothesis_hmm, u, 0));
        // }
        // for (int u = 0; u < N; u++) {
        //     for (int t = 0; t < T-1; t++) {
        //         for (int v = 0; v < N; v++) {
        //             if ((v == N-1) && (u == N-1)){

        //             } else {
        //                 printf("\t%f \n", get_prob_AS_encoded(hypothesis_hmm, u, v, 1));
        //                 printf("\t%f \n\n", get_prob_AS_encoded(hypothesis_hmm, u, v, 0));
        //             }
        //         }
        //     }
        // }
        // for (int u = 0; u < N; u++) {
        //     for (int t = 0; t < T; t++) {
        //         for (int v = 0; v < M; v++) {
        //             if ((v == M-1) && (u == N-1)){

        //             } else {
        //                 printf("\t%f \n", get_prob_AO_encoded(hypothesis_hmm, u, v, 1));
        //                 printf("\t%f \n\n", get_prob_AO_encoded(hypothesis_hmm, u, v, 0));
        //             }
        //         }
        //     }
        // }

    // Step 2: initilize M (= some random HMM) 
        // ToDo currently input




    // Step 3: E-step
    
    // Step 3 (a) : Backward
        
    // Initialize the lookup table with -1.0 (indicating no value computed)

    // lookupTableBackward = (double **)malloc(Cudd_ReadSize(manager) * sizeof(double *));
    // for (int id = 0; id <Cudd_ReadSize(manager); id++) {
    //     (lookupTableBackward)[id] = (double *)malloc(2 * sizeof(double));
    // }


    // for (int i = 0; i < Cudd_ReadSize(manager); i++) {
    //     lookupTableBackward[i][0] = -1.0;
    //     lookupTableBackward[i][1] = -1.0;
    // }
    // double backward_probs = Backward(manager, F_all[0], hypothesis_hmm, 0);
    // printf("Backward Probabilities: %f\n", backward_probs);
    // backward_probs = Backward(manager, F_all[0], hypothesis_hmm, 1);
    // printf("Backward Probabilities: %f\n", backward_probs);

    // Step 3 (b) : Forward
    
    // CalculateForward(manager, F_all, hypothesis_hmm, (int) pow(M, T));

    // Step 3 (c) : Conditional Expectations


    // Step 4: M-step
    // Step 4 (a) : update M
        

    // Step 5: Calculate the log-likelyhood of M









    // ToDo: remove (For Debuging):
    printf("N : %d | ", N );
    printf("M : %d | ", M );
    printf("T : %d | ", T );
    printf("Encode: TRUE | ");
    printf("DdManager vars: %d | ", Cudd_ReadSize(manager) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(manager) ); // countUniqueNodes(manager, pow(M,T), F_all) );/*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(manager) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager) ); /*Returns the memory in use by the manager measured in bytes*/


    char filename[30];
    sprintf(filename, "./graphs/graph.dot"); /*Write .dot filename to a string*/

    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    printf("%f\n", pow(M, T));
    Cudd_DumpDot(manager, pow(M, T), F_all, NULL, NULL, outfile);
    fclose(outfile);


    Cudd_Quit(manager);

    // Step 6: Return the learned model

    return hypothesis_hmm;
}
