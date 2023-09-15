#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
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
        FO[T][i] = Cudd_ReadOne(manager);
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
void encode_variables(DdManager *manager, int N, int M, int T, DdNode *AS1[N], DdNode *AS[N][T-1][N], DdNode *AO[N][T][M]) {
    DdNode *AO_enc[N][T][M-1];
    DdNode *AS1_enc[N - 1];
    DdNode *AS_enc[N][T - 1][N-1];


    // Encode AO
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M; o++) {
                if ((o == M-1) && (u == N-1)){

                } else {
                    AO_enc[u][t][o] = Cudd_bddNewVar(manager);
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
                    AO[v][t][o] = Cudd_ReadLogicZero(manager);
                    Cudd_Ref(AO[v][t][o]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int o0 = 0; o0 <= M-1; o0++){
                            DdNode *temp = Cudd_bddOr(manager, AO[v][t][o], Cudd_Not(AO_enc[v0][t][o0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AO[v][t][o]);
                            AO[v][t][o] = temp;
                        }
                    }
                    for (int o0 = 0; o0 < o; o0++){
                        DdNode *temp = Cudd_bddOr(manager, AO[v][t][o], Cudd_Not(AO_enc[v][t][o0]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AO[v][t][o]);
                        AO[v][t][o] = temp;
                    }
                    
                } else {

                    AO[v][t][o] = Cudd_ReadLogicZero(manager);
                    Cudd_Ref(AO[v][t][o]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int o0 = 0; o0 <= M-1; o0++){
                            DdNode *temp = Cudd_bddOr(manager, AO[v][t][o], Cudd_Not(AO_enc[v0][t][o0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AO[v][t][o]);
                            AO[v][t][o] = temp;
                        }
                    }
                    for (int o0 = 0; o0 < o; o0++){
                        DdNode *temp = Cudd_bddOr(manager, AO[v][t][o], Cudd_Not(AO_enc[v][t][o0]));
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
        AS1_enc[u] = Cudd_bddNewVar(manager);
    }

    for (int u = 0; u < N; u++) {
        if (u == 0) {
            AS1[0] = AS1_enc[0];
            Cudd_Ref(AS1[0]);
        } 
        else if (u < N-1) {
            AS1[u] = Cudd_ReadLogicZero(manager);
            Cudd_Ref(AS1[u]);
            for (int j = 0; j < u; j++){
                DdNode *temp = Cudd_bddOr(manager, AS1[u], Cudd_Not(AS1_enc[j]));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, AS1[u]);
                AS1[u] = temp;
            }
            DdNode *temp = Cudd_bddAnd(manager, AS1[u], AS1_enc[u]);
            Cudd_Ref(temp);
            Cudd_RecursiveDeref(manager, AS1[u]);
            AS1[u] = temp;
            
        } else {
            AS1[u] = Cudd_ReadLogicZero(manager);
            Cudd_Ref(AS1[u]);
            for (int j = 0; j < u; j++){
                DdNode *temp = Cudd_bddOr(manager, AS1[u], Cudd_Not(AS1_enc[j]));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(manager, AS1[u]);
                AS1[u] = temp;
            }
            
        }
    }


    // Encode AS
    // for (int u = 0; u < N; u++) {
    //     for (int t = 0; t < T-1; t++) {
    //         for (int v = 0; v < N - 1; v++) {
    //             AS_enc[u][t][v] = Cudd_bddNewVar(manager);
    //         }
    //     }
    // }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N; v++) {
                if ((v == N-1) && (u == N-1)){

                } else {
                    AS_enc[u][t][v] = Cudd_bddNewVar(manager);
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
                    AS[u][t][v] = Cudd_ReadLogicZero(manager);
                    Cudd_Ref(AS[u][t][v]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int u0 = 0; u0 <= N-1; u0++){
                            DdNode *temp = Cudd_bddOr(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AS[u][t][v]);
                            AS[u][t][v] = temp;
                        }
                    }
                    for (int u0 = 0; u0 < u; u0++){
                        DdNode *temp = Cudd_bddOr(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v]));
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(manager, AS[u][t][v]);
                        AS[u][t][v] = temp;
                    }
                    
                } else {

                    AS[u][t][v] = Cudd_ReadLogicZero(manager);
                    Cudd_Ref(AS[u][t][v]);
                    for (int v0 = 0; v0 <= v-1; v0++){
                        for (int u0 = 0; u0 <= N-1; u0++){
                            DdNode *temp = Cudd_bddOr(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v0]));
                            Cudd_Ref(temp);
                            Cudd_RecursiveDeref(manager, AS[u][t][v]);
                            AS[u][t][v] = temp;
                        }
                    }
                    for (int u0 = 0; u0 < u; u0++){
                        DdNode *temp = Cudd_bddOr(manager, AS[u][t][v], Cudd_Not(AS_enc[u0][t][v]));
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
    
    
    DdNode* C_A= build_C_A(manager, N, M, T, AS1, AS, AO);
    


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
        F_all[i] =  Cudd_bddAnd(manager, C_A, temp); 
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

    Cudd_RecursiveDeref(manager, C_A);
    return F_all;
}


/**
 * @brief Backward prosedure where
 * 
 *      /Beta(x, n) = probability that paths logicallyreach the terminal node x from node n
 * 
 * @param manager 
 * @param bdd 
 * @param M 
 * @return struct DoubleArray 
 */
struct DoubleArray Backward(DdManager *manager, DdNode *bdd, struct HMM M)
{
    // initilize beta (B)
    
    /*
    for level in vars:
        for node n in level:
            for x in {0, 1}
                B(x)(n) = ???
    */

    // DdGen *tmp = Cudd_FirstNode(manager, )
    return (struct DoubleArray){};
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
struct HMM learn(const int N, const int M, int T, int O[T])
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

    // Step 1 build (S)BDD
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0);  
    // Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);
    

    DdNode **F_all = build_F_all_seq(manager, N, M, T);

    // ChatGPT, what sould I put as the number of minterms here?

    // printf("DdNodes: %ld \n", Cudd_ReadNodeCount(manager) );

    // printf("DdManager nodes: %d | ", Cudd_DagSize(F_all)); /*Reports the number of live nodes in BDDs and ADDs*/
    printf("N : %d | ", N );
    printf("M : %d | ", M );
    printf("T : %d | ", T );
    printf("Encode: TRUE | ");
    printf("DdManager vars: %d | ", Cudd_ReadSize(manager) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(manager) ); // countUniqueNodes(manager, pow(M,T), F_all) );/*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(manager) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager) ); /*Returns the memory in use by the manager measured in bytes*/


    char filename[30];
    sprintf(filename, "./graph/graph.dot"); /*Write .dot filename to a string*/

    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    Cudd_DumpDot(manager, pow(M, T), F_all, NULL, NULL, outfile);
    fclose(outfile);

    // bdd = Cudd_BddToAdd(manager, bdd); 
    // write_dd(manager, bdd, filename);
    Cudd_Quit(manager);

    return (struct HMM){};
}