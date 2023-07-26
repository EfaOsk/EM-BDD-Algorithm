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
                Cudd_RecursiveDeref(manager, disjunction);
                Cudd_Ref(FO[t][i]);
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
        Cudd_Ref(FO[0][i]);

        Cudd_RecursiveDeref(manager, temp);
        Cudd_RecursiveDeref(manager, recursive);
    }

    DdNode *FO_ = Cudd_ReadLogicZero(manager);
    Cudd_Ref(FO_);
    for (int i = 0; i < N; i++) {
        FO_ = Cudd_bddOr(manager, FO_, FO[0][i]);
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
    DdNode *AO_enc[N][T][M - 1];
    DdNode *AS1_enc[N - 1];
    DdNode *AS_enc[N][T - 1][N - 1];


    // Encode AO
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M - 1; o++) {
                AO_enc[u][t][o] = Cudd_bddNewVar(manager);
                Cudd_Ref(AO_enc[u][t][o]);
            }
        }
    }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M; o++) {
                if (o == 0) {
                    AO[u][t][0] = AO_enc[u][t][0];
                    Cudd_Ref(AO[u][t][0]);
                } 
                else if (o < M-1)
                {
                    AO[u][t][o] = Cudd_ReadLogicZero(manager);
                    for (int j = 0; j < o; j++){
                        AO[u][t][o] = Cudd_bddOr(manager, AO[u][t][o], Cudd_Not(AO_enc[u][t][j]));
                    }
                    AO[u][t][o] = Cudd_bddAnd(manager, AO[u][t][o], AO_enc[u][t][o]);
                    Cudd_Ref(AO[u][t][o]);
                } else {
                    AO[u][t][o] = Cudd_ReadLogicZero(manager);
                    for (int j = 0; j < o; j++){
                        AO[u][t][o] = Cudd_bddOr(manager, AO[u][t][o], Cudd_Not(AO_enc[u][t][j]));
                    }
                    Cudd_Ref(AO[u][t][o]);
                }
            }
        }
    }

    // for (int u = 0; u < N; u++) {
    //     for (int t = 0; t < T; t++) {
    //         for (int o = 0; o < M; o++) {
    //             Cudd_RecursiveDeref(manager, AO[u][t][o]);
    //         }
    //     }
    // }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T; t++) {
            for (int o = 0; o < M-1; o++) {
                Cudd_RecursiveDeref(manager, AO_enc[u][t][o]);
            }
        }
    }
    // Encode AS1
    for (int u = 0; u < N - 1; u++) {
        AS1_enc[u] = Cudd_bddNewVar(manager);
        Cudd_Ref(AS1_enc[u]);
    }

    for (int u = 0; u < N; u++) {
        if (u == 0) {
            AS1[0] = AS1_enc[0];
            Cudd_Ref(AS1[0]);
        } 
        else if (u < N-1) {
            AS1[u] = Cudd_ReadLogicZero(manager);
            for (int j = 0; j < u; j++){
                AS1[u] = Cudd_bddOr(manager, AS1[u], Cudd_Not(AS1_enc[j]));
            }
            AS1[u] = Cudd_bddAnd(manager, AS1[u], AS1_enc[u]);
            Cudd_Ref(AS1[u]);
        } else {
            AS1[u] = Cudd_ReadLogicZero(manager);
            for (int j = 0; j < u; j++){
                AS1[u] = Cudd_bddOr(manager, AS1[u], Cudd_Not(AS1_enc[j]));
            }
            Cudd_Ref(AS1[u]);
        }
    }
    // for (int u = 0; u < N; u++) {
    //     Cudd_RecursiveDeref(manager, AS1[u]);
    // }
    for (int u = 0; u < N-1; u++) {
        Cudd_RecursiveDeref(manager, AS1_enc[u]);
    }

    // Encode AS
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N - 1; v++) {
                AS_enc[u][t][v] = Cudd_bddNewVar(manager);
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
                else if (v < N-1)
                {
                    AS[u][t][v] = Cudd_ReadLogicZero(manager);
                    for (int j = 0; j < v; j++){
                        AS[u][t][v] = Cudd_bddOr(manager, AS[u][t][v], Cudd_Not(AS_enc[u][t][j]));
                    }
                    AS[u][t][v] = Cudd_bddAnd(manager, AS[u][t][v], AS_enc[u][t][v]);
                    Cudd_Ref(AS[u][t][v]);
                } else {
                    AS[u][t][v] = Cudd_ReadLogicZero(manager);
                    for (int j = 0; j < v; j++){
                        AS[u][t][v] = Cudd_bddOr(manager, AS[u][t][v], Cudd_Not(AS_enc[u][t][j]));
                    }
                    Cudd_Ref(AS[u][t][v]);
                }
            }
        }
    }

    // for (int u = 0; u < N; u++) {
    //     for (int t = 0; t < T-1; t++) {
    //         for (int v = 0; v < N; v++) {
    //             Cudd_RecursiveDeref(manager, AS[u][t][v]);
    //         }
    //     }
    // }
    for (int u = 0; u < N; u++) {
        for (int t = 0; t < T-1; t++) {
            for (int v = 0; v < N-1; v++) {
                Cudd_RecursiveDeref(manager, AS_enc[u][t][v]);
            }
        }
    }
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


    // encode_variables(manager, N, M, T, AS1, AS, AO);

    // If not encode, set the initial varables
        for (int u = 0; u < N; u++){
            for (int t = 0; t < T; t++){
                for (int o = 0; o < M; o++){
                    AO[u][t][o] = Cudd_bddNewVar(manager);
                }
            }
        }

        for (int u = 0; u < N; u++){
            AS1[u] = Cudd_bddNewVar(manager);
        }

        for (int u = 0; u < N; u++){
            for (int t = 0; t < T-1; t++){
                for (int v = 0; v < N; v++){
                    AS[u][t][v] = Cudd_bddNewVar(manager);
                }
            }
        }
    
    // encode_variables(manager, N, M, T, AS1, AS, AO);


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
        F_all[i] = Cudd_BddToAdd(manager, build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, sequence));
        // F_all[i] = build_F_single_seq_O(manager, N, M, T, AS1, AS, AO, sequence);
        // Cudd_RecursiveDeref(manager, F_all[i]);

    }


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

unsigned int countUniqueNodes(DdManager *manager, int n, DdNode **bdds) {
    DdGen *gen;
    DdNode *node;
    int count = 0;

    // Create a set to store unique node addresses
    int initialNodeCount = Cudd_ReadNodeCount(manager);
    int maxNodeCount = initialNodeCount * 2; // To avoid resizing
    void **nodeSet = (void **)malloc(maxNodeCount * sizeof(void *));
    assert(nodeSet != NULL);

    for (int j = 0; j<n; j++){
        // Iterate over nodes and count unique nodes
        Cudd_ForeachNode(manager, bdds[j], gen, node) {
            // Check if the node address is already in the set
            int i;
            for (i = 0; i < count; i++) {
                if (node == nodeSet[i])
                    break;
            }

            // If the node is not in the set, add it and increment count
            if (i == count) {
                nodeSet[count] = node;
                count++;
            }
        }
    }

    // Free the memory used by the node set
    free(nodeSet);

    return count;
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
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 100000);    
    DdNode **F_all = build_F_all_seq(manager, N, M, T);

    // ChatGPT, what sould I put as the number of minterms here?

    // printf("DdManager nodes: %d | ", Cudd_DagSize(F_all)); /*Reports the number of live nodes in BDDs and ADDs*/
    printf("N = %d | ", N );
    printf("M = %d | ", M );
    printf("Encode: TRUE | ");
    printf("DdManager vars: %d | ", Cudd_ReadSize(manager) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager nodes: %d | ", countUniqueNodes(manager, pow(M,T), F_all)); // Cudd_ReadNodeCount(manager) );/*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(manager) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager) ); /*Returns the memory in use by the manager measured in bytes*/
    // Cudd_PrintDebug(manager, F_all[0], 2, 4);

    char filename[30];
    sprintf(filename, "./graph/graph.dot"); /*Write .dot filename to a string*/

    sprintf(filename, "graphs/test_.dot"); /*Write .dot filename to a string*/
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    Cudd_DumpDot(manager, pow(M, T), F_all, NULL, NULL, outfile);
    fclose(outfile);

    // bdd = Cudd_BddToAdd(manager, bdd); 
    // write_dd(manager, bdd, filename);
    Cudd_Quit(manager);

    return (struct HMM){};
}
