#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <time.h>
#include "BDD.h"
#include "cudd.h"

int main(int argc, char *argv[]) {
   
    int T = 2;
    int num_sequences = 4;
    double epsilon =  0.01;
    char logs_folder[256]; 
    int **obs_seq;
    int N = 2;
    int M = 2;
    obs_seq = malloc(num_sequences * sizeof(int*));

    
    for (int j = 0; j < num_sequences; j++) {
        obs_seq[j] = generate_sequence(M, T);
    } 
    obs_seq[0][0] = 0;
    obs_seq[0][1] = 3;
    obs_seq[0][2] = 2;
    obs_seq[0][3] = 3;
    obs_seq[0][4] = 2;

    printf("%d, ", obs_seq[0][0]);
    printf("%d, ", obs_seq[0][1]);
    printf("%d, ", obs_seq[0][2]);
    // printf("%d, ", obs_seq[0][3]);
    // printf("%d\n", obs_seq[0][4]);
    
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0);  
    Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);
    
    int **lookup_table_variables;
    lookup_table_variables = (int **)malloc((N*(M-1)+(N-1)+(N*(N-1))) * sizeof(int *));
    for (int id = 0; id <(N*(M-1)+(N-1)+(N*(N-1))); id++) {
        (lookup_table_variables)[id] = (int *)malloc(4 * sizeof(int));
    }

    printf("Building the BDD ...\n");
    DdNode **F_obs = build_F_seq(manager, N, M, num_sequences, T, obs_seq, lookup_table_variables);


    printf("N : %d | ", N );
    printf("M : %d | ", M );
    printf("T : %d | ", T );
    printf("Encode: TRUE | ");
    printf("DdManager vars: %d | ", Cudd_ReadSize(manager) ); // Returns the number of BDD variables in existence
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(manager) ); // Reports the number of live nodes in BDDs and ADDs
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(manager) ); // Returns the number of times reordering has occurred
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager) ); // Returns the memory in use by the manager measured in bytes

    Cudd_PrintDebug(manager, F_obs[0], 0, 3);
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen("test.dot","w");
    // DdNode **ddnodearray = (DdNode**)malloc(sizeof(DdNode*)); // initialize the function array
    // ddnodearray[0] = dd;
    Cudd_DumpDot(manager, 4, F_obs, NULL, NULL, outfile);
    for (int i = 0; i < num_sequences; i++) {
        free(obs_seq[i]);
    }
    // Free the allocated memory after use
    free(obs_seq);
    return 0;
}
