#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "BDD.h"
#include "cudd.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <T>\n", argv[0]);
        return 1;
    }

    int T = 6;  // Get T from command line arguments
    int num_sequences = 100;
    double epsilon = 0.01;
    char logs_folder[256];
    int **obs_seq;
    int N = atoi(argv[1]);
    int M = 4;
    obs_seq = malloc(num_sequences * sizeof(int*));

    for (int j = 0; j < num_sequences; j++) {
        obs_seq[j] = generate_sequence(M, T);
    }

    DdManager *manager = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);

    int **lookup_table_variables;
    lookup_table_variables = (int **)malloc((N*T*(M-1)+(N-1)+T*(N*(N-1))) * sizeof(int *));
    for (int id = 0; id < (N*T*(M-1)+(N-1)+T*(N*(N-1))); id++) {
        lookup_table_variables[id] = (int *)malloc(4 * sizeof(int));
    }

    // printf("Building the BDD ...\n");
    DdNode **F_obs = build_F_seq(manager, N, M, num_sequences, T, obs_seq, lookup_table_variables);

    printf("N : %d | ", N);
    printf("M : %d | ", M);
    printf("T : %d | ", T);
    printf("Encode: TRUE | ");
    printf("DdManager vars: %d | ", Cudd_ReadSize(manager));
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(manager));
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(manager));
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(manager));

    // Cudd_PrintDebug(manager, F_obs[0], 0, 3);
    // FILE *outfile;
    // outfile = fopen("test.dot", "w");
    // Cudd_DumpDot(manager, 4, F_obs, NULL, NULL, outfile);

    for (int i = 0; i < num_sequences; i++) {
        free(obs_seq[i]);
    }
    free(obs_seq);
    return 0;
}
