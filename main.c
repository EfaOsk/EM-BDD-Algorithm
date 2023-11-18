// #include <stdio.h>
// #include <stdlib.h>
// #include <stdbool.h>

#include "exampleHMM.h"
#include "HMM.h"

#include "learn.h"

int main(int argc, char *argv[]) {
    // Check for the correct number of arguments
    // if (argc != 4) {
    //     fprintf(stderr, "Usage: %s N M T\n where:\n\tN = number of states\n\tM = alphabet size\n\tT = length of observation sequence\n", argv[0]);
    //     return 1; // Return an error code
    // }

    // // Parse command-line arguments and validate them
    // int N, M, T;
    // if (sscanf(argv[1], "%d", &N) != 1 || sscanf(argv[2], "%d", &M) != 1 || sscanf(argv[3], "%d", &T) != 1) {
    //     fprintf(stderr, "Invalid arguments. N, M, and T must be integers.\n");
    //     return 1; // Return an error code
    // }

    // // The observation O
    // int O[T];
    // O[0] = 0;
    // O[1] = 1;

    HMM **example_models = initialize_example_models();
    
    // for (int i = 0; i < NUM_MODELS; i++) {
    //     if (example_models[i] == NULL) {
    //         fprintf(stderr, "Failed to create HMM model\n");
    //         return 1;
    //     }

    //     // Print the HMM model
    //     HMM_print(example_models[i] );
    //     validate_hmm(example_models[i]);
    //     // Destroy the HMM model to free memory
    //     HMM_destroy(example_models[i] );

    // }
    int T = 8;
    // int *observations = HMM_generate_sequence(example_models[7], T);

    // Calculate the probability of the sequence
    // double sequence_probability = probability_single_sequence(example_models[0], observations, 3);

    // Print the result
    // printf("Probability of observing the sequence: %f\n", sequence_probability);

    char result_folder[256];


    for (int i = 10; i < 11; i++) {
        int *observations = HMM_generate_sequence(example_models[i], T);
        for (int o = 0; o<T;o++) {
            printf("%d, ", observations[o]);
        }
        printf("\n");

        sprintf(result_folder, "experiments/experiments_0/model_%d/logs", i);
        HMM* hypothesis_hmm = HMM_random_create(example_models[i]->N, example_models[i]->M, "hypothesis model");
        HMM *learned_model = learn(hypothesis_hmm, T, observations, 0.001, result_folder, "experiments/experiments_0/results.txt"); 
        
        FILE *result_fp = fopen("experiments/experiments_0/results.txt", "a");
        fprintf(result_fp, ", %f\n", log_likelihood_forward(example_models[i], observations, T));
        fclose(result_fp);
    }

    free_example_models(example_models);
    return 0;
}
