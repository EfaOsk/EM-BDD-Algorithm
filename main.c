// #include <stdio.h>
// #include <stdlib.h>
// #include <stdbool.h>

#include "exampleHMM.h"
#include "HMM.h"

#include "learn.h"

int main(int argc, char *argv[]) {
    // Check for the correct number of arguments
    if (argc != 4) {
        fprintf(stderr, "Usage: %s N M T\n where:\n\tN = number of states\n\tM = alphabet size\n\tT = length of observation sequence\n", argv[0]);
        return 1; // Return an error code
    }

    // Parse command-line arguments and validate them
    int N, M, T;
    if (sscanf(argv[1], "%d", &N) != 1 || sscanf(argv[2], "%d", &M) != 1 || sscanf(argv[3], "%d", &T) != 1) {
        fprintf(stderr, "Invalid arguments. N, M, and T must be integers.\n");
        return 1; // Return an error code
    }

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
    int observations[] = {0, 1, 2};

    // Calculate the probability of the sequence
    // double sequence_probability = probability_single_sequence(example_models[0], observations, 3);

    // Print the result
    // printf("Probability of observing the sequence: %f\n", sequence_probability);


    HMM *learned_model = learn(example_models[1], 1, NULL); 
    return 0;
}
