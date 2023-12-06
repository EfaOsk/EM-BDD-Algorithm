// #include <stdio.h>
// #include <stdlib.h>
// #include <stdbool.h>

#include "exampleHMM.h"
#include "HMM.h"
#include <math.h>
#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>

#include "learn.h"

/*
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
    //     HMM_validate(example_models[i]);
    //     // Destroy the HMM model to free memory
    //     HMM_destroy(example_models[i] );

    // }
    int T = 4;
    // int *observations = HMM_generate_sequence(example_models[7], T);

    // Calculate the probability of the sequence
    // double sequence_probability = probability_single_sequence(example_models[0], observations, 3);

    // Print the result
    // printf("Probability of observing the sequence: %f\n", sequence_probability);

    char result_folder[256];


    for (int i = 0; i < 11; i++) {
        int *observations = HMM_generate_sequence(example_models[i], T);
        for (int o = 0; o<T;o++) {
            printf("%d, ", observations[o]);
        }
        printf("\n");

        sprintf(result_folder, "experiments/experiments_test/model_%d/logs", i);
        HMM* hypothesis_hmm = HMM_random_create(example_models[i]->N, example_models[i]->M, "hypothesis model");
        HMM_validate(hypothesis_hmm);
        HMM *learned_model = learn(hypothesis_hmm, T, observations, 0.0001, result_folder, "experiments/experiments_0/results.txt"); 
        
        // FILE *result_fp = fopen("experiments/experiments_0/results.txt", "a");
        // fprintf(result_fp, ", %f\n", log_likelihood_forward(example_models[i], observations, T));
        // fclose(result_fp);
    }

    free_example_models(example_models);
    return 0;
}
*/



int main(int argc, char *argv[]) {

    HMM **example_models = initialize_example_models();
    
    int T = 6;
    int NO = 10;
    int *test;
    char experiment_name[19];
    char obs_seq_file[256];
    char hypo_hmm_file[256];
    char result_folder[256];
    char meta_data_file[256];
    sprintf(experiment_name, "experiments_3_BDD");

    sprintf(meta_data_file, "experiments/%s", experiment_name);
    mkdir(meta_data_file, 0777);
    sprintf(meta_data_file, "experiments/%s/experement_info.txt", experiment_name);
    FILE *metadata = fopen(meta_data_file, "a");
    fprintf(metadata, "%s\n", experiment_name);
    fprintf(metadata, "Original models: example models 0-10\n");
    fprintf(metadata, "Hypothesis models: data/obs_seq_Example Model 0-10_T6/\n");
    fprintf(metadata, "Observations: random from original models\n");
    fprintf(metadata, "epsilon: %f\n", 0.001);
    fprintf(metadata, "date: %d-%02d-%02d\n", 2023, 11, 29);
    fclose(metadata);


    int **obs_seq;
    obs_seq = malloc(NO * sizeof(int*));

    for (int m = 0; m < 11; m++) {
        // sprintf(obs_seq_file, "data/obs_seq_Example Model 0-10_T6/obs_seq_%s_T%d.txt", example_models[m]->name, T);
        // int **obs_seq = read_list(obs_seq_file);
        // sprintf(hypo_hmm_file, "data/obs_seq_Example Model 0-10_T6/obs_seq_%s_T%d.txt", example_models[m]->name, T);
        // HMM **hypo_models = HMM_load(hypo_hmm_file);
        for (int obs = 0; obs < 10; obs++) {
            for (int j = 0; j < NO; j++) {
                obs_seq[j] = HMM_generate_sequence(example_models[m], T);
            } 
            
            for (int i = 0; i < 10; i++) {
                sprintf(result_folder, "experiments/%s/%s_Obs%d_Hypo%d", experiment_name, example_models[m]->name, obs, i);
                mkdir(result_folder, 0777);
                sprintf(result_folder, "experiments/%s/%s_Obs%d_Hypo%d/logs", experiment_name, example_models[m]->name, obs, i);
                mkdir(result_folder, 0777);

                FILE *result_fp = fopen("experiments/experiments_3_BDD/results.txt", "a");
                fprintf(result_fp, "Model: %s, Obs: %d, Run: %d, L(Org): %f, ", example_models[m]->name, obs, i, log_likelihood_forward_multiple(example_models[m], obs_seq, NO, T));
                fclose(result_fp);
                HMM* hypothesis_hmm = HMM_random_create(example_models[m]->N, example_models[m]->M, "hypothesis model");
                HMM *learned_model = learn(hypothesis_hmm, T, NO, obs_seq, 0.001, result_folder, "experiments/experiments_3_BDD/results.txt");
            }
        }
    }


    // Free the allocated memory after use
    for (int i = 0; i < NO; i++) {
        free(obs_seq[i]);
    }
    free_example_models(example_models);
    return 0;
}

