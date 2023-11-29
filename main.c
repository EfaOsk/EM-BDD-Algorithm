// #include <stdio.h>
// #include <stdlib.h>
// #include <stdbool.h>

#include "exampleHMM.h"
#include "HMM.h"
#include <math.h>
#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>

// #include "learn.h"
/**
int main_BDD(int argc, char *argv[]) {
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
    int T = 3;
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

        sprintf(result_folder, "experiments/experiments_0/model_%d/logs", i);
        HMM* hypothesis_hmm = HMM_random_create(example_models[i]->N, example_models[i]->M, "hypothesis model");
        HMM_validate(hypothesis_hmm);
        HMM *learned_model = learn(hypothesis_hmm, T, observations, 0.001, result_folder, "experiments/experiments_0/results.txt"); 
        
        FILE *result_fp = fopen("experiments/experiments_0/results.txt", "a");
        fprintf(result_fp, ", %f\n", log_likelihood_forward(example_models[i], observations, T));
        fclose(result_fp);
    }

    free_example_models(example_models);
    return 0;
}

*/

// int main(int argc, char *argv[]) {

//     HMM **example_models = initialize_example_models();
    
//     int T = 6;
//     char model_folder[256];
//     char result_folder[256];
//     char log_filename[256];
//     char obs_seq_file[256];
//     char experiment_name[256];
//     sprintf(experiment_name, "experiments_0_BW");

//     sprintf(file_name, "experiments/%s/experement_info.txt", experiment_name);
//     FILE *metadata = fopen(filename, "a");
    
//     fprintf(metadata, "%s\n", experiment_name);
//     fprintf(metadata, "Original models: example models 0-10\n");
//     fprintf(metadata, "Hypothesis models: data/random_models\n");
//     fprintf(metadata, "Observations: data/obs_seq_Example Model 0-10_T6\n");
//     fprintf(metadata, "epsilon: %f\n", 0.0001);
//     fprintf(metadata, "date: %d-%02d-%02d\n", 2023, 11, 28);
//     fclose(metadata);

//     for (int m = 0; m < 11; m++) {
//         sprintf(model_folder, "experiments/%s/%s", experiment_name, example_models[m]->name);
//         mkdir(model_folder)
//         sprintf(obs_seq_file, "data/obs_seq_Example Model 0-10_T6/obs_seq_%s_T%d.txt", example_models[m]->name, T);
//         int **obs_seq = read_list(obs_seq_file);
//         for (int obs = 0; obs < 10; obs++) {
//             for (int i = 0; i < 10; i++) {
//                 sprintf(result_folder, "%s/logs", model_folder);
//                 sprintf(log_filename, "%s/log.txt", result_folder);
//                 FILE *log_file = fopen(log_filename, "a");
//                 fprintf("%s, obs %d, run %d\n");
//                 fclose(log_file);
//                 HMM *learned_model = HMM_learn(hypothesis_hmm, T, observations, 0.0001, result_folder, "experiments/experiments_0_BW/results.txt");
//             }
//         }

//         sprintf(result_folder, "experiments/experiments_0/model_%d/logs", i);
//         HMM* hypothesis_hmm = HMM_random_create(example_models[i]->N, example_models[i]->M, "hypothesis model");
//         HMM_validate(hypothesis_hmm);
//         HMM *learned_model = HMM_learn(hypothesis_hmm, T, observations, 0.0001, result_folder, "experiments/experiments_0/results.txt"); 
        
//         // FILE *result_fp = fopen("experiments/experiments_0/results.txt", "a");
//         // fprintf(result_fp, ", %f\n", log_likelihood_forward(example_models[i], observations, T));
//         // fclose(result_fp);

//         printf("\n");
//     }

//     free_example_models(example_models);
//     return 0;
// }


int main(int argc, char *argv[]) {

    HMM **example_models = initialize_example_models();
    
    int T = 6;
    int *test;
    char experiment_name[18];
    char obs_seq_file[256];
    char hypo_hmm_file[256];
    char result_folder[256];
    char meta_data_file[256];
    sprintf(experiment_name, "experiments_0_BW");

    sprintf(meta_data_file, "experiments/%s", experiment_name);
    mkdir(meta_data_file, 0777);
    sprintf(meta_data_file, "experiments/%s/experement_info.txt", experiment_name);
    FILE *metadata = fopen(meta_data_file, "a");
    fprintf(metadata, "%s\n", experiment_name);
    fprintf(metadata, "Original models: example models 0-10\n");
    fprintf(metadata, "Hypothesis models: data/random_models\n");
    fprintf(metadata, "Observations: data/obs_seq_Example Model 0-10_T6\n");
    fprintf(metadata, "epsilon: %f\n", 0.0001);
    fprintf(metadata, "date: %d-%02d-%02d\n", 2023, 11, 29);
    fclose(metadata);

    for (int m = 0; m < 11; m++) {
        sprintf(obs_seq_file, "data/obs_seq_Example Model 0-10_T6/obs_seq_%s_T%d.txt", example_models[m]->name, T);
        int **obs_seq = read_list(obs_seq_file);
        sprintf(hypo_hmm_file, "data/obs_seq_Example Model 0-10_T6/obs_seq_%s_T%d.txt", example_models[m]->name, T);
        // HMM **hypo_models = HMM_load(hypo_hmm_file);
        for (int obs = 0; obs < 10; obs++) {
            for (int i = 0; i < 10; i++) {
                sprintf(result_folder, "experiments/%s/%s_Obs%d_Hypo%d", experiment_name, example_models[m]->name, obs, i);
                mkdir(result_folder, 0777);
                sprintf(result_folder, "experiments/%s/%s_Obs%d_Hypo%d/logs", experiment_name, example_models[m]->name, obs, i);
                mkdir(result_folder, 0777);

                FILE *result_fp = fopen("experiments/experiments_0_BW/results.txt", "a");
                fprintf(result_fp, "Model: %s, Obs: %d, Run: %d, L(Org): %f, ", example_models[m]->name, obs, i, log_likelihood_forward(example_models[m], obs_seq[obs], T));
                fclose(result_fp);
                HMM* hypothesis_hmm = HMM_random_create(example_models[m]->N, example_models[m]->M, "hypothesis model");
                HMM *learned_model = HMM_learn(hypothesis_hmm, T, obs_seq[obs], 0.0001, result_folder, "experiments/experiments_0_BW/results.txt");
            }
        }
    }

    free_example_models(example_models);
    return 0;
}

