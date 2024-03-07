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
#include <time.h>
#include "BDD.h"



int main(int argc, char *argv[]) {

    HMM **example_models = initialize_example_models_large();
    
    int T = 50;
    int num_sequences = 10000;
    double epsilon = 0.01;
    HMM *hypothesis_hmm, *learned_model_BW, *learned_model_EMBDD;

    int **obs_seq;
    obs_seq = malloc(num_sequences * sizeof(int*));

    for (int m = 0; m < 10; m++) {
        for (int obs = 0; obs < 1; obs++) {
            for (int j = 0; j < num_sequences; j++) {
                obs_seq[j] = HMM_generate_sequence(example_models[m], T);
            } 
            hypothesis_hmm = HMM_random_create(example_models[m]->N, example_models[m]->M, "hypothesis model");
            printf("Baum-Welch %s\n", example_models[m]->name);
            char bw_filename[256];
            sprintf(bw_filename, "%s_BW_learned_model.txt", example_models[m]->modelName);
            HMM_save(learned_model_BW, bw_filename);
            learned_model_BW = BW_learn_multiple(hypothesis_hmm, num_sequences, obs_seq, T, epsilon);
            printf("EM-BDD %s\n", example_models[m]->name);
            learned_model_EMBDD = EMBDD_learn(hypothesis_hmm, num_sequences, obs_seq, T, epsilon);
            char embdd_filename[256];
            sprintf(embdd_filename, "%s_EMBDD_learned_model.txt", example_models[m]->modelName);
            HMM_save(learned_model_EMBDD, embdd_filename);
        }
    }


    // Free the allocated memory after use
    for (int i = 0; i < num_sequences; i++) {
        free(obs_seq[i]);
    }
    free_example_models(example_models);
    return 0;
}

