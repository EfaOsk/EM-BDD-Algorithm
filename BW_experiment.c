#include "exampleHMM.h"
#include "HMM.h"
#include <math.h>
#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>


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