#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "HMM.h"
#include "helpers.h"
#include "BDD.h"

// Function to evaluate models and log results
void evaluate_model(HMM* model, HMM* original_model, int** sequences, int num_sequences, int sequence_length, const char* log_filename) {
    double log_likelihood_learned = log_likelihood_forward_multiple(model, sequences, num_sequences, sequence_length);
    double log_likelihood_original = log_likelihood_forward_multiple(original_model, sequences, num_sequences, sequence_length);

    FILE* file = fopen(log_filename, "w");
    if (file) {
        fprintf(file, "Log-Likelihood of Learned Model: %f\n", log_likelihood_learned);
        fprintf(file, "Log-Likelihood of Original Model: %f\n", log_likelihood_original);
        fprintf(file, "Improvement: %f\n", log_likelihood_learned - log_likelihood_original);
        fclose(file);
    } else {
        perror("Error writing log file");
    }
}

int main(int argc, char* argv[])  {
    // Check if all necessary arguments are provided
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <dataset_file> <num_states> <epsilon>\n", argv[0]);
        return 1;
    }

    // Get parameters from command-line arguments
    const char* dataset_file = argv[1];
    int num_states = atoi(argv[2]);
    double epsilon = atof(argv[3]);

    int num_observations;
    int sequence_length;
    int num_sequences;

    // Load dataset
    int** sequences = load_dataset(dataset_file, &num_observations, &num_sequences, &sequence_length);
    if (!sequences) {
        return 1; // Exit if dataset loading fails
    }

    // Create and initialize models
    HMM* original_model = HMM_random_create(num_states, num_observations, "Original Model");
    HMM* hypothesis_model = HMM_random_create(num_states, num_observations, "Hypothesis Model");

    /**
    // Run Baum-Welch learning
    printf("Running Baum-Welch...\n");
    HMM* learned_model_bw = BW_learn(hypothesis_model, num_sequences, sequences, sequence_length, epsilon, "logs", "BW_learned_model.txt");
    evaluate_model(learned_model_bw, original_model, sequences, num_sequences, sequence_length, "logs/BW_evaluation.txt");
    HMM_save(learned_model_bw, "BW_learned_model.txt");
    HMM_destroy(learned_model_bw);
    */

    // Run EM-BDD learning
    printf("Running EM-BDD...\n");
    HMM* learned_model_embdd = EMBDD_learn(hypothesis_model, num_sequences, sequences, sequence_length, epsilon, "logs", "EMBDD_learned_model.txt");
    evaluate_model(learned_model_embdd, original_model, sequences, num_sequences, sequence_length, "logs/EMBDD_evaluation.txt");
    HMM_save(learned_model_embdd, "EMBDD_learned_model.txt");
    HMM_destroy(learned_model_embdd);

    // Clean up
    HMM_destroy(original_model);
    HMM_destroy(hypothesis_model);
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);

    printf("Execution completed. Results stored in logs folder.\n");
    return 0;
}
