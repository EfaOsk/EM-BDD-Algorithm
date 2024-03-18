#include "HMM.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>

/**
 * Perform the forward algorithm to compute the alpha matrix (probabilities of sequences up to time t).
 * The caller is responsible for freeing the allocated matrix by calling free_matrix.
 */
double **forward(const HMM *hmm, const int *observations, int T) {
    int N = hmm->N;

    double **alpha = allocate_matrix(N, T, -1);

    // Initialize the first row of alpha with initial probabilities and observations
    for (int s = 0; s < N; s++) {
        alpha[s][0] = hmm->C[s] * hmm->B[s][observations[0]];
    }

    // Calculate the forward probabilities for the rest of the sequence
    for (int t = 1; t < T; t++) {
        for (int s = 0; s < N; s++) {
            alpha[s][t] = 0.0;
            for (int s0 = 0; s0 < N; s0++) {
                alpha[s][t] += alpha[s0][t-1] * hmm->A[s0][s];
            }
            alpha[s][t] *= hmm->B[s][observations[t]];
        }
    }

    return alpha;
}

// ... [other functions remain unchanged] ...

/**
 * Updates the HMM model based on the forward and backward algorithms across multiple sequences.
 * This function takes ownership of the new_hmm, which must be destroyed by the caller.
 */
HMM* HMM_update_multiple(HMM *hmm, int **observations, int num_sequences, int T) {
    // ... [rest of the HMM_update_multiple function remains unchanged] ...
}

/**
 * The Baum-Welch learning algorithm for a single observation sequence.
 * The function creates and returns a new HMM model, which the caller is responsible for destroying.
 */
HMM* BW_learn(HMM *hypothesis_hmm, int **observations, int num_sequences, int T, double epsilon, const char *logs_folder, const char *result_file) {
    // ... [initialize variables and log file] ...

    double prob_priv = log_likelihood_forward_multiple(model, observations, num_sequences, T);
    double prob_new;
    
    // Learning loop
    while (!has_converged) {
        // ... [perform the E-step and M-step] ...

        // Destroy the old HMM model and replace it with the updated one
        HMM_destroy(new_hmm);

        // ... [check for convergence] ...
    }

    // Clean up and close log file
    fclose(log_file);

    // Write the results to the result file and clean up
    FILE *result_fp = fopen(result_file, "a");
    if (!result_fp) {
        perror("Error opening result file");
        HMM_destroy(model);
        return NULL;
    }

    fprintf(result_fp, "%d, %f, %f\n", iteration, prob_new, prob_new - prob_priv);
    fclose(result_fp);

    // Return the learned HMM model
    return model;
}
