#include "HMM.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>


void forward(const HMM *hmm, const int *observations, int T, double **alpha)
{
    int N = hmm->N;

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
}


double **forward_log(const HMM *hmm, const int *observations, int T)
{
    int N = hmm->N;

    double **alpha = allocate_matrix(N, T, -1);

    // Initialize the first row of alpha with initial probabilities and observations
    for (int i = 0; i < N; i++) {
        if ( (hmm->C[i] > 0.0) && (hmm->B[i][observations[0]] > 0.0) ) {
            alpha[i][0] = log(hmm->C[i]) + log(hmm->B[i][observations[0]]);
        } else {
            alpha[i][0] = -50;
        }
    }

    // Calculate the forward probabilities for the rest of the sequence
    for (int t = 1; t < T; t++) {
        for (int j = 0; j < N; j++) {
            alpha[j][t] = -INFINITY;
            for (int i = 0; i < N; i++) {
                if (hmm->A[i][j] > 0.0) {
                    alpha[j][t] = log_sum_exp(alpha[j][t], alpha[i][t-1]+log(hmm->A[i][j]));
                } else {
                    alpha[j][t] = log_sum_exp(alpha[j][t], alpha[i][t-1]-50);
                }
            }
            if (hmm->B[j][observations[t]] > 0.0) {
                alpha[j][t] += log(hmm->B[j][observations[t]]);
            } else {
                alpha[j][t] += -50;
            }
        }

    }

    return alpha;
}


double log_likelihood_forward(const HMM *hmm, const int *observations, int T)
{
    double **alpha_log = forward_log(hmm, observations, T);
    // Calculate the log-likelihood
    double log_likelihood = -INFINITY;
    for (int j = 0; j < hmm->N; j++) {
        log_likelihood = log_sum_exp(log_likelihood, alpha_log[j][T-1]);
    }

    free_matrix(alpha_log, hmm->N);
    return log_likelihood;
}


double log_likelihood_forward_multiple(const HMM *hmm, int **observations, int num_sequences, int T) {
    double total_log_likelihood = 0;

    for (int seq = 0; seq < num_sequences; seq++) {
        double seq_log_likelihood = log_likelihood_forward(hmm, observations[seq], T);
        total_log_likelihood += seq_log_likelihood;
    }

    return total_log_likelihood;
}


double **backward_log(const HMM *hmm, const int *observations, int T) {
    int N = hmm->N;

    double **beta = allocate_matrix(N, T, -1);

    // Initialize the first row of beta with initial probabilities and observations
    for (int i = 0; i < N; i++) {
        beta[i][T-1] = 0.0;
    }

    // Calculate the backward probabilities for the rest of the sequence
    for (int t = T-2; t >= 0; t--) {
        for (int s = 0; s < N; s++) {
            beta[s][t] = -INFINITY;
            for (int s0 = 0; s0 < N; s0++) {
                if (hmm->A[s][s0] > 0.0 && hmm->B[s0][observations[t+1]] > 0.0) {
                    beta[s][t] = log_sum_exp(beta[s][t], log(hmm->A[s][s0]) + log(hmm->B[s0][observations[t+1]]) + beta[s0][t+1]);
                } else {
                    beta[s][t] = -50;
                }
                
            }
        }
    }

    return beta;
}


void backward(const HMM *hmm, const int *observations, int T, double **beta) {
    int N = hmm->N;

    // Initialize the first row of beta with initial probabilities and observations
    for (int s = 0; s < N; s++) {
        beta[s][T-1] = 1.0;
    }

    // Calculate the scaled forward probabilities for the rest of the sequence
    for (int t = T-2; t >= 0; t--) {
        for (int s = 0; s < N; s++) {
            beta[s][t] = 0.0;
            for (int s0 = 0; s0 < N; s0++) {
                beta[s][t] += hmm->A[s][s0] * hmm->B[s0][observations[t+1]] * beta[s0][t+1];
            }
        }

    }

}


void calculate_Xi(HMM *hmm, double ***Xi, double **alpha, double **beta, int *observations, int T) {
    int N = hmm->N;
    
    for (int t = 0; t < T-1; t++) {
        double PrO = 0.0; // Pr(observations | hmm)
        for (int u = 0; u < N; u++) {
            PrO += alpha[u][t] * beta[u][t];
        }

        for (int u = 0; u < N; u++) {
            for (int u0 = 0; u0 < N; u0++) {
                Xi[u][u0][t] = alpha[u][t] * hmm->A[u][u0] * hmm->B[u0][observations[t+1]] * beta[u0][t+1];
                Xi[u][u0][t] /= PrO;
            }
        }
    }
}


void calculate_gamma(HMM *hmm, double **gamma, double **alpha, double **beta, int *observations, int T) {   
    int N = hmm->N;
    
    for (int t = 0; t < T; t++) {
        double PrO = 0.0; // Pr(observations | hmm)
        for (int u = 0; u < N; u++) {
            PrO += alpha[u][t]* beta[u][t];
        }

        for (int u = 0; u < N; u++) {
            gamma[u][t] = alpha[u][t] * beta[u][t];
            gamma[u][t] /= PrO;
        }
    }
}



HMM* HMM_update(HMM *hmm, int **observations, int num_sequences, int T, double **alpha, double **beta, double ****Xi, double ***gamma) {
    int N = hmm->N;
    int M = hmm->M;
    double min_p_f = 0.00001;

    HMM *new_hmm = HMM_create(N, M, "Updated model");

    // double ****Xi = allocate_4D_matrix(num_sequences, N, N, T-1, 0.0);
    // double ***gamma = allocate_3D_matrix(num_sequences, N, T, 0.0);


    // accumulated Xi and gamma over all sequences
    for (int seq = 0; seq < num_sequences; seq++) {
        forward(hmm, observations[seq], T, alpha);
        backward(hmm, observations[seq], T, beta);
        calculate_Xi(hmm, Xi[seq], alpha, beta, observations[seq], T);
        calculate_gamma(hmm, gamma[seq], alpha, beta, observations[seq], T);
    }


    for (int s = 0; s < N; s++) {
        double gamma_sum = 0.0; 
        for (int k = 0; k < num_sequences; k++) {
            gamma_sum += gamma[k][s][0];
        }
        new_hmm->C[s] = (gamma_sum + min_p_f) / (num_sequences+N*min_p_f);
    }

    for (int s = 0; s < N; s++) {
        double gamma_sum = 0.0; // Pr(observations | hmm)
        for (int k = 0; k < num_sequences; k++) {
            for (int t = 0; t < T-1; t++) {
                gamma_sum += gamma[k][s][t];
            }
        }

        // update Transition probability
        for (int s0 = 0; s0 < N; s0++) {
            double xi_sum  = 0.0;
            for (int k = 0; k < num_sequences; k++) {
                for (int t = 0; t < T-1; t++) {
                    xi_sum  += Xi[k][s][s0][t];
                }
            }
            new_hmm->A[s][s0] = (xi_sum + min_p_f)/ ( gamma_sum + min_p_f*N);
        }
        for (int k = 0; k < num_sequences; k++) {
            gamma_sum += gamma[k][s][T-1];
        }
        // update Observation probability
        for (int o0 = 0; o0 < M; o0++) {
            double xi_sum  = 0.0;
            for (int k = 0; k < num_sequences; k++) {
                for (int t = 0; t < T; t++) {
                    if (o0 == observations[k][t]){
                        xi_sum += gamma[k][s][t];
                    }
                }
            }
            new_hmm->B[s][o0] = (xi_sum + min_p_f)/ ( gamma_sum + min_p_f*M);
        }

    }
    // free_4D_matrix(Xi, num_sequences, N, N);
    // free_3D_matrix(gamma, num_sequences, N);
    return new_hmm;
}


HMM* BW_learn(HMM *hypothesis_hmm, int num_sequences, int **observations, int T, double epsilon, const char *logs_folder, const char *result_file) {
    int N = hypothesis_hmm->N;
    int M = hypothesis_hmm->M;
    HMM *model = HMM_create(N, M, "model");
    HMM_copy(model, hypothesis_hmm);
    HMM_validate(model);

    char log_filename[256];
    char model_filename[256];
    // Construct the log filename
    sprintf(log_filename, "%s/log.txt", logs_folder);

    // Open the log file
    FILE *log_file = fopen(log_filename, "w");
    if (log_file == NULL) {
        HMM_destroy(model);
        perror("Error opening log file");
        return NULL;
    }

    sprintf(model_filename, "%s/models", logs_folder);
    mkdir(model_filename, 0777);


    double prob_priv, prob_original, prob_new;
    prob_original = log_likelihood_forward_multiple(model, observations, num_sequences, T);
    prob_priv = prob_original;
    int converged = 0;
    int iteration = 0;

    double **alpha = allocate_matrix(N, T, -1);
    double **beta = allocate_matrix(N, T, -1);

    double ****Xi = allocate_4D_matrix(num_sequences, N, N, T-1, 0.0);
    double ***gamma = allocate_3D_matrix(num_sequences, N, T, 0.0);

    printf("Starting Baum-Welch learning process...\n");
    while (!converged) {   
        clock_t start_time = clock();
        
        // M-step: Update HMM parameters using forward and backward probabilities
        HMM *new_hmm = HMM_update(model, observations, num_sequences, T, alpha, beta, Xi, gamma);

        prob_new = log_likelihood_forward_multiple(new_hmm, observations, num_sequences, T);

        HMM_copy(model, new_hmm); 
        HMM_destroy(new_hmm);

        // HMM_print(model);
        HMM_validate(model);

        clock_t end_time = clock();
        double iteration_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        fprintf(log_file, "Iteration %d: Log-likelihood = %f, Improvement = %f, Time = %f seconds\n", 
            ++iteration, prob_new, prob_new - prob_priv, iteration_time);
        fflush(log_file); 

    	// ToDo remove:
        printf("Iteration %d: Log-likelihood = %f, Improvement = %f, Time = %f seconds\n", 
            iteration, prob_new, prob_new - prob_priv, iteration_time);
        
        sprintf(model_filename, "%s/models/model_%d", logs_folder, iteration);
        HMM_save(model, model_filename); 

        if (prob_new <= prob_priv+epsilon) {
            converged = 1;
            printf("Convergence achieved after %d iterations.\n", iteration);
        }
        prob_priv = prob_new;

        reset_matrix(alpha, N, T, -1);
        reset_matrix(beta, N, T, -1);

        reset_4D_matrix(Xi, num_sequences, N, N, T-1, 0.0);
        reset_3D_matrix(gamma, num_sequences, N, T, 0.0);
    }

    free_matrix(alpha, N);
    free_matrix(beta, N);

    free_4D_matrix(Xi, num_sequences, N, N);
    free_3D_matrix(gamma, num_sequences, N);

    fclose(log_file);

    // Open the result file in append mode
    FILE *result_fp = fopen(result_file, "a");
    if (result_fp == NULL) {
        perror("Error opening result file");
        HMM_destroy(model);
        return NULL;
    }

    fprintf(result_fp, "%d, %f, %f\n", iteration, prob_new, prob_new - prob_original);
    fclose(result_fp);

    // Return the learned model
    return model;
}
