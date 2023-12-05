#include "HMM.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "helpers.h"
#include <sys/stat.h>
#include <sys/types.h>


double **forward(const HMM *hmm, const int *observations, int T)
{
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


double **forward_log(const HMM *hmm, const int *observations, int T)
{
    int N = hmm->N;

    double **alpha = allocate_matrix(N, T, -1);

    // Initialize the first row of alpha with initial probabilities and observations
    for (int i = 0; i < N; i++) {
        alpha[i][0] = log(hmm->C[i]) + log(hmm->B[i][observations[0]]);
    }

    // Calculate the forward probabilities for the rest of the sequence
    for (int t = 1; t < T; t++) {
        for (int j = 0; j < N; j++) {
            alpha[j][t] = -INFINITY;
            for (int i = 0; i < N; i++) {
                alpha[j][t] = log_sum_exp(alpha[j][t], alpha[i][t-1]+log(hmm->A[i][j]));
            }
            alpha[j][t] += log(hmm->B[j][observations[t]]);
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


double **backward_log(const HMM *hmm, const int *observations, int T)
{
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
                beta[s][t] = log_sum_exp(beta[s][t], log(hmm->A[s][s0]) + log(hmm->B[s0][observations[t+1]]) + beta[s0][t+1]);
            }
        }
    }

    return beta;
}


double **backward(const HMM *hmm, const int *observations, int T)
{
    int N = hmm->N;

    double **beta = allocate_matrix(N, T, -1);

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

    return beta;
}


void calculate_Xi(HMM *hmm, double ***Xi, double **alpha, double **beta, int *observations, int T)
{
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


void calculate_gamma(HMM *hmm, double **gamma, double **alpha, double **beta, int *observations, int T)
{   
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


HMM* HMM_update(HMM *hmm, double **alpha, double **beta, int *observations, int T) 
{
    int N = hmm->N;
    int M = hmm->M;

    double min_p_f = 0.00001;

    HMM *new_hmm = HMM_create(N, M, "Updated model");

    double ***Xi = allocate_3D_matrix(N, N, T-1, -1);

    double **gamma = allocate_matrix(N, T, -1);


    for (int i = 0; i < N; i++) {
        gamma[i] = malloc(T * sizeof(double));
        if (gamma[i] == NULL) {
            perror("Failed to allocate memory for gamma");
        }
    }

    calculate_Xi(hmm, Xi, alpha, beta, observations, T);
    calculate_gamma(hmm, gamma, alpha, beta, observations, T);

    for (int s = 0; s < N; s++) {
        double gamma_sum = 0.0; // Pr(observations | hmm)
        for (int t = 0; t < T-1; t++) {
            gamma_sum += gamma[s][t];
        }

        // update Transition probability
        for (int s0 = 0; s0 < N; s0++) {
            double xi_sum  = 0.0;
            for (int t = 0; t < T-1; t++) {
                xi_sum  += Xi[s][s0][t];
            }
            new_hmm->A[s][s0] = (xi_sum + min_p_f)/ ( gamma_sum + min_p_f*N);
        }

        gamma_sum += gamma[s][T-1] + min_p_f*M;
        // update Observation probability
        for (int o = 0; o < M; o++) {
            double gamma_obs_sum  = 0.0;
            for (int t = 0; t < T; t++) {
                if (o == observations[t]) {
                    gamma_obs_sum  += gamma[s][t];
                }
            }
            new_hmm->B[s][o] = (gamma_obs_sum + min_p_f)  / gamma_sum ;
        }
    }

    // update Initial state probability
    for (int s = 0; s < N; s++) {
        new_hmm->C[s] = (gamma[s][0] + min_p_f) / (1 + N*min_p_f);
    }

    // Free Xi
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            free(Xi[i][j]);
        }
        free(Xi[i]);
    }
    free(Xi);

    // Free gamma
    for (int i = 0; i < N; i++) {
        free(gamma[i]);
    }
    free(gamma);
    return new_hmm;
}



HMM* HMM_learn(HMM *hypothesis_hmm, int T, int O[T], double epsilon, const char *logs_folder, const char *result_file)
{
    /*

    EM on BDD

        (repeat steps 3-5 until converged)

        (3) E-step
            (a) Backward
            (b) Forward
        
        (4) M-step
            (a) update M
        
        (5) Calculate the log-likelyhood of M

        (6) retrun M

    */
    int N = hypothesis_hmm->N;
    int M = hypothesis_hmm->M;
    HMM *model = HMM_create(N, M, "model");
    HMM_copy(model, hypothesis_hmm);
    HMM_validate(model);

    // for loging
    int iteration = 0;
    char log_filename[256];
    char model_filename[256];

    // Construct the log filename
    sprintf(log_filename, "%s/log.txt", logs_folder);

    // Open the log file
    FILE *log_file = fopen(log_filename, "w");
    if (log_file == NULL) {
        perror("Error opening log file");
        return NULL;
    }

    sprintf(model_filename, "%s/models", logs_folder);
    mkdir(model_filename, 0777);
    double prob_priv, prob_original, prob_new;
    prob_original = log_likelihood_forward(model, O, T);
    prob_priv = prob_original;
    int converged = 0;


    double **beta = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        beta[i] = (double *)malloc(T * sizeof(double));
    }

    while (!converged)
    {   
        clock_t start_time = clock();

        // Step 3: E-step
    
        // Step 3 (a) : Backward

        // Step 3 (b) : Forward

        double **alpha = forward(model, O, T);

        double **beta = backward(model, O, T);


        // Step 3 (c) : Conditional Expectations

        // Step 4: M-step
        // Step 4 (a) : update M

        const HMM *new_hmm = HMM_update(model, alpha, beta, O, T);

        prob_new = log_likelihood_forward(new_hmm, O, T);

        // Free allocated memory
        for (int i = 0; i < N; i++) {
            free(alpha[i]);
            free(beta[i]);
        }
        free(alpha);
        free(beta);
        
        HMM_copy(model, new_hmm); 
        HMM_destroy(new_hmm);

        // HMM_print(model);
        HMM_validate(model);
        clock_t end_time = clock();
        double iteration_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        // printf("\timprovement: %f\n", prob_new-prob_priv);
        fprintf(log_file, "Iteration: %d, Log Likelihood: %f, Improvement: %f, Time: %f\n",
                iteration, prob_new, prob_new-prob_priv, iteration_time);
        fflush(log_file); 

        sprintf(model_filename, "%s/models/model_%d", logs_folder, iteration);
        HMM_save(model, model_filename); 
        if (prob_new <= prob_priv+epsilon) {
            converged = 1;
        }

        prob_priv = prob_new;
        iteration++;


    }

    fclose(log_file);

    // Open the result file in append mode
    FILE *result_fp = fopen(result_file, "a");
    if (result_fp == NULL) {
        perror("Error opening result file");
        return NULL;
    }

    // Append the results to the result file
    fprintf(result_fp, "%d, %f, %f\n", iteration, prob_new, prob_new - prob_original);


    // Close the result file
    fclose(result_fp);


    // Step 6: Return the learned model
    
    return model;
}