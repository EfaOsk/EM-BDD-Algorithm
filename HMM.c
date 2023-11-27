#include "HMM.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "helpers.h"

/**
 * @brief Create a Hidden Markov Model (HMM) with the given parameters.
 *
 * @param N     The number of states in the HMM.
 * @param M     The number of observations in the HMM.
 * @param name  The name or identifier for the HMM.
 * @return      A pointer to the created HMM.
 */
HMM* HMM_create(int N, int M, const char *name) 
{
    HMM *hmm = (HMM *)malloc(sizeof(HMM));
    if (hmm == NULL) {
        perror("Failed to allocate memory for HMM structure.\n");
        return NULL;
    }

    hmm->N = N;
    hmm->M = M;
    hmm->name = strdup(name);  // Allocate and copy name

    hmm->A = allocate_matrix(N, N, 0.0);
    hmm->B = allocate_matrix(N, M, 0.0);
    hmm->C = (double *)malloc(N * sizeof(double));
    if (hmm->A == NULL || hmm->B == NULL || hmm->C == NULL) {
        perror("Failed to allocate memory for HMM matrices/vectors.\n");
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        hmm->C[i] = 0.0;
    }

    return hmm;
}


/**
 * @brief Destroy a Hidden Markov Model (HMM) and free its memory.
 * 
 * @param hmm A pointer to the HMM to be destroyed.
 */
void HMM_destroy(const HMM *hmm) {
    HMM *modifiable_hmm = (HMM *)hmm;

    if (modifiable_hmm != NULL) {
        for (int i = 0; i < modifiable_hmm->N; i++) {
            free(modifiable_hmm->A[i]);
            free(modifiable_hmm->B[i]);
        }
        free(modifiable_hmm->A);
        free(modifiable_hmm->B);
        free(modifiable_hmm->C);
        free(modifiable_hmm);
    }

}

void HMM_copy(HMM* dest, const HMM* src) {
        
    // Copy basic fields
    dest->N = src->N;
    dest->M = src->M;
    
    // Copy name
    dest->name = strdup(src->name);  // Remember to free the old name in dest if it was dynamically allocated

    // Copy A matrix
    for (int i = 0; i < src->N; ++i) {
        memcpy(dest->A[i], src->A[i], src->N * sizeof(double));
    }
    
    // Copy B matrix
    for (int j = 0; j < src->N; ++j) {
        memcpy(dest->B[j], src->B[j], src->M * sizeof(double));
    }
    
    // Copy C array
    memcpy(dest->C, src->C, src->N * sizeof(double));
}


/**
 * @brief Print the details of a Hidden Markov Model (HMM) to the console.
 * 
 * @param hmm A pointer to the HMM
 */
void HMM_print(const HMM *hmm) {
    if (hmm == NULL) {
        perror("Invalid hmm\n");
        return;
    }

    printf("HMM Name: %s\n", hmm->name);
    printf("Number of States (N): %d\n", hmm->N);
    printf("Number of Observations (M): %d\n", hmm->M);

    // Print A, B, and C 
    printf("Transition Probability Matrix (A):\n");
    for (int i = 0; i < hmm->N; i++) {
        for (int j = 0; j < hmm->N; j++) {
            printf("%lf ", hmm->A[i][j]);
        }
        printf("\n");
    }

    printf("Observation Probability Matrix (B):\n");
    for (int i = 0; i < hmm->N; i++) {
        for (int j = 0; j < hmm->M; j++) {
            printf("%lf ", hmm->B[i][j]);
        }
        printf("\n");
    }

    printf("Initial State Probability (C):\n");
    for (int i = 0; i < hmm->N; i++) {
        printf("%lf ", hmm->C[i]);
    }
    printf("\n\n");
}


/**
 * Saves an HMM model to a file.
 * 
 * This function writes the HMM's parameters (number of states, number of observations,
 * transition probabilities, observation probabilities, initial state probabilities,
 * and name) to a file in a plain text format.
 * 
 * @param hmm A pointer to the HMM structure to be saved.
 * @param filename The name of the file to which the HMM data will be written.
 *                 If the file already exists, it will be overwritten.
 */
void HMM_save(const HMM *hmm, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    // Write the number of states and number of observations
    fprintf(file, "%d %d\n", hmm->N, hmm->M);

    // Write the name of the HMM
    fprintf(file, "%s\n", hmm->name);

    // Write the initial state probability vector
    for (int i = 0; i < hmm->N; i++) {
        fprintf(file, "%.6f ", hmm->C[i]);
    }
    fprintf(file, "\n");

    // Write the state transition probability matrix
    for (int i = 0; i < hmm->N; i++) {
        for (int j = 0; j < hmm->N; j++) {
            fprintf(file, "%.6f ", hmm->A[i][j]);
        }
        fprintf(file, "\n");
    }

    // Write the observation probability matrix
    for (int i = 0; i < hmm->N; i++) {
        for (int j = 0; j < hmm->M; j++) {
            fprintf(file, "%.6f ", hmm->B[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}


/**
 * Loads an HMM model from a file.
 * 
 * @param filename The name of the file from which the HMM data will be read.
 * @return A pointer to the newly created HMM structure, or NULL if the file
 *         could not be opened, the data could not be read, or memory allocation
 *         failed.
 */
HMM* HMM_load(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    int N, M;
    // Read the number of states and number of observations
    if (fscanf(file, "%d %d\n", &N, &M) != 2) {
        fclose(file);
        return NULL;
    }

    // Read the name of the HMM
    char nameBuffer[256]; // Assuming the name will not exceed 255 characters
    if (fgets(nameBuffer, sizeof(nameBuffer), file) == NULL) {
        fclose(file);
        return NULL;
    }
    // Remove possible newline character read by fgets
    nameBuffer[strcspn(nameBuffer, "\r\n")] = 0;

    // Allocate memory for HMM
    HMM* hmm = HMM_create(N, M, nameBuffer);
    if (hmm == NULL) {
        fclose(file);
        return NULL;
    }

    // Read the initial state probability vector
    for (int i = 0; i < N; i++) {
        if (fscanf(file, "%lf", &(hmm->C[i])) != 1) {
            HMM_destroy(hmm);
            fclose(file);
            return NULL;
        }
    }

    // Read the state transition probability matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fscanf(file, "%lf", &(hmm->A[i][j])) != 1) {
                HMM_destroy(hmm);
                fclose(file);
                return NULL;
            }
        }
    }

    // Read the observation probability matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (fscanf(file, "%lf", &(hmm->B[i][j])) != 1) {
                HMM_destroy(hmm);
                fclose(file);
                return NULL;
            }
        }
    }

    fclose(file);
    return hmm;
}

/**
 * @brief Validate a Hidden Markov Model (HMM) to ensure it meets certain criteria.
 * This function raises an error using the assert macro if the HMM is not valid.
 *
 * @param hmm   A pointer to the HMM.
 */
void HMM_validate(const HMM *hmm) {
    int N = hmm->N;
    int M = hmm->M;

    // Check the sum of each row in matrix A
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += hmm->A[i][j];
        }
        assert(fabs(sum - 1.0) < 1e-6); // Raise an error if the sum is not close to 1.0
    }

    // Check the sum of each row in matrix B
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < M; j++) {
            sum += hmm->B[i][j];
        }
        assert(fabs(sum - 1.0) < 1e-6); // Raise an error if the sum is not close to 1.0
    }

    // Check the sum of elements in vector C
    double c_sum = 0.0;
    for (int i = 0; i < N; i++) {
        c_sum += hmm->C[i];
    }
    assert(fabs(c_sum - 1.0) < 1e-6); // Raise an error if the sum is not close to 1.0

    // Check that all probabilities in matrices A and B are in [0, 1]
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            assert(hmm->A[i][j] >= 0.0-1e-6 && hmm->A[i][j] <= 1.0+1e-6); // Raise an error if not in [0, 1]
        }
        for (int j = 0; j < M; j++) {
            assert(hmm->B[i][j] >= 0.0-1e-6 && hmm->B[i][j] <= 1.0+1e-6); // Raise an error if not in [0, 1]
        }
        assert(hmm->C[i] >= 0.0-1e-6 && hmm->C[i] <= 1.0+1e-6); // Raise an error if not in [0, 1]
    }
}


/**
 * Generate a Graphviz DOT file representation of a Hidden Markov Model (HMM).
 *
 * @param hmm         A pointer to the HMM structure to be represented.
 * @param dotFileName The name of the output DOT file.
 */
void draw_hmm(const HMM *hmm, const char *dotFileName) {
    FILE *dotFile = fopen(dotFileName, "w");

    if (dotFile == NULL) {
        perror("Error opening DOT file");
        return;
    }

    // Write DOT file header
    fprintf(dotFile, "digraph HMM {\n");

    // Write nodes for states
    fprintf(dotFile, "  node [shape=circle];\n");
    for (int i = 0; i < hmm->N; i++) {
        fprintf(dotFile, "  S%d [label=\"S%d\"];\n", i, i);
    }

    // Write initial state probabilities
    fprintf(dotFile, "  node [shape=box, style=filled, color=lightblue];\n");
    fprintf(dotFile, "  I [label=\"Initial Probabilities\"];\n");
    for (int i = 0; i < hmm->N; i++) {
        fprintf(dotFile, "  I -> S%d [label=\"%.2f\"];\n", i, hmm->C[i]);
    }

    // Write transition probabilities
    fprintf(dotFile, "  node [shape=box, style=filled, color=lightgreen];\n");
    fprintf(dotFile, "  T [label=\"Transition Probabilities\"];\n");
    for (int i = 0; i < hmm->N; i++) {
        for (int j = 0; j < hmm->N; j++) {
            fprintf(dotFile, "  S%d -> S%d [label=\"%.2f\"];\n", i, j, hmm->A[i][j]);
        }
    }

    // Write observation probabilities
    fprintf(dotFile, "  node [shape=box, style=filled, color=lightcoral];\n");
    fprintf(dotFile, "  O [label=\"Observation Probabilities\"];\n");
    for (int i = 0; i < hmm->N; i++) {
        for (int j = 0; j < hmm->M; j++) {
            fprintf(dotFile, "  S%d -> O%d [style=\"dashed\", label=\"%.2f\", dir=\"none\", color=\"gray\"];\n", i, j, hmm->B[i][j]);
        }
    }

    // Write DOT file footer
    fprintf(dotFile, "}\n");

    fclose(dotFile);
}


/**
 * Generate a Random Hidden Markov Model (HMM) for a given number of states (N).
 *
 * @param N     The number of states in the HMM.
 * @param M     The number of observations in the HMM.
 * @param name  The name or identifier for the HMM.
 * @return      A pointer to the randomly created HMM.
 */
HMM* HMM_random_create(int N, int M, const char *name) {

    HMM* hmm = HMM_create(N, M, name);
    if (hmm == NULL) {
        return NULL;
    }

    // Initialize the transition probability matrix A
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            hmm->A[i][j] = (double)rand() / RAND_MAX;
            sum += hmm->A[i][j];
        }
        // Normalize to make the sum of probabilities equal to 1
        for (int j = 0; j < N; ++j) {
            hmm->A[i][j] /= sum;
        }
    }

    // Initialize the observation probability matrix B
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < M; ++j) {
            hmm->B[i][j] = (double)rand() / RAND_MAX;
            sum += hmm->B[i][j];
        }
        // Normalize
        for (int j = 0; j < M; ++j) {
            hmm->B[i][j] /= sum;
        }
    }

    // Initialize the initial state probability vector C
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        hmm->C[i] = (double)rand() / RAND_MAX;
        sum += hmm->C[i];
    }
    // Normalize
    for (int i = 0; i < N; ++i) {
        hmm->C[i] /= sum;
    }

    return hmm;
}


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


double log_likelihood(const HMM *hmm, const int *observations, int T)
{
    // todo
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


// Utility function to choose a state or observation based on a probability distribution
static int choose_from_distribution(double *probabilities, int size)
{
    double r = (double)rand() / (double)RAND_MAX;
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += probabilities[i];
        if (r <= sum) {
            return i;
        }
    }
    // In case of floating point arithmetic issues
    return size - 1; 
}


// Function to generate a sequence of observations from a HMM model
int* HMM_generate_sequence(const HMM *hmm, int T) 
{
    if (hmm == NULL || T <= 0) {
        return NULL;
    }

    // Allocate memory for the sequence
    int *sequence = (int*)malloc(T * sizeof(int));
    if (sequence == NULL) {
        fprintf(stderr, "Memory allocation failed for sequence.\n");
        return NULL;
    }

    // Start from the initial state based on initial state probabilities C
    int state = choose_from_distribution(hmm->C, hmm->N);

    for (int t = 0; t < T; t++) {
        // Generate the observation based on the emission probabilities B
        sequence[t] = choose_from_distribution(hmm->B[state], hmm->M);
        
        // Transition to the next state based on the transition probabilities A
        state = choose_from_distribution(hmm->A[state], hmm->N);
    }

    return sequence;
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

        printf("\timprovement: %f\n", prob_new-prob_priv);
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
    

    // Open the result file in append mode
    FILE *result_fp = fopen(result_file, "a");
    if (result_fp == NULL) {
        perror("Error opening result file");
        // Handle the error, possibly by cleaning up and returning
        HMM_destroy(model);
        fclose(log_file);
        return NULL;
    }

    // Append the results to the result file
    fprintf(result_fp, "%d, %f, %f", iteration, prob_new, prob_new - prob_original);


    // Close the result file
    fclose(result_fp);


    // Step 6: Return the learned model
    
    return model;
}