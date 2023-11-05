#include "HMM.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

/**
 * @brief Create a Hidden Markov Model (HMM) with the given parameters.
 *
 * @param N     The number of states in the HMM.
 * @param M     The number of observations in the HMM.
 * @param name  The name or identifier for the HMM.
 * @return      A pointer to the created HMM.
 */
HMM* HMM_create(int N, int M, const char *name) {
    HMM *hmm = malloc(sizeof(HMM));
    if (hmm == NULL) {
        // Handle memory allocation error
        return NULL;
    }

    hmm->N = N;
    hmm->M = M;
    hmm->name = name;

    // Allocate memory for A, B, and C matrices/vectors
    hmm->A = (double **)malloc(N * sizeof(double *));
    hmm->B = (double **)malloc(N * sizeof(double *));
    hmm->C = (double *)malloc(N * sizeof(double));

    if (hmm->A == NULL || hmm->B == NULL || hmm->C == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(hmm->A);
        free(hmm->B);
        free(hmm->C);
        free(hmm);
        return NULL;
    }

    // ToDo: Remove this
    // Initialize A B and C 
    for (int i = 0; i < N; i++) {
        hmm->A[i] = (double *)malloc(N * sizeof(double));
        hmm->B[i] = (double *)malloc(M * sizeof(double));
        for (int j = 0; j < N; j++) {
            hmm->A[i][j] = 0.0;
        }
        for (int j = 0; j < M; j++) {
            hmm->B[i][j] = 0.0;
        }
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
        fprintf(stderr, "Invalid hmm\n");
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
void validate_hmm(const HMM *hmm) {
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
            assert(hmm->A[i][j] >= 0.0 && hmm->A[i][j] <= 1.0); // Raise an error if not in [0, 1]
        }
        for (int j = 0; j < M; j++) {
            assert(hmm->B[i][j] >= 0.0 && hmm->B[i][j] <= 1.0); // Raise an error if not in [0, 1]
        }
        assert(hmm->C[i] >= 0.0 && hmm->C[i] <= 1.0); // Raise an error if not in [0, 1]
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

double probability_single_sequence(const HMM *hmm, const int *observations, int T) {
    int N = hmm->N;
    int M = hmm->M;

    // Allocate memory for the forward probabilities matrix
    double **alpha = (double **)malloc(T * sizeof(double *));
    for (int t = 0; t < T; t++) {
        alpha[t] = (double *)malloc(N * sizeof(double));
    }

    // Initialize the first row of alpha with initial probabilities and observations
    for (int i = 0; i < N; i++) {
        alpha[0][i] = hmm->C[i] * hmm->B[i][observations[0]];
    }

    // Calculate the forward probabilities for the rest of the sequence
    for (int t = 1; t < T; t++) {
        for (int j = 0; j < N; j++) {
            alpha[t][j] = 0.0;
            for (int i = 0; i < N; i++) {
                alpha[t][j] += alpha[t - 1][i] * hmm->A[i][j];
            }
            alpha[t][j] *= hmm->B[j][observations[t]];
        }
    }

    // Calculate the total probability of the sequence
    double p = 0.0;
    for (int i = 0; i < N; i++) {
        p += alpha[T - 1][i];
    }

    // Free allocated memory
    for (int t = 0; t < T; t++) {
        free(alpha[t]);
    }
    free(alpha);

    return p;
}


double log_likelihood_forward(const HMM *hmm, const int *observations, int T) {
    int N = hmm->N;

    // Allocate memory for the forward probabilities matrix and scaling factors
    double **alpha = (double **)malloc(T * sizeof(double *));
    double *scale_factors = (double *)malloc(T * sizeof(double));
    for (int t = 0; t < T; t++) {
        alpha[t] = (double *)malloc(N * sizeof(double));
    }

    // Initialize the first row of alpha with initial probabilities and observations
    double scale = 0.0;
    for (int i = 0; i < N; i++) {
        alpha[0][i] = hmm->C[i] * hmm->B[i][observations[0]];
        scale += alpha[0][i];
    }
    scale_factors[0] = scale;

    // Scale the alpha values at time t=0
    for (int i = 0; i < N; i++) {
        alpha[0][i] /= scale;
    }

    // Calculate the scaled forward probabilities for the rest of the sequence
    for (int t = 1; t < T; t++) {
        scale = 0.0;
        for (int j = 0; j < N; j++) {
            alpha[t][j] = 0.0;
            for (int i = 0; i < N; i++) {
                alpha[t][j] += alpha[t - 1][i] * hmm->A[i][j];
            }
            alpha[t][j] *= hmm->B[j][observations[t]];
            scale += alpha[t][j];
        }
        scale_factors[t] = scale;

        // Scale the alpha values at time t
        for (int j = 0; j < N; j++) {
            alpha[t][j] /= scale;
        }
    }

    // Calculate the log-likelihood
    double log_likelihood = 0.0;
    for (int t = 0; t < T; t++) {
        log_likelihood += log(scale_factors[t]);
    }

    // Free allocated memory
    for (int t = 0; t < T; t++) {
        free(alpha[t]);
    }
    free(alpha);
    free(scale_factors);

    return log_likelihood;
}


// Utility function to choose a state or observation based on a probability distribution
static int choose_from_distribution(double *probabilities, int size) {
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
int* HMM_generate_sequence(const HMM *hmm, int T) {
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