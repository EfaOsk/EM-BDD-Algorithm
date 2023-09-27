#include "HMM.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>


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
void HMM_destroy(HMM *hmm) {
    if (hmm != NULL) {
        for (int i = 0; i < hmm->N; i++) {
            free(hmm->A[i]);
            free(hmm->B[i]);
        }
        free(hmm->A);
        free(hmm->B);
        free(hmm->C);
        free(hmm);
    }
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
