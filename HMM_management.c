#include "HMM.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
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
void HMM_destroy(const HMM *hmm) 
{
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


void HMM_copy(HMM* dest, const HMM* src)
{
        
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
void HMM_print(const HMM *hmm) 
{
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
void HMM_save(const HMM *hmm, const char *filename)
{
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
HMM* HMM_load(const char *filename)
{
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
void HMM_validate(const HMM *hmm)
{
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

