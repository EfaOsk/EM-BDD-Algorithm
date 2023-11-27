#include "HMM.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>


/**
 * Generate a Graphviz DOT file representation of a Hidden Markov Model (HMM).
 *
 * @param hmm         A pointer to the HMM structure to be represented.
 * @param dotFileName The name of the output DOT file.
 */
void HMM_draw(const HMM *hmm, const char *dotFileName)
{
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
HMM* HMM_random_create(int N, int M, const char *name)
{

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
