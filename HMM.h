#ifndef HMM_H
#define HMM_H

#include <stdio.h>


typedef struct HMM {
    int N;      // Number of states
    int M;      // Number of observations
    double **A; // Transition probability matrix
    double **B; // Observation probability matrix
    double *C;  // Initial state probability vector
    const char *name;
} HMM;


// Prototypes for HMM management functions (from HMM_management.c)

HMM* HMM_create(int N, int M, const char *name);
void HMM_destroy(const HMM *hmm);
void HMM_copy(HMM* dest, const HMM* src);
void HMM_print(const HMM *hmm);
void HMM_save(const HMM *hmm, const char *filename);
HMM* HMM_load(const char *filename);
void HMM_validate(const HMM *hmm);



// Prototypes for algorithm functions (from HMM_algorithms.c)
double **forward(const HMM *hmm, const int *observations, int T);
double **backward(const HMM *hmm, const int *observations, int T);
double **forward_log(const HMM *hmm, const int *observations, int T);
double **backward_log(const HMM *hmm, const int *observations, int T);
double log_likelihood(const HMM *hmm, const int *observations, int T);
double log_likelihood_forward(const HMM *hmm, const int *observations, int T);
void calculate_Xi(HMM *hmm, double ***Xi, double **alpha, double **beta, int *observations, int T);
void calculate_gamma(HMM *hmm, double **gamma, double **alpha, double **beta, int *observations, int T);
HMM* HMM_update(HMM *hmm, double **alpha, double **beta, int *observations, int T);
HMM* HMM_learn(HMM *hypothesis_hmm, int T, int O[T], double epsilon, const char *logs_folder, const char *result_file);


// Prototypes for utility functions (from HMM_utils.c)

int* HMM_generate_sequence(const HMM *hmm, int T);
static int choose_from_distribution(double *probabilities, int size);
HMM* HMM_random_create(int N, int M, const char *name);
void HMM_draw(const HMM *hmm, const char *dotFileName);


#endif // HMM_H