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

HMM* HMM_create(int N, int M, const char *name);
void HMM_destroy(const HMM *hmm);
void HMM_copy(HMM* dest, const HMM* src);
void HMM_print(const HMM *hmm);
void HMM_save(const HMM *hmm, const char *filename);
HMM* HMM_load(const char *filename);
void validate_hmm(const HMM *hmm);
HMM* HMM_random_create(int N, int M, const char *name);
double probability_single_sequence(const HMM *hmm, const int *observations, int T);
double log_likelihood(const HMM *hmm, const int *observations, int T);
HMM* HMM_learn(HMM *hypothesis_hmm, int T, int O[T], double epsilon, const char *logs_folder, const char *result_file);
void draw_hmm(const HMM *hmm, const char *dotFileName);
int* HMM_generate_sequence(const HMM *hmm, int T);


#endif // HMM_H
