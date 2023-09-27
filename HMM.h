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
void HMM_destroy(HMM *hmm);
void HMM_print(const HMM *hmm);
void validate_hmm(const HMM *hmm);
void draw_hmm(const HMM *hmm, const char *dotFileName);

#endif // HMM_H
