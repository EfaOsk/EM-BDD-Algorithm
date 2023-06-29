#ifndef HMM_h
#define HMM_h


#include <stdio.h>
#include "helpers.h"


struct HMM {
    const int N;      // Number of states         The set of states S= {0, 1 ... n_s-1}
    const int M;      // Number of observations   The set of observations \Sigma = {0, 1 ... n_o-1}
    struct DoubleMatrix A;  // Transition probability   A[i][j] = propabailty of transition from state i to state j
    struct DoubleMatrix B;  // Observation probability  B[i][j] = propabailty of observing observation j in state i
    struct DoubleArray C;   // Initialstate probability C[i] =  probability of starting in state i ( \pi (i) )

	const char *name;   // name of the model
	// const char *(*print)(struct HMM *self, size_t bufsize,char buf[bufsize]);
};


extern const struct HMMClass {
	struct HMM (*init)(const int N, const int M, struct DoubleMatrix A, struct DoubleMatrix B, struct DoubleArray C, const char *name);
} HMM;



#endif