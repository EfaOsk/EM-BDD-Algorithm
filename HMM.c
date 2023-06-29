#include "HMM.h"
#include<string.h>

/**
 * Initializes a Hidden Markov Model (HMM) structure.
 *
 * @param N     The number of states.
 * @param M     The number of observations.
 * @param A     The state transition matrix.
 * @param B     The observation emission matrix.
 * @param C     The initial state distribution.
 * @param name  The name of the HMM.
 *
 * @return The initialized HMM structure.
 */
static struct HMM init(const int N, const int M, struct DoubleMatrix A, struct DoubleMatrix B, struct DoubleArray C, const char *name)
{

    struct HMM ret=(struct HMM){
        .N=N,
        .M=M,
        .A=A,
        .B=B,
        .C=C,
		.name=strdup(name)
        // .print=&print
	};
    return ret;
}
const struct HMMClass HMM={.init=&init};