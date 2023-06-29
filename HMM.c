#include "HMM.h"
#include<string.h>


static struct HMM init(const int n_s, const int n_o, struct DoubleMatrix A, struct DoubleMatrix B, struct DoubleArray C, const char *name)
{

    struct HMM ret=(struct HMM){
        .n_s=n_s,
        .n_o=n_o,
        .A=A,
        .B=B,
        .C=C,
		.name=strdup(name)
        // .print=&print
	};
    return ret;
}
const struct HMMClass HMM={.init=&init};