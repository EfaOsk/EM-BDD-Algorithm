#include "HMM.h"
#include "helpers.h"
#include <assert.h>


struct HMM HMM1()
{
    //  basic info about the HMM
    const int n_s= 2;
    const int n_o= 2;
    
    double a0[]= {(double) 1/3, (double) 2/3};
    double a1[]= {(double) 2/3, (double) 1/3};
    struct DoubleArray a[2]=
        {(struct DoubleArray){
            .n=n_s,
            .data=a0
        },
        (struct DoubleArray){
            .n=n_s,
            .data=a1
        }
    };
    struct DoubleMatrix A=(struct DoubleMatrix){
        .n=n_s,
        .data=a
	};

    double b0[]= {(double) 1/4, (double) 3/4};
    double b1[]= {(double) 3/4, (double) 1/4};
    struct DoubleArray b[2]=
        {(struct DoubleArray){
            .n=n_o,
            .data=b0
        },
        (struct DoubleArray){
            .n=n_o,
            .data=b1
        }
    };
    struct DoubleMatrix B=(struct DoubleMatrix){
        .n=n_s,
        .data=b
	};


    double c[]= {(double) 1/2, (double) 1/2};
    struct DoubleArray C= (struct DoubleArray){
            .n=n_s,
            .data=c
    };
    
    struct HMM M = HMM.init(n_s, n_o, A, B, C, "HMM1");

    return M;
}