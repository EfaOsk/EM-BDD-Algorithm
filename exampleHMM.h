#ifndef example_HMM
#define example_HMM


#include <stdio.h>
#include "HMM.h"

#define NUM_MODELS 11 // Number of HMM models to create

extern HMM* example_models[NUM_MODELS];

HMM** initialize_example_models();

#endif