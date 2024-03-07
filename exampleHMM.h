#ifndef example_HMM
#define example_HMM


#include <stdio.h>
#include "HMM.h"

#define NUM_MODELS 10 // Number of HMM models to create
#define NUM_LARGE_MODELS 10

extern HMM* example_models[NUM_MODELS];

HMM** initialize_example_models_large();
HMM** initialize_example_models_small();
HMM** initialize_example_models_old();
void free_example_models(HMM **models);

#endif