#include "HMM.h"
#include "helpers.h"
#include "exampleHMM.h"
#include "stdlib.h"
#include <string.h>

void initialize_large_model_with_varying_structure(HMM *model, int numStates, int numObservations) {
    // Zero-initialize matrices
    for (int i = 0; i < numStates; i++) {
        memset(model->A[i], 0, numStates * sizeof(double));
        memset(model->B[i], 0, numObservations * sizeof(double));
    }
    memset(model->C, 0, numStates * sizeof(double));

    // Transition matrix with varying out-degree
    for (int i = 0; i < numStates; i++) {
        // Determine out-degree for this state (at least 1, up to numStates)
        int outDegree = (rand() % (numStates - 1)) + 1; // Ensures at least 1 out-degree
        double totalProb = 0.0, prob;
        for (int j = 0; j < outDegree; j++) {
            // Randomly pick a state to transition to, ensuring it's not a self-transition
            int transitionState = (i + j + 1) % numStates;
            prob = (1.0 - totalProb) / (outDegree - j); // Distribute remaining probability among remaining transitions
            model->A[i][transitionState] = prob;
            totalProb += prob;
        }
        model->A[i][i] = 1.0 - totalProb; // Remaining probability assigned to self-transition
    }

    // Observation matrix with varying probabilities
    for (int i = 0; i < numStates; i++) {
        double totalProb = 0.0;
        int significantObservation = rand() % numObservations; // Choose one observation to be more likely than others
        for (int j = 0; j < numObservations; j++) {
            if (j == significantObservation) {
                model->B[i][j] = 0.5 + ((double)rand() / RAND_MAX) * 0.5; // Assign a higher probability to one observation
            } else {
                double prob = ((double)rand() / RAND_MAX) * (0.5 / (numObservations - 1)); // Distribute the rest among others
                model->B[i][j] = prob;
            }
            totalProb += model->B[i][j];
        }
        // Normalize the probabilities to ensure they sum up to 1
        for (int j = 0; j < numObservations; j++) {
            model->B[i][j] /= totalProb;
        }
    }

    // Initial state distribution with varying probabilities
    double totalCProb = 0.0;
    for (int i = 0; i < numStates; i++) {
        double prob;
        if (i < numStates / 2) {
            // Assign higher probabilities to the first half of states
            prob = ((double)rand() / RAND_MAX) * (2.0 / numStates);
        } else {
            // Lower probabilities for the second half
            prob = ((double)rand() / RAND_MAX) * (1.0 / numStates);
        }
        model->C[i] = prob;
        totalCProb += prob;
    }
    // Normalize the initial state probabilities
    for (int i = 0; i < numStates; i++) {
        model->C[i] /= totalCProb;
    }
}


HMM** initialize_example_models_large() {
    HMM **models = (HMM **)malloc(NUM_LARGE_MODELS * sizeof(HMM *));
    // Assuming larger models have more states and observations
    int numStatesArr[NUM_LARGE_MODELS+1] = {10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65};
    int numObservationsArr[NUM_LARGE_MODELS+1] = {5, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

    for (int i = 0; i < NUM_LARGE_MODELS; i++) {
        char modelName[50];
        sprintf(modelName, "Example Large Model %d", i);
        models[i] = HMM_create(numStatesArr[i], numObservationsArr[i], modelName);
        initialize_large_model_with_varying_structure(models[i], numStatesArr[i], numObservationsArr[i]);
    }

    return models;
}

/**
 * @brief Initilizes ten example HMM.
 * 
 * @return HMM** 
 */
HMM** initialize_example_models_small() {
    // Initilizes ten example HMM. 
    
    HMM **models = (HMM **)malloc(NUM_MODELS * sizeof(HMM *)); // Allocate memory for 10 HMM models

    double A0[6][6] = { {0.0, 0.5, 0.5, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.2, 0.8, 0.0},
                        {0.0, 0.0, 0.0, 0.3, 0.7, 0.0},
                        {0.0, 0.0, 0.0, 0.8, 0.2, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
    double B0[6][4] = { {0.2, 0.4, 0.2, 0.2},
                        {0.4, 0.1, 0.3, 0.2},
                        {0.4, 0.1, 0.4, 0.1},
                        {0.1, 0.1, 0.6, 0.2},
                        {0.3, 0.3, 0.2, 0.2},
                        {0.0, 0.3, 0.2, 0.5}};
    double C0[6] = {0.5, 0.1, 0.1, 0.1, 0.1, 0.1};


    models[0] = HMM_create(6, 4, "Example Model 0");

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            models[0]->A[i][j] = A0[i][j];
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            models[0]->B[i][j] = B0[i][j];
        }
    }

    for (int i = 0; i < 6; i++) {
        models[0]->C[i] = C0[i];
    }




    double A1[6][6] = { {0.0, 0.5, 0.5, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.6, 0.4, 0.0},
                        {0.0, 0.2, 0.0, 0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.3, 0.7},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
    double B1[6][5] = { {1.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.8, 0.2, 0.0, 0.0},
                        {0.0, 0.0, 0.2, 0.8, 0.0},
                        {0.1, 0.9, 0.0, 0.0, 0.0},
                        {0.1, 0.0, 0.0, 0.9, 0.0},
                        {0.0, 0.0, 0.0, 0.1, 0.9}};
    double C1[6] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    models[1] = HMM_create(6, 5, "Example Model 1");
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            models[1]->A[i][j] = A1[i][j];
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 5; j++) {
            models[1]->B[i][j] = B1[i][j];
        }
    }

    for (int i = 0; i < 6; i++) {
        models[1]->C[i] = C1[i];
    }



    double A2[4][4] = { {0.0, 0.8, 0.2, 0.0},
                        {0.0, 0.0, 0.0, 1.0},
                        {0.0, 0.8, 0.0, 0.2},
                        {0.0, 0.0, 0.1, 0.9}};
    double B2[4][5] = { {0.5, 0.5, 0.0, 0.0, 0.0},
                        {0.0, 0.1, 0.8, 0.1, 0.0},
                        {0.0, 0.8, 0.1, 0.1, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 1.0}};
    double C2[4] = {1.0, 0.0, 0.0, 0.0};

    models[2] = HMM_create(4, 5, "Example Model 2");

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            models[2]->A[i][j] = A2[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            models[2]->B[i][j] = B2[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        models[2]->C[i] = C2[i];
    }



    double A3[4][4] = {{0.5, 0.1, 0.2, 0.2}, {0.2, 0.3, 0.2, 0.3}, {0.1, 0.2, 0.4, 0.3}, {0.3, 0.2, 0.1, 0.4}};
    double B3[4][5] = {{0.1, 0.2, 0.3, 0.2, 0.2}, {0.3, 0.1, 0.2, 0.2, 0.2}, {0.2, 0.2, 0.1, 0.3, 0.2}, {0.2, 0.1, 0.2, 0.2, 0.3}};
    double C3[4] = {0.3, 0.2, 0.2, 0.3};
    models[3] = HMM_create(4, 5, "Example Model 3");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            models[3]->A[i][j] = A3[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            models[3]->B[i][j] = B3[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        models[3]->C[i] = C3[i];
    }

    
    double A4[4][4] = {{0.4, 0.3, 0.2, 0.1}, {0.2, 0.4, 0.2, 0.2}, {0.1, 0.2, 0.5, 0.2}, {0.3, 0.1, 0.3, 0.3}};
    double B4[4][5] = {{0.2, 0.3, 0.2, 0.2, 0.1}, {0.3, 0.1, 0.3, 0.2, 0.1}, {0.1, 0.3, 0.1, 0.2, 0.3}, {0.2, 0.2, 0.2, 0.3, 0.1}};
    double C4[4] = {0.3, 0.2, 0.1, 0.4};
    models[4] = HMM_create(4, 5, "Example Model 4");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            models[4]->A[i][j] = A4[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            models[4]->B[i][j] = B4[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        models[4]->C[i] = C4[i];
    }


    double A5[5][5] = { {0.2, 0.0, 0.8, 0.0, 0.0},
                        {0.0, 0.1, 0.9, 0.0, 0.0},
                        {0.1, 0.1, 0.6, 0.1, 0.1},
                        {0.0, 0.0, 0.8, 0.2, 0.0},
                        {0.0, 0.0, 0.9, 0.0, 0.1}};
    double B5[5][5] = { {0.9, 0.0, 0.1, 0.0, 0.0},
                        {0.0, 0.1, 0.9, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 1.0},
                        {0.1, 0.0, 0.0, 0.9, 0.0},
                        {0.0, 0.9, 0.0, 0.1, 0.0}};
    double C5[5] = {0.4, 0.4, 0.0, 0.1, 0.1};
    models[5] = HMM_create(5, 5, "Example Model 5");

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            models[5]->A[i][j] = A5[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            models[5]->B[i][j] = B5[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        models[5]->C[i] = C5[i];
    }



    double A6[6][6] = { {0.0, 0.0, 0.9, 0.1, 0.0, 0.0},
                        {0.0, 0.5, 0.0, 0.5, 0.0, 0.0},
                        {0.0, 0.0, 0.4, 0.3, 0.3, 0.0},
                        {0.0, 0.0, 0.2, 0.1, 0.7, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
    double B6[6][5] = { {1.0, 0.0, 0.0, 0.0, 0.0},
                        {0.8, 0.1, 0.1, 0.0, 0.0},
                        {0.0, 0.9, 0.1, 0.0, 0.0},
                        {0.0, 0.5, 0.5, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 1.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 1.0}};
    double C6[6] = {0.7, 0.3, 0.0, 0.0, 0.0, 0.0};
    models[6] = HMM_create(6, 5, "Example Model 6");
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            models[6]->A[i][j] = A6[i][j];
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 5; j++) {
            models[6]->B[i][j] = B6[i][j];
        }
    }

    for (int i = 0; i < 6; i++) {
        models[6]->C[i] = C6[i];
    }

    double A7[5][5] = { {0.4, 0.5, 0.0, 0.0, 0.1},
                        {0.1, 0.4, 0.5, 0.0, 0.0},
                        {0.0, 0.1, 0.4, 0.5, 0.0},
                        {0.0, 0.0, 0.1, 0.4, 0.5},
                        {0.5, 0.0, 0.0, 0.1, 0.4}};
    double B7[5][5] = { {0.4, 0.3, 0.0, 0.0, 0.3},
                        {0.3, 0.4, 0.3, 0.0, 0.0},
                        {0.0, 0.3, 0.4, 0.3, 0.0},
                        {0.0, 0.0, 0.3, 0.4, 0.3},
                        {0.3, 0.0, 0.0, 0.3, 0.4}};
    double C7[5] = {0.7, 0.3, 0.0, 0.0, 0.0};

    models[7] = HMM_create(5, 5, "Example Model 7");

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            models[7]->A[i][j] = A7[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            models[7]->B[i][j] = B7[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        models[7]->C[i] = C7[i];
    }

    double A8[5][5] = { {0.2, 0.4, 0.4, 0.0, 0.0},
                        {0.0, 0.0, 0.5, 0.5, 0.0},
                        {0.0, 0.0, 0.0, 0.4, 0.6},
                        {0.0, 0.0, 0.0, 1.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 1.0}};
    double B8[5][6] = { {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.3, 0.4, 0.3, 0.0, 0.0},
                        {0.0, 0.6, 0.2, 0.2, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.8, 0.2},
                        {0.0, 0.0, 0.0, 0.0, 0.2, 0.8}};
    double C8[5] = {0.8, 0.2, 0.0, 0.0, 0.0};

    models[8] = HMM_create(5, 6, "Example Model 8");

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            models[8]->A[i][j] = A8[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 6; j++) {
            models[8]->B[i][j] = B8[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        models[8]->C[i] = C8[i];
    }


    double A9[5][5] = { {0.1, 0.45, 0.45, 0.0, 0.0},
                        {0.0, 0.0, 0.1, 0.8, 0.1},
                        {0.0, 0.1, 0.0, 0.0, 0.9},
                        {0.0, 0.0, 0.0, 1.0, 0.0},
                        {0.0, 0.1, 0.8, 0.0, 0.1}};
    double B9[5][6] = { {0.3, 0.3, 0.4, 0.0, 0.0, 0.0},
                        {0.0, 0.3, 0.3, 0.4, 0.0, 0.0},
                        {0.0, 0.0, 0.3, 0.3, 0.4, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.4, 0.6},
                        {0.2, 0.1, 0.0, 0.0, 0.0, 0.7}};
    double C9[5] = {1.0, 0.0, 0.0, 0.0, 0.0};

    models[9] = HMM_create(5, 6, "Example Model 9");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            models[9]->A[i][j] = A9[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 6; j++) {
            models[9]->B[i][j] = B9[i][j];
        }
    }

    for (int i = 0; i < 5; i++) {
        models[9]->C[i] = C9[i];
    }

    return models;
}



/**
 * @brief Initilizes ten example HMM. The models are designed by ChatGPT. 
 * 
 * @return HMM** 
 */
HMM** initialize_example_models_old() {
    // Initilizes ten example HMM. The models are designed by ChatGPT. 
    
    HMM **models = (HMM **)malloc(NUM_MODELS * sizeof(HMM *)); // Allocate memory for 10 HMM models


    double A0[3][3] = {{0.7, 0.2, 0.1}, {0.3, 0.5, 0.2}, {0.2, 0.3, 0.5}};
    double B0[3][4] = {{0.2, 0.3, 0.3, 0.2}, {0.5, 0.2, 0.1, 0.2}, {0.3, 0.2, 0.3, 0.2}};
    double C0[3] = {0.4, 0.3, 0.3};
    models[0] = HMM_create(3, 4, "Example Model 0");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            models[0]->A[i][j] = A0[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            models[0]->B[i][j] = B0[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        models[0]->C[i] = C0[i];
    }



    double A1[2][2] = {{0.7, 0.3}, {0.4, 0.6}};
    double B1[2][3] = {{0.1, 0.4, 0.5}, {0.7, 0.2, 0.1}};
    double C1[2] = {0.5, 0.5};
    models[1] = HMM_create(2, 3, "Example Model 1");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            models[1]->A[i][j] = A1[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            models[1]->B[i][j] = B1[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        models[1]->C[i] = C1[i];
    }

    double A2[3][3] = {{0.6, 0.2, 0.2}, {0.3, 0.4, 0.3}, {0.2, 0.5, 0.3}};
    double B2[3][4] = {{0.1, 0.3, 0.4, 0.2}, {0.6, 0.1, 0.1, 0.2}, {0.2, 0.2, 0.2, 0.4}};
    double C2[3] = {0.4, 0.3, 0.3};
    models[2] = HMM_create(3, 4, "Example Model 2");

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            models[2]->A[i][j] = A2[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            models[2]->B[i][j] = B2[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        models[2]->C[i] = C2[i];
    }


    double A3[2][2] = {{0.6, 0.4}, {0.3, 0.7}};
    double B3[2][3] = {{0.2, 0.5, 0.3}, {0.7, 0.1, 0.2}};
    double C3[2] = {0.5, 0.5};
    models[3] = HMM_create(2, 3, "Example Model 3");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            models[3]->A[i][j] = A3[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            models[3]->B[i][j] = B3[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        models[3]->C[i] = C3[i];
    }


    double A4[3][3] = {{0.8, 0.1, 0.1}, {0.2, 0.7, 0.1}, {0.3, 0.2, 0.5}};
    double B4[3][4] = {{0.2, 0.5, 0.2, 0.1}, {0.3, 0.2, 0.1, 0.4}, {0.5, 0.3, 0.1, 0.1}};
    double C4[3] = {0.4, 0.3, 0.3};
    models[4] = HMM_create(3, 4, "Example Model 4");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            models[4]->A[i][j] = A4[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            models[4]->B[i][j] = B4[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        models[4]->C[i] = C4[i];
    }


    double A5[2][2] = {{0.9, 0.1}, {0.2, 0.8}};
    double B5[2][3] = {{0.4, 0.3, 0.3}, {0.1, 0.2, 0.7}};
    double C5[2] = {0.6, 0.4};
    models[5] = HMM_create(2, 3, "Example Model 5");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            models[5]->A[i][j] = A5[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            models[5]->B[i][j] = B5[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        models[5]->C[i] = C5[i];
    }


    double A6[4][4] = {{0.5, 0.1, 0.2, 0.2}, {0.2, 0.3, 0.2, 0.3}, {0.1, 0.2, 0.4, 0.3}, {0.3, 0.2, 0.1, 0.4}};
    double B6[4][5] = {{0.1, 0.2, 0.3, 0.2, 0.2}, {0.3, 0.1, 0.2, 0.2, 0.2}, {0.2, 0.2, 0.1, 0.3, 0.2}, {0.2, 0.1, 0.2, 0.2, 0.3}};
    double C6[4] = {0.3, 0.2, 0.2, 0.3};
    models[6] = HMM_create(4, 5, "Example Model 6");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            models[6]->A[i][j] = A6[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            models[6]->B[i][j] = B6[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        models[6]->C[i] = C6[i];
    }


    double A7[3][3] = {{0.8, 0.1, 0.1}, {0.2, 0.7, 0.1}, {0.1, 0.2, 0.7}};
    double B7[3][4] = {{0.2, 0.3, 0.2, 0.3}, {0.3, 0.2, 0.3, 0.2}, {0.1, 0.4, 0.1, 0.4}};
    double C7[3] = {0.4, 0.3, 0.3};
    models[7] = HMM_create(3, 4, "Example Model 7");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            models[7]->A[i][j] = A7[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            models[7]->B[i][j] = B7[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        models[7]->C[i] = C7[i];
    }


    double A8[2][2] = {{0.6, 0.4}, {0.3, 0.7}};
    double B8[2][3] = {{0.3, 0.2, 0.5}, {0.6, 0.1, 0.3}};
    double C8[2] = {0.5, 0.5};
    models[8] = HMM_create(2, 3, "Example Model 8");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            models[8]->A[i][j] = A8[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            models[8]->B[i][j] = B8[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        models[8]->C[i] = C8[i];
    }


    double A9[4][4] = {{0.4, 0.3, 0.2, 0.1}, {0.2, 0.4, 0.2, 0.2}, {0.1, 0.2, 0.5, 0.2}, {0.3, 0.1, 0.3, 0.3}};
    double B9[4][5] = {{0.2, 0.3, 0.2, 0.2, 0.1}, {0.3, 0.1, 0.3, 0.2, 0.1}, {0.1, 0.3, 0.1, 0.2, 0.3}, {0.2, 0.2, 0.2, 0.3, 0.1}};
    double C9[4] = {0.3, 0.2, 0.2, 0.3};
    models[9] = HMM_create(4, 5, "Example Model 9");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            models[9]->A[i][j] = A9[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            models[9]->B[i][j] = B9[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        models[9]->C[i] = C9[i];
    }

    double A10[2][2] = {{0.5, 0.5}, {0.5, 0.5}};
    double B10[2][2] = {{0.5, 0.5}, {0.5, 0.5}};
    double C10[2] = {0.5, 0.5};
    models[10] = HMM_create(2, 2, "Example Model 10");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            models[10]->A[i][j] = A10[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            models[10]->B[i][j] = B10[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        models[10]->C[i] = C10[i];
    }


    return models;
}


void free_example_models(HMM **models) {
    if (models != NULL) {
        for (int i = 0; i < NUM_MODELS; ++i) {
            if (models[i] != NULL) {
                HMM_destroy(models[i]);
            }
        }
        free(models);
    }
}