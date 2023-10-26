#include "HMM.h"
#include "helpers.h"


#include "exampleHMM.h"
#include "stdlib.h"

/**
 * @brief Initilizes ten example HMM. The models are designed by ChatGPT. 
 * 
 * @return HMM** 
 */
HMM** initialize_example_models() {
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
    double C10[2] = {1, 0};
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
