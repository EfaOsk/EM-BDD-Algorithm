#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cudd.h"
#include "HMM.h"
#include "BDD.h"
#include "helpers.h"



double get_theta(const HMM *hmm, int x, int i, int j) { 
    if (x == 0) {
        return hmm->B[i][j];
    } else if (x == 1) {
        return hmm->C[j];
    } else if (x == 2) {
        return hmm->A[i][j];
    } else {
        perror("Error: Invalid value of x");
        return -1.0; 
    }
}


double get_sigma(const HMM *hmm, int x, int i, int j) { 
    double sum = 0.0;
    if (x == 0) {
        for (int j0 = j; j0 < hmm->M; j0++) {
            sum += hmm->B[i][j0];
        }
    } else if (x == 1) {
        for (int j0 = j; j0 < hmm->N; j0++){
            sum += hmm->C[j0];
        }
    } else if (x == 2) {
        for (int j0 = j; j0 < hmm->N; j0++) {
            sum += hmm->A[i][j0];
        }
    } else {
        perror("Error: Invalid value of x");
        return -1.0; 
    }
    return sum;
}


double get_prob_encoded(DdManager* manager, const HMM *hmm, DdNode *n, int b, int **lookup_table_variables) { 
    int index = Cudd_NodeReadIndex(n);
    int x = lookup_table_variables[index][0];
    int i = lookup_table_variables[index][1];
    // int t = lookup_table_variables[index][2];
    int j = lookup_table_variables[index][3];

    if (b == 0) {
        return get_sigma(hmm, x, i, j+1) / get_sigma(hmm, x, i, j);
    } else if (b == 1) {
        return get_theta(hmm, x, i, j) / get_sigma(hmm, x, i, j);
    } else {
        perror("Error: Invalid value of b");
        return -1.0; 
    }
}

NodeDataList** nodeData;

NodeDataNode* FindTargetNodeAtLevel(DdManager* manager, int targetLevel, DdNode* targetNode) {
    if (targetLevel == -1) {
        int numVars = Cudd_ReadSize(manager);
        targetLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(targetNode));
        if (targetLevel > numVars) {
            targetLevel = numVars;
        }
    }
    NodeDataNode* current = nodeData[targetLevel]->head;
    while (current != NULL) {
        if (current->node == targetNode) {
            return current; // Return the NodeDataNode containing the target node
        }
        current = current->next;
    }
    printf("Notfound");
    return NULL; // Node not found at the specified level
}



/**
 * @brief Backward procedure where
 * 
 *      /Beta(x, n) = probability that paths logically reach the terminal node x from node n
 * 
 *  for level in vars: // from lowest level to highest
 *      for node n in level:
 *          B_1[n] = 0.5* B_1[Child(True)] + 0.5* B_1[Child(False)]
 * 
 * @param manager 
 * @param node 
 * @param M 
 * @return double
 */
double Backward(DdManager* manager, DdNode* node, const HMM *hmm, int** lookup_table_variables) {
    NodeDataNode *targetNodeData = FindTargetNodeAtLevel(manager, -1, node);
    
    if (targetNodeData->backward >= 0) {
        return targetNodeData->backward; 
    }


    DdNode* high = Cudd_T(node);
    DdNode* low = Cudd_E(node);
    
    double prob_high = get_prob_encoded(manager, hmm, node, 1, lookup_table_variables) * Backward(manager, high, hmm, lookup_table_variables);
    double prob_low;

    if (!Cudd_IsComplement(low)) {
        prob_low = get_prob_encoded(manager, hmm, node, 0, lookup_table_variables) * Backward(manager, low, hmm, lookup_table_variables);
    } else {
        // Adjusting for negative (complemented) edge
        prob_low = get_prob_encoded(manager, hmm, node, 0, lookup_table_variables) * (1.0 - Backward(manager, Cudd_Regular(low), hmm, lookup_table_variables));
    }
    
    // Calculate the probability for the current node
    double prob = prob_low + prob_high;
    targetNodeData->backward = prob; // Store in the lookup table

    return prob;
}


/**
 * Calculate forward values for nodes in a Binary Decision Diagram (BDD).
 *
 * @param manager   The CUDD manager for the BDD operations.
 * @param F_all     An array of BDD nodes representing the roots of the BDDs.
 * @param hmm       A pointer to a Hidden Markov Model (HMM) or a similar probabilistic model.
 * @param num_roots The number of BDD roots in the 'F_all' array.
 */
void CalculateForward(DdManager* manager, DdNode** F_seq, const HMM *hmm, int T, int num_sequences, int** lookup_table_variables) {
    int numVars = Cudd_ReadSize(manager); // Cudd_ReadNodeCount(manager)+(hmm->N)*T+2;
 
    for (int r = 0; r < num_sequences; r++) {
        DdNode *targetNode = F_seq[r];
        int isNegated = Cudd_IsComplement(targetNode);
        if (isNegated) {
            targetNode = Cudd_Regular(targetNode);
        }
        int level = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(targetNode));
        // Find target node
        NodeDataNode* targetNodeData = FindTargetNodeAtLevel(manager, level, targetNode);
        if (targetNodeData == NULL) {
            // raise error, nod not found
            printf("ERROR!");
            return;
        }
        double backwardVal = Backward(manager, targetNode, hmm, lookup_table_variables);
        if (!isNegated) {
            // even number of comple edges
            targetNodeData->forward[1] += 1 / (1-backwardVal);
        } else {
            // odd number of comple edges
            targetNodeData->forward[0] += 1 / (1-backwardVal);
        }
    }


    for (int level = 0 ; level < numVars; level++) {
        // int index = Cudd_ReadInvPerm(manager, level);
        NodeDataNode* targetNode = nodeData[level]->head;
        while (targetNode != NULL) {
            // Find children:
            DdNode *lowChild = Cudd_E(targetNode->node);
            DdNode *highChild = Cudd_T(targetNode->node);
            
            int lowLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(lowChild));
            if (lowLevel > numVars) {
                lowLevel = numVars;
            }
            int highLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(highChild));
            if (highLevel > numVars) {
                highLevel = numVars;
            }
            NodeDataNode* lowNode = FindTargetNodeAtLevel(manager, lowLevel, Cudd_Regular(lowChild));
            NodeDataNode* highNode = FindTargetNodeAtLevel(manager, highLevel, highChild);

            double ProbLowEdge = get_prob_encoded(manager, hmm, targetNode->node, 0, lookup_table_variables);
            double ProbHighEdge = get_prob_encoded(manager, hmm, targetNode->node, 1, lookup_table_variables);

            highNode->forward[0] += targetNode->forward[0]* ProbHighEdge; 
            highNode->forward[1] += targetNode->forward[1]* ProbHighEdge; 

            if (!Cudd_IsComplement(lowChild)) {
                lowNode->forward[0] += targetNode->forward[0]* ProbLowEdge; 
                lowNode->forward[1] += targetNode->forward[1]* ProbLowEdge; 
            } else {
                lowNode->forward[0] += targetNode->forward[1]* ProbLowEdge; 
                lowNode->forward[1] += targetNode->forward[0]* ProbLowEdge; 
            }

            // next node
            targetNode = targetNode->next;
        }
    }
}


void InitNodeData(DdManager* manager, DdNode** F_seq, int T, int num_sequences) {
    int numVars = Cudd_ReadSize(manager); 
    nodeData = (NodeDataList **)malloc((numVars + 1) * sizeof(NodeDataList*));

    if (nodeData == NULL) {
        // Handle memory allocation failure
        perror("Memory allocation failed");
    }
    for (int i = 0; i < numVars + 1; i++) {
        nodeData[i] = (NodeDataList *)malloc(sizeof(NodeDataList));
        nodeData[i]->head = NULL;
        nodeData[i]->tail = NULL;
    }

    DdNode *node;
    DdGen *gen;


    for (int r = 0; r < num_sequences; r++){
        Cudd_ForeachNode(manager, F_seq[r], gen, node) {
            if (node == Cudd_ReadLogicZero(manager) || node == Cudd_Not(Cudd_ReadLogicZero(manager))) {
                // Terminal node
            } else {
                // Non-terminal nod
                // Get level
                int level = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(node));

                // Check if the node is already in nodeData[level]
                NodeDataNode* current = nodeData[level]->head;
                int nodeExists = 0;

                while (current != NULL) {
                    if (current->node == node) {
                        nodeExists = 1;
                        break;
                    }
                    current = current->next;
                }

                if (!nodeExists) {
                    NodeDataNode* newNode = (NodeDataNode*)malloc(sizeof(NodeDataNode));
                    newNode->node = node;
                    newNode->forward[0] = 0.0; // Inizilize forward with 0
                    newNode->forward[1] = 0.0; // Inizilize forward with 0
                    newNode->backward = -1.0; // Inizilize backward with 0
                    
                    
                    
                    newNode->next = NULL;

                    // Append the new node to the linked list for the level
                    if (nodeData[level]->tail == NULL) {
                        nodeData[level]->head = newNode;
                        nodeData[level]->tail = newNode;
                    } else {
                        nodeData[level]->tail->next = newNode;
                        nodeData[level]->tail = newNode;
                    }
                }
            }
        }
    }

    NodeDataNode* trueTerminal = (NodeDataNode*)malloc(sizeof(NodeDataNode));
    trueTerminal->node = Cudd_Not(Cudd_ReadLogicZero(manager));
    trueTerminal->forward[0] = 0.0;
    trueTerminal->forward[1] = 0.0;
    trueTerminal->backward = 1.0;
    trueTerminal->next = NULL;
    nodeData[numVars]->head = trueTerminal;
    nodeData[numVars]->tail = trueTerminal;
    
}


void CleanNodeData(int numVars) {
    for (int level = 0; level <= numVars; level++) {
        NodeDataNode* node = nodeData[level]->head;
        while (node != NULL) {
            node->forward[0] = 0.0; // Inizilize forward with 0
            node->forward[1] = 0.0; // Inizilize forward with 0
            node->backward = -1.0; // Inizilize backward with 0
            // next node
            node = node->next;
        }
    }
    nodeData[numVars]->head->backward = 1.0; // Inizilize forward with 0
}


void FreeNodeData(int numVars) {
    // Clean up the allocated memory
    for (int i = 0; i < numVars + 1; i++) {
        NodeDataNode* current = nodeData[i]->head;
        while (current != NULL) {
            NodeDataNode* temp = current;
            current = current->next;
            free(temp);
        }
        free(nodeData[i]);
    }
    free(nodeData);
}


void computeConditionalExpectations(DdManager *manager, const HMM *hmm, int T, double ***eta,  double ***gamma,  double *D, int **lookup_table_variables) {

    int N = hmm->N;
    int NX = (hmm->N > hmm->M) ? hmm->N : hmm->M;
    int numVars = Cudd_ReadSize(manager);

    double e0, e1, temp;

    // Reset gamma and eta
    for (int x = 0; x < 3; x++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < NX; j++) {
                gamma[x][i][j] = 0.0;
                eta[x][i][j] = 0.0;
            }
        }
    }

    // Reset D
    for (int i = 0; i < numVars+1; i++) {
        D[i] = 0.0;
    }

    for (int level = 0; level < numVars; level++) {
        NodeDataNode* node = nodeData[level]->head;
        int index = Cudd_ReadInvPerm(manager, level);
        int x = lookup_table_variables[index][0];
        int i = lookup_table_variables[index][1];
        int t = lookup_table_variables[index][2];
        int j = lookup_table_variables[index][3];
        while (node != NULL) {
            DdNode *lowChild = Cudd_E(node->node);
            DdNode *highChild = Cudd_T(node->node);
                
            int lowLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(Cudd_Regular(lowChild)));
            if (lowLevel > numVars) {
                lowLevel = numVars;
            }
            int highLevel = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(highChild));
            if (highLevel > numVars) {
                highLevel = numVars;
            }
            NodeDataNode *lowChildData = FindTargetNodeAtLevel(manager, lowLevel, Cudd_Regular(lowChild));
            NodeDataNode *highChildData = FindTargetNodeAtLevel(manager, highLevel, highChild);
            double PrLow = get_prob_encoded(manager, hmm, node->node, 0, lookup_table_variables);
            double PrHigh = get_prob_encoded(manager, hmm, node->node, 1, lookup_table_variables);

            e1 = node->forward[0]*PrHigh*(1-highChildData->backward) + node->forward[1]*PrHigh*highChildData->backward;
            
            if (!Cudd_IsComplement(lowChild)) {
                e0 = node->forward[0]*PrLow*(1-lowChildData->backward) + node->forward[1]*PrLow*lowChildData->backward;
            } else {
                e0 = node->forward[0]*PrLow*(lowChildData->backward) + node->forward[1]*PrLow*(1-lowChildData->backward);
            }
            eta[x][i][j] += e1;

            gamma[x][i][j+1] += e0;
            gamma[x][i][j] -= e0 + e1;

            
            D[index+1] += e0 + e1;
            int highIndex = Cudd_ReadInvPerm(manager, highLevel);
            if (highIndex == -1) {  highIndex = numVars;    }
            D[highIndex] -= e1;
            int lowIndex = Cudd_ReadInvPerm(manager, lowLevel);
            if (lowIndex == -1) {  lowIndex = numVars;    }
            D[lowIndex] -= e0;

            // next node
            node = node->next;
        }
    }

    temp = 0;
    for (int index = 0; index < numVars; index++) {
        int x = lookup_table_variables[index][0];
        int i = lookup_table_variables[index][1];
        int t = lookup_table_variables[index][2];
        int j = lookup_table_variables[index][3];
        temp += D[index];
        if (j == 0) {
            gamma[x][i][0] = temp;
        }
    }

    for (int x = 0; x < 3; x++) { 
        int tempN = N;
        if (x == 1) { tempN = 1; }
        for (int i = 0; i < tempN; i++) {
            temp = 0;
            if (x > 0) { NX = N; } 
            else { NX = hmm->M; }
            for (int j = 0; j < NX; j++) {
                temp += gamma[x][i][j] / get_sigma(hmm, x, i, j); // TODO: / sigma[mu(i)][j];
                eta[x][i][j] += temp * get_theta(hmm, x, i, j);  // TODO: * thetda[mu(i)][j];
            }
        }
    }

}


HMM* BDD_update(HMM *hmm, double ***eta) {
    HMM *new_hmm = HMM_create(hmm->N, hmm->M, "Updated model");
    double sum;
    double min_p_f = 0.00001;
    // update b
    for (int u = 0; u < hmm->N; u++) {
        sum = min_p_f*hmm->M;
        // ToDo: remove this quick fix, handles when obs is in alphabet, but not data to learn
        for (int o = 0; o < hmm->M; o++) {
            if (eta[0][u][o]<=0){
                eta[0][u][o] = 0.0;
            }
            sum += eta[0][u][o];
        }
        for (int o = 0; o < hmm->M; o++) {
            new_hmm->B[u][o] = (eta[0][u][o]+min_p_f) / sum;
        }
    }
    // update pi
    sum = hmm->N*min_p_f;
    for (int v = 0; v < hmm->N; v++) {
        sum += eta[1][0][v];
    }
    for (int v = 0; v < hmm->N; v++) {
        new_hmm->C[v] = (eta[1][0][v]+min_p_f) / sum;
    }
    // update a
    for (int u = 0; u < hmm->N; u++) {
        sum = hmm->N*min_p_f;
        for (int v = 0; v < hmm->N; v++) {
            sum += eta[2][u][v];
        }
        for (int v = 0; v < hmm->N; v++) {
            new_hmm->A[u][v] = (eta[2][u][v]+min_p_f) / sum;
        }
    }
    return new_hmm;
}


/**
 * @brief This function builds an HMM that "best" fits a given observation sequence O with the given number of states and observations.
 * 
 *      Uses the cour structure of the Baum-Welch algorithm, and BDDs, in order to improve the time and memory complexity
 *       
 * @param hypothesis_hmm  Pointer to an HMM structure, which is the initial hypothesis for the HMM.
 * @param T               The length of the observation sequence.
 * @param observations    Observation sequence, an array of integers of length T.
 * @param epsilon         The convergence threshold; when the change in log likelihood is less than or equal to this value, learning stops.
 * @return HMM*           Pointer to the learned HMM structure.
 */
HMM* EMBDD_learn( HMM *hypothesis_hmm, int num_sequences, int **observations, int T, double epsilon, const char *logs_folder, const char *result_file) {
    /*

    EM on BDD
        (1) Build (S)BDD
        (1.5) Ordered encoding
        
        (2) initilize M (= some random HMM) 

        (repeat steps 3-5 until converged)

        (3) E-step
            (a) Backward
            (b) Forward
            (c) Conditional expectations
        
        (4) M-step
            (a) update M
        
        (5) Calculate the log-likelyhood of M

        (6) retrun M

    */
    
    int N = hypothesis_hmm->N;
    int M = hypothesis_hmm->M;
    HMM *model = HMM_create(N, M, "model");

    HMM_copy(model, hypothesis_hmm);

    char log_filename[256];
    char model_filename[256];
    // Construct the log filename
    sprintf(log_filename, "%s/log.txt", logs_folder);

    // Open the log file
    FILE *log_file = fopen(log_filename, "w");
    if (log_file == NULL) {
        perror("Error opening log file");
        return NULL;
    }

    sprintf(model_filename, "%s/models", logs_folder);
    mkdir(model_filename, 0777);


    // Step 1 build (S)BDD
    DdManager *manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0);  
    Cudd_AutodynEnable(manager, CUDD_REORDER_SAME);
    
    int **lookup_table_variables;
    lookup_table_variables = (int **)malloc((N*T*(M-1)+(N-1)+(N*(T-1)*(N-1))) * sizeof(int *));
    for (int id = 0; id <(N*T*(M-1)+(N-1)+(N*(T-1)*(N-1))); id++) {
        (lookup_table_variables)[id] = (int *)malloc(4 * sizeof(int));
    }

    printf("Building the BDD ...\n");
    DdNode **F_obs = build_F_seq(manager, N, M, num_sequences, T, observations, lookup_table_variables);

    int numVars = Cudd_ReadSize(manager);


    char BDDfilename[526];
    sprintf(BDDfilename, "%s/BDD.dot", logs_folder); // Write .dot filename to a string
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(BDDfilename,"w");
    Cudd_DumpDot(manager, (1), F_obs, NULL, NULL, outfile);
    fclose(outfile);


    // Step 2: initilize M (= some random HMM) 
        // ToDo currently input


    double prob_priv, prob_original, prob_new;
    prob_original = log_likelihood_forward_multiple(model, observations, num_sequences, T);
    prob_priv = prob_original;
    int converged = 0;
    int iteration = 0;

    int tmp = (N > M) ? N : M;
    double ***eta = allocate_3D_matrix(3, N, tmp, 0.0);
    double ***gamma = allocate_3D_matrix(3, N, tmp, 0.0);
    double *D = malloc((numVars+1) * sizeof(double));

    InitNodeData(manager, F_obs, T, num_sequences);

    printf("Starting EM-BDD learning process...\n");
    while (!converged)
    {
        clock_t start_time = clock();

        // Step 3: E-step
        
        CleanNodeData(numVars);
    
        // Step 3 (a) : Backward

        // Step 3 (b) : Forward
    
        CalculateForward(manager, F_obs, model, T, num_sequences, lookup_table_variables);

        // Step 3 (c) : Conditional Expectations

        computeConditionalExpectations(manager, model, T, eta, gamma, D, lookup_table_variables);
    
        // Step 4: M-step
        // Step 4 (a) : update M

        HMM *new_hmm = BDD_update(model, eta);


        prob_new = log_likelihood_forward_multiple(new_hmm, observations, num_sequences, T);
        
        
        HMM_copy(model, new_hmm); 
        HMM_destroy(new_hmm);

        HMM_validate(model);

        clock_t end_time = clock();
        double iteration_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        fprintf(log_file, "Iteration %d: Log-likelihood = %f, Improvement = %f, Time = %f seconds\n", 
            ++iteration, prob_new, prob_new - prob_priv, iteration_time);
        fflush(log_file); 
        
        sprintf(model_filename, "%s/models/model_%d", logs_folder, iteration);
        HMM_save(model, model_filename); 
        if (prob_new <= prob_priv+epsilon) {
            converged = 1;
            printf("Convergence achieved after %d iterations.\n", iteration);
        }
        prob_priv = prob_new;
    }

    fclose(log_file);

    // Open the result file in append mode
    FILE *result_fp = fopen(result_file, "a");
    if (result_fp == NULL) {
        perror("Error opening result file");
        HMM_destroy(model);
        fclose(log_file);
    }

    fprintf(result_fp, "%d, %f, %f\n", iteration, prob_new, prob_new - prob_original);

    // Clean up
    for (int x = 0; x < 3; x++) {
        for (int i = 0; i < N; i++) {
            free(gamma[x][i]);
            free(eta[x][i]);
        }
        free(gamma[x]);
        free(eta[x]);
    }
    free(gamma);
    free(eta);

    free(D);


    free(F_obs);
    FreeNodeData(numVars);
    free_lookup_table_variables(numVars, lookup_table_variables);
    Cudd_Quit(manager);

    // Step 6: Return the learned model

    return model;
}
