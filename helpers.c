#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "cudd.h"
#include <math.h>


/**
 * Print a dd summary
 * pr = 0 : prints nothing
 * pr = 1 : prints counts of nodes and minterms
 * pr = 2 : prints counts + disjoint sum of product
 * pr = 3 : prints counts + list of nodes
 * pr > 3 : prints counts + disjoint sum of product + list of nodes
 * @param the dd node
 */
void print_dd (DdManager *gbm, DdNode *dd, int n, int pr )
{
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(gbm)); /*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager vars: %d | ", Cudd_ReadSize(gbm) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(gbm) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(gbm) ); /*Returns the memory in use by the manager measured in bytes*/
    // Cudd_PrintDebug(gbm, dd, n, pr);  // Prints to the standard output a DD and its statistics: number of nodes, number of leaves, number of minterms.
}


/**
 * Writes a dot file representing the argument DDs
 * @param the node object
 */
void write_dd (DdManager *gbm, DdNode *dd, char* filename)
{
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    DdNode **ddnodearray = (DdNode**)malloc(sizeof(DdNode*)); // initialize the function array
    ddnodearray[0] = dd;
    Cudd_DumpDot(gbm, 1, ddnodearray, NULL, NULL, outfile); // dump the function to .dot file
    free(ddnodearray);
    fclose (outfile); // close the file */
}

/**
 * Writes a dot file representing the argument DDs
 * @param the node object
 */
void write_dds (DdManager *gbm, DdNode **ddnodearray, char* filename)
{
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    // DdNode **ddnodearray = (DdNode**)malloc(sizeof(DdNode*)); // initialize the function array
    // ddnodearray[0] = dd;
    Cudd_DumpDot(gbm, 2, ddnodearray, NULL, NULL, outfile); // dump the function to .dot file
    free(ddnodearray);
    fclose (outfile); // close the file */
}

double** allocate_matrix(int rows, int cols, double initialValue) 
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    if (matrix == NULL) {
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));
        if (matrix[i] == NULL) {
            // Free previously allocated rows
            for (int j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = initialValue;
        }
    }
    return matrix;
}

double*** allocate_3D_matrix(int depth, int rows, int cols, double initialValue) {
    // Allocate memory for the array of pointers to 2D matrices
    double ***matrix = (double ***)malloc(depth * sizeof(double **));
    if (matrix == NULL) {
        return NULL;
    }

    // Allocate each 2D matrix
    for (int d = 0; d < depth; d++) {
        matrix[d] = allocate_matrix(rows, cols, initialValue);
    }

    return matrix;
}


double log_sum_exp(double a, double b) {
    double max_val = (a > b) ? a : b;
    return max_val + log(exp(a - max_val) + exp(b - max_val));
}
