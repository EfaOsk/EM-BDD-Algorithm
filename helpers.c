#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "cudd.h"
#include <math.h>
#include <errno.h>

#include <sys/resource.h>

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
    Cudd_DumpDot(gbm, 4, ddnodearray, NULL, NULL, outfile); // dump the function to .dot file
    free(ddnodearray);
    fclose (outfile); // close the file */
}

/**
 * @brief Generates a random sequence of length T with (whole) numbers in range 0-M
 * 
 * @param M int
 * @param T int
 * @return int* 
 */
int* generate_sequence(int M, int T) 
{
    int* sequence = malloc(T * sizeof(int));
    if (sequence == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < T; i++) {
        sequence[i] = rand() % M;
    }

    return sequence;
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


double**** allocate_4D_matrix(int height, int depth, int rows, int cols, double initialValue) {
    double**** matrix = (double****)malloc(height * sizeof(double***));
    if (matrix == NULL) {
        return NULL;
    }

    for (int h = 0; h < height; h++) {
        matrix[h] = allocate_3D_matrix(depth, rows, cols, initialValue);
        if (matrix[h] == NULL) {
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}


void free_matrix(double **matrix, int rows) {
    if (matrix != NULL) {
        for (int i = 0; i < rows; i++) {
            if (matrix[i] != NULL) {
                free(matrix[i]); 
            }
        }
        free(matrix);
    }
}


void free_3D_matrix(double ***matrix, int rows, int cols) {
    if (matrix != NULL) {
        for (int i = 0; i < rows; i++) {
            if (matrix[i] != NULL) {
                for (int j = 0; j < cols; j++) {
                    if (matrix[i][j] != NULL) {
                        free(matrix[i][j]);
                    }
                }
                free(matrix[i]);
            }
        }
        free(matrix);
    }
}


void free_4D_matrix(double**** matrix, int height, int depth, int rows) {
    for (int h = 0; h < height; h++) {
        free_3D_matrix(matrix[h], depth, rows);
    }
    free(matrix);
}


void reset_matrix(double **matrix, int rows, int cols, double value) {
    if (matrix != NULL) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = value;
            }
        }
    }
}

void reset_3D_matrix(double ***matrix, int depth, int rows, int cols, double value) {
    for (int d = 0; d < depth; d++) {        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[d][i][j] = value;
            }
        }
    }
}


void reset_4D_matrix(double ****matrix, int height, int depth, int rows, int cols, double value) {
    for (int h = 0; h < height; h++) {
        for (int d = 0; d < depth; d++) {        
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    matrix[h][d][i][j] = value;
                }
            }
        }
    }
}



double log_sum_exp(double a, double b) {
    double max_val = (a > b) ? a : b;
    if (a == -INFINITY) {
        return b;
    } 
    if (b == -INFINITY) {
        return a;
    }
    return max_val + log(exp(a - max_val) + exp(b - max_val));
}


// Function to save a list of integers to a file
void save_list(const int *list, int size, const char *filename) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file for writing");
        return;
    }

    for (int i = 0; i < size; i++) {
        fprintf(file, "%d", list[i]);
        if (i < size - 1) {
            fprintf(file, ", ");
        }
    }
    fprintf(file, "\n");

    fclose(file);
}


// Function to read a list of lists of integers from a file
int **read_list(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file for reading");
        return NULL;
    }

    char line[1024];
    int **lists = NULL;
    int size = 0;

    while (fgets(line, sizeof(line), file)) {
        int *sublist = NULL;
        int subsize = 0;
        char *token = strtok(line, ", ");
        while (token != NULL) {
            sublist = realloc(sublist, (subsize + 1) * sizeof(int));
            sublist[subsize++] = atoi(token);
            token = strtok(NULL, ", ");
        }

        lists = realloc(lists, (size + 1) * sizeof(int *));
        lists[size] = sublist;
        size++;
    }

    fclose(file);
    return lists;
}


long get_mem_usage() {
    struct rusage myusage;
    getrusage(RUSAGE_SELF, &myusage);
    return myusage.ru_maxrss;
}

