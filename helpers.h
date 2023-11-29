#ifndef HELPERS  
#define HELPERS


#include <stdio.h>
#include "cudd.h"


void write_dd (DdManager *gbm, DdNode *dd, char* filename );
void print_dd (DdManager *gbm, DdNode *dd, int n, int pr );

struct DoubleArray
{
    size_t n;
    double *data;
};

struct DoubleMatrix
{
    size_t n;
    struct DoubleArray *data;
};

double** allocate_matrix(int rows, int cols, double initialValue);
double*** allocate_3D_matrix(int depth, int rows, int cols, double initialValue);
double log_sum_exp(double a, double b);
void save_list(const int *list, int size, const char *filename);
int **read_list(const char *filename);


#endif