#ifndef HELPERS  
#define HELPERS


#include <stdio.h>
#include "cudd.h"
#include <stdbool.h>


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

#endif