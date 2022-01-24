#ifndef MATRIX_INCLUDED
#define MATRIX_INCLUDED

typedef struct {
    int n;
    int nz;
    int *restrict row;
    int *restrict column;
    double *restrict A;
} raw_matrix;

void read_matrix(
    const char *restrict path,
    raw_matrix *restrict m
);

#endif
