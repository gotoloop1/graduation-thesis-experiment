#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

void read_matrix(
    const char *restrict path,
    raw_matrix *restrict m
) {
    FILE* f = fopen(path, "r");
    char buffer[100];
    bool first = true;
    int count = 0;
    while(1) {
        char *_ = fgets(buffer, 100, f);
        if(buffer[0] != '%') {
            if(first) {
                int a, b, c;
                sscanf(buffer, "%d %d %d", &a, &b, &c);
                if(a != b) {
                    printf("not square matrix\n");
                    exit(1);
                }
                m->n = a;
                m->nz = c * 2 - m->n;
                m->row = (int*)malloc(sizeof(int) * m->nz);
                m->column = (int*)malloc(sizeof(int) * m->nz);
                m->A = (double*)malloc(sizeof(double) * m->nz);
                first = false;
            } else {
                int a, b;
                double c;
                sscanf(buffer, "%d %d %lf", &a, &b, &c);
                m->row[count] = a - 1;
                m->column[count] = b - 1;
                m->A[count] = c;
                count++;
                if(a != b) {
                    m->row[count] = b - 1;
                    m->column[count] = a - 1;
                    m->A[count] = c;
                    count++;
                }
                if(count == m->nz) {
                    break;
                }
            }
        }
    }
    fclose(f);
}
