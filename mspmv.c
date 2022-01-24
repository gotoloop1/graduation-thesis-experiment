#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <omp.h>
#include "lib.h"
#include "matrix.h"

#ifndef MSPMV_NUM
#define MSPMV_NUM 32
#endif

typedef struct {
    int n;

    int spmv_block_per_thread;
    int spmv_max_threads;
    int *restrict spmv_block_separator;

    double *restrict D;
    int *restrict index;
    int *restrict item;
    double *restrict A;
} reordered_matrix;

void reorder_matrix(
    const raw_matrix *restrict um,
    const int *restrict order,
    reordered_matrix *restrict m
) {
    m->n = um->n;

    m->spmv_block_per_thread = 25;
    m->spmv_max_threads = omp_get_max_threads();
    const int max_blocks = m->spmv_block_per_thread * m->spmv_max_threads;
    m->spmv_block_separator = (int*)malloc(sizeof(int) * (max_blocks + 1));

    m->D = (double*)malloc(sizeof(double) * um->n);
    m->index = (int*)malloc(sizeof(int) * (um->n + 1));
    m->item = (int*)malloc(sizeof(int) * um->nz);
    m->A = (double*)malloc(sizeof(double) * um->nz);

    int *order_inverse = (int*)malloc(sizeof(int) * um->n);
    for(int i = 0; i < um->n; i++) {
        order_inverse[order[i]] = i;
    }

    int *count = (int*)malloc(sizeof(int) * (um->n + 1));
    for(int i = 0; i <= um->n; i++) {
        count[i] = 0;
    }
    for(int i = 0; i < um->nz; i++) {
        count[order_inverse[um->row[i]] + 1]++;
    }
    for(int i = 0; i < um->n; i++) {
        count[i + 1] += count[i];
    }
    for(int i = um->n; i > 0; i--) {
        count[i] = count[i - 1];
    }
    int *buffer = (int*)malloc(sizeof(int) * um->nz);
    for(int i = 0; i < um->nz; i++) {
        const int p = order_inverse[um->row[i]];
        buffer[count[p + 1]] = i;
        count[p + 1]++;
    }

    int *tmp = (int*)malloc(sizeof(int) * 2 * um->n);
    for(int i = 0; i < um->n; i++) {
        const int b = count[i];
        const int e = count[i + 1];
        for(int j = 0; j < e - b; j++) {
            const int p = buffer[b + j];
            tmp[j * 2] = order_inverse[um->column[p]];
            tmp[j * 2 + 1] = p;
        }
        first_element_sort(e - b, tmp, um->n - 1);
        for(int j = 0; j < e - b; j++) {
            buffer[b + j] = tmp[j * 2 + 1];
        }
    }

    const double num_of_element_per_block = (double)um->nz / max_blocks;
    int block_num = 0;
    m->spmv_block_separator[0] = 0;
    for(int i = 0; i < um->n; i++) {
        if(i + 1 - m->spmv_block_separator[block_num] >= 8 && count[i + 1] >= (block_num + 1) * num_of_element_per_block) {
            block_num++;
            m->spmv_block_separator[block_num] = i + 1;
        }
    }
    if(um->n - m->spmv_block_separator[block_num] < 8) {
        m->spmv_block_separator[block_num] = um->n;
    }
    for(int i = block_num; i < max_blocks; i++) {
        m->spmv_block_separator[i + 1] = um->n;
    }

    #pragma omp parallel default(none), shared(um, m, order_inverse, count, buffer)
    {
        const int thread_num = omp_get_thread_num();
        if(thread_num == 0) {
            m->index[0] = 0;
        }
        for(int bi = 0; bi < m->spmv_block_per_thread; bi++) {
            const int bb = m->spmv_block_separator[bi * m->spmv_max_threads + thread_num];
            const int be = m->spmv_block_separator[bi * m->spmv_max_threads + thread_num + 1];
            for(int i = bb; i < be; i++) {
                const int b = count[i];
                const int e = count[i + 1];
                m->index[i + 1] = count[i];
                bool diag_updated = false;
                for(int j = 0; j < e - b; j++) {
                    const int p = buffer[b + j];
                    const int new_c = order_inverse[um->column[p]];
                    const double a = um->A[p];
                    m->item[m->index[i + 1]] = new_c;
                    m->A[m->index[i + 1]] = a;
                    m->index[i + 1]++;
                    if(i == new_c) {
                        diag_updated = true;
                        m->D[i] = a;
                    }
                }
                if(!diag_updated) {
                    printf("invalid matrix: %d\n", i);
                    exit(1);
                }
            }
        }
    }

    free(order_inverse);
    free(count);
    free(buffer);
    free(tmp);
}

void MSpMV(
    const reordered_matrix *restrict m,
    double *const restrict X[MSPMV_NUM],
    double *restrict Y[MSPMV_NUM]
) {
    #pragma omp parallel default(none), shared(m, X, Y)
    {
        const int thread_num = omp_get_thread_num();
        for(int bi = 0; bi < m->spmv_block_per_thread; bi++) {
            const int bb = m->spmv_block_separator[bi * m->spmv_max_threads + thread_num];
            const int be = m->spmv_block_separator[bi * m->spmv_max_threads + thread_num + 1];
            if(bb != be) {
                double *restrict local_Y[MSPMV_NUM];
                for(int i = 0; i < MSPMV_NUM; i++) {
                    local_Y[i] = Y[i] + bb;
                }

                #pragma loop count min(8)
                for(int i = bb; i < be; i++) {
                    for(int j = 0; j < MSPMV_NUM; j++) {
                        double y = 0.0;
                        const int b = m->index[i];
                        const int e = m->index[i + 1];
                        for(int k = b; k < e; k++) {
                            const double a = m->A[k];
                            const int p = m->item[k];
                            const double x = X[j][p];
                            y += a * x;
                        }
                        local_Y[j][i - bb] = y;
                    }
                }
            }
        }
    }
}

void main(int argn, char** args) {
    assert(argn == 4);
    const char* matrix_file = args[1];
    const char* order_file = args[2];
    const int mspmv_repeat = atoi(args[3]);

    raw_matrix rm;
    read_matrix(matrix_file, &rm);

    int *order = (int*)malloc(sizeof(int) * rm.n);
    FILE *f = fopen(order_file, "r");
    for(int i = 0; i < rm.n; i++) {
        int _ = fscanf(f, "%d", order + i);
    }
    fclose(f);

    double *restrict x[MSPMV_NUM];
    double *restrict y[MSPMV_NUM];
    for(int i = 0; i < MSPMV_NUM; i++) {
        x[i] = (double*)malloc(sizeof(double) * rm.n);
        y[i] = (double*)malloc(sizeof(double) * rm.n);
    }

    reordered_matrix m;
    reorder_matrix(&rm, order, &m);

    double c = 0.0;
    for(int i = 0; i < mspmv_repeat; i++) {
        #pragma omp parallel default(none), shared(m, x)
        {
            const int thread_num = omp_get_thread_num();
            for(int bi = 0; bi < m.spmv_block_per_thread; bi++) {
                const int bb = m.spmv_block_separator[bi * m.spmv_max_threads + thread_num];
                const int be = m.spmv_block_separator[bi * m.spmv_max_threads + thread_num + 1];
                for(int j = 0; j < MSPMV_NUM; j++) {
                    for(int k = bb; k < be; k++) {
                        x[j][k] = 1.1;
                    }
                }
            }
        }
        double t = omp_get_wtime();
        MSpMV(&m, x, y);
        c += omp_get_wtime() - t;
    }
    printf("%le %le %le\n", c, c / mspmv_repeat, c / mspmv_repeat / MSPMV_NUM);
}
