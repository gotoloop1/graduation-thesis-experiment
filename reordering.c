#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "lib.h"
#include "matrix.h"
#include "reordering.h"

typedef struct {
    int n;
    int *restrict edge_separator;
    int *restrict edge;
    int *restrict edge_weight;
} graph;

void generate_graph(
    const raw_matrix *restrict m,
    graph *restrict g
) {
    int *count = (int*)malloc(sizeof(int) * m->n);
    int *buffer_separator = (int*)malloc(sizeof(int) * (m->n + 1));
    int *buffer = (int*)malloc(sizeof(int) * 2 * m->nz);

    for(int i = 0; i < m->n; i++) {
        count[i] = 0;
    }
    for(int i = 0; i < m->nz; i++) {
        if(m->row[i] != m->column[i]) {
            count[m->row[i]]++;
            count[m->column[i]]++;
        }
    }

    buffer_separator[0] = 0;
    for(int i = 0; i < m->n; i++) {
        buffer_separator[i + 1] = buffer_separator[i] + count[i];
        count[i] = 0;
    }

    for(int i = 0; i < m->nz; i++) {
        if(m->row[i] != m->column[i]) {
            buffer[buffer_separator[m->row[i]] + count[m->row[i]]] = m->column[i];
            count[m->row[i]]++;
            buffer[buffer_separator[m->column[i]] + count[m->column[i]]] = m->row[i];
            count[m->column[i]]++;
        }
    }
    int all_count = 0;
    for(int i = 0; i < m->n; i++) {
        const int b = buffer_separator[i];
        const int e = buffer_separator[i + 1];
        single_sort(e - b, buffer + b, m->n - 1);
        int pre = -1;
        for(int j = b; j < e; j++) {
            if(buffer[j] != pre) {
                all_count++;
            }
            pre = buffer[j];
        }
    }

    g->n = m->n;
    g->edge_separator = (int*)malloc(sizeof(int) * (m->n + 1));
    g->edge = (int*)malloc(sizeof(int) * all_count);
    g->edge_weight = (int*)malloc(sizeof(int) * all_count);
    g->edge_separator[0] = 0;
    for(int i = 0; i < m->n; i++) {
        g->edge_separator[i + 1] = g->edge_separator[i];
        const int b = buffer_separator[i];
        const int e = buffer_separator[i + 1];
        int pre = -1;
        for(int j = b; j < e; j++) {
            if(buffer[j] == pre) {
                g->edge_weight[g->edge_separator[i + 1] - 1]++;
            } else {
                g->edge[g->edge_separator[i + 1]] = buffer[j];
                g->edge_weight[g->edge_separator[i + 1]] = 1;
                g->edge_separator[i + 1]++;
            }
            pre = buffer[j];
        }
    }

    free(count);
    free(buffer_separator);
    free(buffer);
}

void main(int argn, char** args) {
    assert(argn == 5);
    const char *matrix_file = args[1];
    const char *order_file = args[2];
    const int ordering_type = atoi(args[3]);
    const unsigned int seed = atoi(args[4]);

    xor128_set_seed(seed);

    raw_matrix m;
    read_matrix(matrix_file, &m);

    int *pre_random = (int*)malloc(sizeof(int) * m.n);
    generate_random_permutation(m.n, pre_random);
    int *pre_random_inverse = (int*)malloc(sizeof(int) * m.n);
    for(int i = 0; i < m.n; i++) {
        pre_random_inverse[pre_random[i]] = i;
    }
    for(int i = 0; i < m.nz; i++) {
        m.row[i] = pre_random[m.row[i]];
        m.column[i] = pre_random[m.column[i]];
    }

    graph g;
    generate_graph(&m, &g);

    int *order = (int*)malloc(sizeof(int) * m.n);

    clock_t c = clock();
    switch(ordering_type) {
    case 0:
        for(int i = 0; i < m.n; i++) {
            order[i] = i;
        }
        break;
    case 1:
        for(int i = 0; i < m.n; i++) {
            order[i] = pre_random[i];
        }
        break;
    case 2:
        reverse_cuthill_mckee(g.n, g.edge_separator, g.edge, order);
        break;
    case 3:
        edge_based_cost_minimization(g.n, g.edge_separator, g.edge, order);
        break;
    case 4:
        node_based_cost_minimization(g.n, g.edge_separator, g.edge, order);
        break;
    case 5:
        hybrid_cost_minimization(g.n, g.edge_separator, g.edge, order);
        break;
    case 6:
        coarsen_refine_hybrid_cost_minimization(g.n, g.edge_separator, g.edge, order);
        break;
    }
    printf("%le\n", (double)(clock() - c) / CLOCKS_PER_SEC);

    for(int i = 0; i < m.n; i++) {
        order[i] = pre_random_inverse[order[i]];
    }

    FILE *f = fopen(order_file, "w");
    for(int i = 0; i < m.n; i++) {
        fprintf(f, "%d\n", order[i]);
    }
    fclose(f);
}
