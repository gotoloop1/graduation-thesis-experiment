#include <stdlib.h>
#include "lib.h"
#include "reordering.h"

int rcm_select_first_vertex(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge
) {
    int v = 0;
    for(int i = 1; i < n; i++) {
        if(edge_separator[i + 1] - edge_separator[i] < edge_separator[v + 1] - edge_separator[v]) {
            v = i;
        }
    }
    return v;
}

void rcm_reverse_order(
    int n,
    int *restrict order
) {
    for(int i = 0; i < n / 2; i++) {
        const int t = order[i];
        order[i] = order[n - 1 - i];
        order[n - 1 - i] = t;
    }
}

void cuthill_mckee_inner(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order,

    int *restrict _2n_buffer, /* longer than 2*n */
    bool *restrict b_n_buffer /* longer than n */
) {
    const int first_v = rcm_select_first_vertex(n, edge_separator, edge);

    bool* used = b_n_buffer;
    int* q = _2n_buffer;

    int order_size = 1;
    int front_pos = 0;
    order[0] = first_v;
    for(int i = 0; i < n; i++) {
        used[i] = false;
    }
    used[first_v] = true;
    while(order_size < n) {
        const int p = order[front_pos];
        front_pos++;
        int q_size = 0;
        for(int j = edge_separator[p]; j < edge_separator[p + 1]; j++) {
            if(!used[edge[j]]) {
                used[edge[j]] = true;
                q[q_size * 2] = edge_separator[edge[j] + 1] - edge_separator[edge[j]];
                q[q_size * 2 + 1] = edge[j];
                q_size++;
            }
        }
        first_element_sort(q_size, q, n - 1);
        for(int i = 0; i < q_size; i++) {
            order[order_size] = q[i * 2 + 1];
            order_size++;
        }
    }
}

void cuthill_mckee(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
) {
    int *sub_vertex = (int*)malloc(sizeof(int) * n);
    int *sub_vertex_inverse = (int*)malloc(sizeof(int) * n);
    int *sub_edge_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *sub_edge = (int*)malloc(sizeof(int) * edge_separator[n]);
    int *stack = (int*)malloc(sizeof(int) * n);

    int *_2n_buffer = (int*)malloc(sizeof(int) * 2 * n);
    bool *b_n_buffer = (bool*)malloc(sizeof(bool) * n);

    for(int i = 0; i < n; i++) {
        sub_vertex_inverse[i] = -1;
    }
    int count = 0;
    for(int i = 0; i < n; i++) {
        if(sub_vertex_inverse[i] < 0) {
            int sub_n = 1;
            sub_vertex[0] = i;
            sub_vertex_inverse[i] = 1;
            int stack_size = 1;
            stack[0] = i;
            while(stack_size > 0) {
                stack_size--;
                const int p = stack[stack_size];
                for(int j = edge_separator[p]; j < edge_separator[p + 1]; j++) {
                    const int t = edge[j];
                    if(sub_vertex_inverse[t] < 0) {
                        sub_vertex_inverse[t] = 1;
                        sub_vertex[sub_n] = t;
                        sub_n++;
                        stack[stack_size] = t;
                        stack_size++;
                    }
                }
            }

            single_sort(sub_n, sub_vertex, n - 1);
            for(int j = 0; j < sub_n; j++) {
                sub_vertex_inverse[sub_vertex[j]] = j;
            }

            sub_edge_separator[0] = 0;
            for(int j = 0; j < sub_n; j++) {
                sub_edge_separator[j + 1] = sub_edge_separator[j];
                const int p = sub_vertex[j];
                for(int k = edge_separator[p]; k < edge_separator[p + 1]; k++) {
                    sub_edge[sub_edge_separator[j + 1]] = sub_vertex_inverse[edge[k]];
                    sub_edge_separator[j + 1]++;
                }
            }

            cuthill_mckee_inner(sub_n, sub_edge_separator, sub_edge, order + count, _2n_buffer, b_n_buffer);

            for(int j = count; j < count + sub_n; j++) {
                order[j] = sub_vertex[order[j]];
            }
            count += sub_n;
        }
    }

    free(sub_vertex);
    free(sub_vertex_inverse);
    free(sub_edge_separator);
    free(sub_edge);
    free(stack);
    free(_2n_buffer);
    free(b_n_buffer);
}

void reverse_cuthill_mckee(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
) {
    cuthill_mckee(n, edge_separator, edge, order);
    rcm_reverse_order(n, order);
}
