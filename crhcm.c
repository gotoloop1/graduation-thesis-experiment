#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"
#include "reordering.h"

#define MAX_SHRINK_DEPTH 10000

double crhcm_calculate_separation_cost(
    const int n,
    const int adjacent_n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict from_adjacent_separator,
    const int *restrict from_adjacent,
    const double *restrict node_weight_first,
    const double *restrict node_weight_second,
    const double *restrict adjacent_node_weight_first,
    const double *restrict adjacent_node_weight_second,
    const double *restrict edge_cost_first,
    const double edge_cost_border,
    const double *restrict edge_cost_second,
    const int *restrict belongs
) {
    double sum = 0.0;
    for(int i = 0; i < n; i++) {
        int first_count = 0;
        int second_count = 0;
        if(belongs[i] == 0) {
            first_count++;
            sum += edge_cost_first[i];
        } else {
            second_count++;
            sum += edge_cost_second[i];
        }
        for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
            const int p = edge[j];
            if(belongs[p] == 0) {
                first_count++;
            } else {
                second_count++;
            }
            
            if(i < p && belongs[i] != belongs[p]) {
                sum += edge_cost_border;
            }
        }
        double f, b, s;
        if(node_weight_first[i] > node_weight_second[i]) {
            f = 0.0;
            b = node_weight_first[i];
            s = node_weight_first[i] - node_weight_second[i];
        } else {
            f = node_weight_second[i] - node_weight_first[i];
            b = node_weight_second[i];
            s = 0.0;
        }
        if(first_count > 0 && second_count > 0) {
            sum += b;
        } else if(first_count > 0) {
            sum += f;
        } else { // second_count > 0
            sum += s;
        }
    }
    for(int i = 0; i < adjacent_n; i++) {
        int first_count = 0;
        int second_count = 0;
        for(int j = from_adjacent_separator[i]; j < from_adjacent_separator[i + 1]; j++) {
            if(belongs[from_adjacent[j]] == 0) {
                first_count++;
            } else {
                second_count++;
            }
        }
        double f, b, s;
        if(adjacent_node_weight_first[i] > adjacent_node_weight_second[i]) {
            f = 0.0;
            b = adjacent_node_weight_first[i];
            s = adjacent_node_weight_first[i] - adjacent_node_weight_second[i];
        } else {
            f = adjacent_node_weight_second[i] - adjacent_node_weight_first[i];
            b = adjacent_node_weight_second[i];
            s = 0.0;
        }
        if(first_count > 0 && second_count > 0) {
            sum += b;
        } else if(first_count > 0) {
            sum += f;
        } else { // second_count > 0
            sum += s;
        }
    }
    return sum;
}

void crhcm_first_separation(
    const bool is_first,
    const int n,
    const double *restrict edge_cost_first,
    const double *restrict edge_cost_second,
    int *restrict belongs,

    int *restrict _2n_buffer
) {
    const int MAX_WEIGHT = 10000;
    const double DISTURB_RATE = 0.25;

    double min_val = edge_cost_first[0] - edge_cost_second[0];
    double max_val = min_val;
    for(int i = 1; i < n; i++) {
        const double d = edge_cost_first[i] - edge_cost_second[i];
        if(min_val > d) {
            min_val = d;
        }
        if(max_val < d) {
            max_val = d;
        }
    }
    int *tmp = _2n_buffer;
    if(max_val - min_val < 0.001) {
        for(int i = 0; i < n; i++) {
            tmp[i * 2] = xor128() % (MAX_WEIGHT + 1);
            tmp[i * 2 + 1] = i;
        }
    } else {
        for(int i = 0; i < n; i++) {
            int d = (edge_cost_first[i] - edge_cost_second[i] - min_val) * MAX_WEIGHT / (max_val - min_val);
            if(!is_first) {
                if(xor128() % 10000 < DISTURB_RATE * 10000) {
                    if(xor128() % 2 == 0) {
                        d = 0;
                    } else {
                        d = MAX_WEIGHT;
                    }
                }
            }
            tmp[i * 2] = d;
            tmp[i * 2 + 1] = i;
        }
    }

    first_element_sort(n, tmp, MAX_WEIGHT);

    for(int i = 0; i < n / 2; i++) {
        belongs[tmp[i * 2 + 1]] = 0;
    }
    for(int i = n / 2; i < n; i++) {
        belongs[tmp[i * 2 + 1]] = 1;
    }
}

#define REGISTER_CHANGE_FIRST(pos, value) if(changed[(pos)] < 0) { changed[(pos)] = changed_first_size; changed_first_pos[changed_first_size] = (pos); changed_first_value[changed_first_size] = (value); changed_first_size++; } else { changed_first_value[changed[(pos)]] += (value); }
#define REGISTER_CHANGE_SECOND(pos, value) if(changed[(pos)] < 0) { changed[(pos)] = changed_second_size; changed_second_pos[changed_second_size] = (pos); changed_second_value[changed_second_size] = (value); changed_second_size++; } else { changed_second_value[changed[(pos)]] += (value); }

double crhcm_fiduccia_mattheyses_inner(
    const int n,
    const int adjacent_n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict to_adjacent_separator,
    const int *restrict to_adjacent,
    const int *restrict from_adjacent_separator,
    const int *restrict from_adjacent,
    const double *restrict node_weight_first,
    const double *restrict node_weight_second,
    const double *restrict adjacent_node_weight_first,
    const double *restrict adjacent_node_weight_second,
    const double *restrict edge_cost_first,
    const double edge_cost_border,
    const double *restrict edge_cost_second,
    const double error_tolerance,
    int *restrict belongs,

    heap *restrict h1,
    heap *restrict h2,
    double *restrict d_n_buffer_1,
    double *restrict d_n_buffer_2,
    double *restrict d_n_buffer_3,
    double *restrict d_n_buffer_4,
    double *restrict d_n_buffer_5,
    int *restrict n_buffer_1,
    int *restrict n_buffer_2,
    int *restrict n_buffer_3,
    int *restrict n_buffer_4,
    int *restrict n_buffer_5,
    int *restrict n_buffer_6,
    int *restrict n_buffer_7,
    int *restrict n_buffer_8,
    int *restrict n_buffer_9,
    int *restrict n_buffer_10,
    int *restrict n_buffer_11,
    bool *restrict b_n_buffer
) {
    const double TORELANCE = 1.3;

    int *first_count = n_buffer_1;
    int *second_count = n_buffer_2;
    for(int i = 0; i < n; i++) {
        first_count[i] = 0;
        second_count[i] = 0;
        for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
            if(belongs[edge[j]] == 0) {
                first_count[i]++;
            } else {
                second_count[i]++;
            }
        }
    }
    int *adjacent_first_count = n_buffer_3;
    int *adjacent_second_count = n_buffer_4;
    for(int i = 0; i < adjacent_n; i++) {
        adjacent_first_count[i] = 0;
        adjacent_second_count[i] = 0;
        for(int j = from_adjacent_separator[i]; j < from_adjacent_separator[i + 1]; j++) {
            if(belongs[from_adjacent[j]] == 0) {
                adjacent_first_count[i]++;
            } else {
                adjacent_second_count[i]++;
            }
        }
    }

    heap *h_first = h1;
    heap *h_second = h2;
    int *bucket_first = n_buffer_5;
    int *bucket_second = n_buffer_6;
    double *bucket_gain_first = d_n_buffer_1;
    double *bucket_gain_second = d_n_buffer_2;
    int *inverse = n_buffer_7;
    int first_size = 0;
    int second_size = 0;
    for(int i = 0; i < n; i++) {
        if(belongs[i] == 0) {
            inverse[i] = first_size;
            bucket_first[first_size] = i;
            bucket_gain_first[first_size] = edge_cost_first[i] - edge_cost_second[i];
            if(first_count[i] + second_count[i] == 0) { // F -> S @ i
                bucket_gain_first[first_size] += node_weight_second[i] - node_weight_first[i];
            } else {
                if(first_count[i] == 0) { // B -> S @ i
                    bucket_gain_first[first_size] += node_weight_second[i];
                } else if(second_count[i] == 0) { // F -> B @ i
                    bucket_gain_first[first_size] -= node_weight_first[i];
                }
            }
            for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
                const int p = edge[j];
                if(belongs[p] == 0) {
                    bucket_gain_first[first_size] -= edge_cost_border;
                    if(second_count[p] == 0) { // F -> B @ p
                        bucket_gain_first[first_size] -= node_weight_first[p];
                    }
                } else {
                    bucket_gain_first[first_size] += edge_cost_border;
                    if(first_count[p] == 1) { // B -> S @ p
                        bucket_gain_first[first_size] += node_weight_second[p];
                    }
                }
            }
            for(int j = to_adjacent_separator[i]; j < to_adjacent_separator[i + 1]; j++) {
                const int p = to_adjacent[j];
                if(adjacent_first_count[p] == 1 && adjacent_second_count[p] == 0) { // F -> S @ p
                    bucket_gain_first[first_size] += adjacent_node_weight_second[p] - adjacent_node_weight_first[p];
                } else if(adjacent_first_count[p] == 1) { // B -> S @ p
                    bucket_gain_first[first_size] += adjacent_node_weight_second[p];
                } else if(adjacent_second_count[p] == 0) { // F -> B @ p
                    bucket_gain_first[first_size] -= adjacent_node_weight_first[p];
                }
            }
            first_size++;
        } else {
            inverse[i] = second_size;
            bucket_second[second_size] = i;
            bucket_gain_second[second_size] = edge_cost_second[i] - edge_cost_first[i];
            if (first_count[i] + second_count[i] == 0) { // S -> F @ i
                bucket_gain_second[second_size] += node_weight_first[i] - node_weight_second[i];
            } else {
                if(first_count[i] == 0) { // S -> B @ i
                    bucket_gain_second[second_size] -= node_weight_second[i];
                } else if(second_count[i] == 0) { // B -> F @ i
                    bucket_gain_second[second_size] += node_weight_first[i];
                }
            }
            for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
                const int p = edge[j];
                if(belongs[p] == 0) {
                    bucket_gain_second[second_size] += edge_cost_border;
                    if(second_count[p] == 1) { // B -> F @ p
                        bucket_gain_second[second_size] += node_weight_first[p];
                    }
                } else {
                    bucket_gain_second[second_size] -= edge_cost_border;
                    if(first_count[p] == 0) { // S -> B @ p
                        bucket_gain_second[second_size] -= node_weight_second[p];
                    }
                }
            }
            for(int j = to_adjacent_separator[i]; j < to_adjacent_separator[i + 1]; j++) {
                const int p = to_adjacent[j];
                if(adjacent_first_count[p] == 0 && adjacent_second_count[p] == 1) { // S -> F @ p
                    bucket_gain_second[second_size] += adjacent_node_weight_first[p] - adjacent_node_weight_second[p];
                } else if(adjacent_first_count[p] == 0) { // S -> B @ p
                    bucket_gain_second[second_size] -= adjacent_node_weight_second[p];
                } else if(adjacent_second_count[p] == 1) { // B -> F @ p
                    bucket_gain_second[second_size] += adjacent_node_weight_first[p];
                }
            }
            second_size++;
        }
    }
    heap_initialize(first_size, bucket_gain_first, h_first);
    heap_initialize(second_size, bucket_gain_second, h_second);

    bool *used = b_n_buffer;
    for(int i = 0; i < n; i++) {
        used[i] = false;
    }
    int moved_size = 0;
    int *moved = n_buffer_8;
    double *moved_gain = d_n_buffer_3;
    int current_first_size = first_size;
    int current_second_size = second_size;

    int *changed = n_buffer_9;
    for(int i = 0; i < n; i++) {
        changed[i] = -1;
    }
    int *changed_first_pos = n_buffer_10;
    double *changed_first_value = d_n_buffer_4;
    int *changed_second_pos = n_buffer_11;
    double *changed_second_value = d_n_buffer_5;
    while(1) {
        const bool c_first = !heap_empty(h_first);
        const bool c_second = !heap_empty(h_second);
        int m;
        if(c_first && c_second) {
            if(current_first_size > current_second_size && (current_first_size + 1) > (current_second_size - 1) * TORELANCE) {
                m = 0;
            } else if(current_second_size > current_first_size && (current_second_size + 1) > (current_first_size - 1) * TORELANCE) {
                m = 1;
            } else {
                const double d_first = heap_front(h_first).data;
                const double d_second = heap_front(h_second).data;
                if(d_first > d_second) {
                    m = 0;
                } else {
                    m = 1;
                }
            }
        } else if(c_first) {
            m = 0;
        } else if(c_second) {
            m = 1;
        } else {
            break;
        }

        int changed_first_size = 0;
        int changed_second_size = 0;
        if(m == 0) {
            const int p = bucket_first[heap_front(h_first).pos];
            moved[moved_size] = p;
            moved_gain[moved_size] = heap_front(h_first).data;
            heap_pop(h_first);
            used[p] = true;

            for(int i = edge_separator[p]; i < edge_separator[p + 1]; i++) {
                const int q = edge[i];
                if(belongs[q] == 0) {
                    if(!used[q]) {
                        double sum = 2.0 * edge_cost_border;
                        if(first_count[p] == 1) {
                            sum += node_weight_second[p];
                        }
                        if(second_count[p] == 0) {
                            sum += node_weight_first[p];
                        }
                        if(first_count[q] == 1) {
                            sum += node_weight_second[q];
                        }
                        if(second_count[q] == 0) {
                            sum += node_weight_first[q];
                        }
                        if(fabs(sum) > 0.000001) {
                            REGISTER_CHANGE_FIRST(q, sum);
                        }
                    }

                    if(second_count[q] == 0) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t]) { // (F -> B) -> (B -> B) @ q -> t(0)
                                REGISTER_CHANGE_FIRST(t, node_weight_first[q]);
                            }
                        }
                    } else if(second_count[q] == 1) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t] && belongs[t] == 1) { // (B -> F) -> (B -> B) @ q -> t(1)
                                REGISTER_CHANGE_SECOND(t, -node_weight_first[q]);
                            }
                        }
                    }
                } else {
                    if(!used[q]) {
                        double sum = -2.0 * edge_cost_border;
                        if(first_count[p] == 0) {
                            sum -= node_weight_second[p];
                        }
                        if(second_count[p] == 1) {
                            sum -= node_weight_first[p];
                        }
                        if(first_count[q] == 1) {
                            sum -= node_weight_second[q];
                        }
                        if(second_count[q] == 0) {
                            sum -= node_weight_first[q];
                        }
                        if(fabs(sum) > 0.000001) {
                            REGISTER_CHANGE_SECOND(q, sum);
                        }
                    }

                    if(first_count[q] == 1) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t]) { // (B -> B) -> (S -> B) @ q -> t(1)
                                REGISTER_CHANGE_SECOND(t, -node_weight_second[q]);
                            }
                        }
                    } else if(first_count[q] == 2) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t] && belongs[t] == 0) { // (B -> B) -> (B -> S) @ q -> t(0)
                                REGISTER_CHANGE_FIRST(t, node_weight_second[q]);
                            }
                        }
                    }
                }
            }

            for(int i = to_adjacent_separator[p]; i < to_adjacent_separator[p + 1]; i++) {
                const int q = to_adjacent[i];
                if(adjacent_first_count[q] == 2 && adjacent_second_count[q] == 0) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (F -> B) -> (B -> S) @ q -> t(0)
                            REGISTER_CHANGE_FIRST(t, adjacent_node_weight_first[q] + adjacent_node_weight_second[q]);
                        }
                    }
                } else if(adjacent_first_count[q] == 1 && adjacent_second_count[q] == 1) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (B -> F) -> (S -> B) @ q -> t(1)
                            REGISTER_CHANGE_SECOND(t, -adjacent_node_weight_first[q] - adjacent_node_weight_second[q]);
                        }
                    }
                } else if(adjacent_first_count[q] == 1) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (B -> B) -> (S -> B) @ q -> t(1)
                            REGISTER_CHANGE_SECOND(t, -adjacent_node_weight_second[q]);
                        }
                    }
                } else if(adjacent_second_count[q] == 0) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (F -> B) -> (B -> B) @ q -> t(0)
                            REGISTER_CHANGE_FIRST(t, adjacent_node_weight_first[q]);
                        }
                    }
                } else if(adjacent_first_count[q] == 2 || adjacent_second_count[q] == 1) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) {
                            if(adjacent_first_count[q] == 2 && belongs[t] == 0) { // (B -> B) -> (B -> S) @ q -> t(0)
                                REGISTER_CHANGE_FIRST(t, adjacent_node_weight_second[q]);
                            } else if(adjacent_second_count[q] == 1 && belongs[t] == 1) { // (B -> F) -> (B -> B) @ q -> t(1)
                                REGISTER_CHANGE_SECOND(t, -adjacent_node_weight_first[q]);
                            }
                        }
                    }
                }
            }

            for(int i = edge_separator[p]; i < edge_separator[p + 1]; i++) {
                const int q = edge[i];
                first_count[q]--;
                second_count[q]++;
            }

            for(int i = to_adjacent_separator[p]; i < to_adjacent_separator[p + 1]; i++) {
                const int q = to_adjacent[i];
                adjacent_first_count[q]--;
                adjacent_second_count[q]++;
            }

            belongs[p] = 1;
            current_first_size--;
            current_second_size++;
        } else {
            const int p = bucket_second[heap_front(h_second).pos];
            moved[moved_size] = p;
            moved_gain[moved_size] = heap_front(h_second).data;
            heap_pop(h_second);
            used[p] = true;

            for(int i = edge_separator[p]; i < edge_separator[p + 1]; i++) {
                const int q = edge[i];
                if(belongs[q] == 0) {
                    if(!used[q]) {
                        double sum = -2.0 * edge_cost_border;
                        if(first_count[p] == 1) {
                            sum -= node_weight_second[p];
                        }
                        if(second_count[p] == 0) {
                            sum -= node_weight_first[p];
                        }
                        if(first_count[q] == 0) {
                            sum -= node_weight_second[q];
                        }
                        if(second_count[q] == 1) {
                            sum -= node_weight_first[q];
                        }
                        if(fabs(sum) > 0.000001) {
                            REGISTER_CHANGE_FIRST(q, sum);
                        }
                    }

                    if(second_count[q] == 1) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t]) { // (B -> B) -> (F -> B) @ q -> t(0)
                                REGISTER_CHANGE_FIRST(t, -node_weight_first[q]);
                            }
                        }
                    } else if(second_count[q] == 2) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t] && belongs[t] == 1) { // (B -> B) -> (B -> F) @ q -> t(1)
                                REGISTER_CHANGE_SECOND(t, node_weight_first[q]);
                            }
                        }
                    }
                } else {
                    if(!used[q]) {
                        double sum = 2.0 * edge_cost_border;
                        if(first_count[p] == 0) {
                            sum += node_weight_second[p];
                        }
                        if(second_count[p] == 1) {
                            sum += node_weight_first[p];
                        }
                        if(first_count[q] == 0) {
                            sum += node_weight_second[q];
                        }
                        if(second_count[q] == 1) {
                            sum += node_weight_first[q];
                        }
                        if(fabs(sum) > 0.000001) {
                            REGISTER_CHANGE_SECOND(q, sum);
                        }
                    }

                    if(first_count[q] == 0) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t]) { // (S -> B) -> (B -> B) @ q -> t(1)
                                REGISTER_CHANGE_SECOND(t, node_weight_second[q]);
                            }
                        }
                    } else if(first_count[q] == 1) {
                        for(int j = edge_separator[q]; j < edge_separator[q + 1]; j++) {
                            const int t = edge[j];
                            if(!used[t] && belongs[t] == 0) { // (B -> S) -> (B -> B) @ q -> t(0)
                                REGISTER_CHANGE_FIRST(t, -node_weight_second[q]);
                            }
                        }
                    }
                }
            }

            for(int i = to_adjacent_separator[p]; i < to_adjacent_separator[p + 1]; i++) {
                const int q = to_adjacent[i];
                if(adjacent_first_count[q] == 1 && adjacent_second_count[q] == 1) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (B -> S) -> (F -> B) @ q -> t(0)
                            REGISTER_CHANGE_FIRST(t, -adjacent_node_weight_first[q] - adjacent_node_weight_second[q]);
                        }
                    }
                } else if(adjacent_first_count[q] == 0 && adjacent_second_count[q] == 2) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (S -> B) -> (B -> F) @ q -> t(1)
                            REGISTER_CHANGE_SECOND(t, adjacent_node_weight_first[q] + adjacent_node_weight_second[q]);
                        }
                    }
                } else if(adjacent_first_count[q] == 0) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (S -> B) -> (B -> B) @ q -> t(1)
                            REGISTER_CHANGE_SECOND(t, adjacent_node_weight_second[q]);
                        }
                    }
                } else if(adjacent_second_count[q] == 1) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) { // (B -> B) -> (F -> B) @ q -> t(0)
                            REGISTER_CHANGE_FIRST(t, -adjacent_node_weight_first[q]);
                        }
                    }
                } else if(adjacent_first_count[q] == 1 || adjacent_second_count[q] == 2) {
                    for(int j = from_adjacent_separator[q]; j < from_adjacent_separator[q + 1]; j++) {
                        const int t = from_adjacent[j];
                        if(!used[t]) {
                            if(adjacent_first_count[q] == 1 && belongs[t] == 0) { // (B -> S) -> (B -> B) @ q -> t(0)
                                REGISTER_CHANGE_FIRST(t, -adjacent_node_weight_second[q]);
                            } else if(adjacent_second_count[q] == 2 && belongs[t] == 1) { // (B -> B) -> (B -> F) @ q -> t(1)
                                REGISTER_CHANGE_SECOND(t, adjacent_node_weight_first[q]);
                            }
                        }
                    }
                }
            }

            for(int i = edge_separator[p]; i < edge_separator[p + 1]; i++) {
                const int q = edge[i];
                first_count[q]++;
                second_count[q]--;
            }

            for(int i = to_adjacent_separator[p]; i < to_adjacent_separator[p + 1]; i++) {
                const int q = to_adjacent[i];
                adjacent_first_count[q]++;
                adjacent_second_count[q]--;
            }

            belongs[p] = 0;
            current_first_size++;
            current_second_size--;
        }

        for(int i = 0; i < changed_first_size; i++) {
            const int p = changed_first_pos[i];
            const double v = changed_first_value[i];
            heap_update(inverse[p], v, h_first);
            changed[p] = -1;
        }
        for(int i = 0; i < changed_second_size; i++) {
            const int p = changed_second_pos[i];
            const double v = changed_second_value[i];
            heap_update(inverse[p], v, h_second);
            changed[p] = -1;
        }

        moved_size++;
    }

    double max_val = 0.0;
    int max_val_pos = -1;
    double current_gain = 0.0;
    current_first_size = first_size;
    current_second_size = second_size;
    for(int i = 0; i < moved_size; i++) {
        if(belongs[moved[i]] == 1) {
            current_first_size--;
            current_second_size++;
        } else {
            current_first_size++;
            current_second_size--;
        }
        current_gain += moved_gain[i];
        if((current_first_size - current_second_size > 1 && current_first_size > current_second_size * TORELANCE) || (current_second_size - current_first_size > 1 && current_second_size > current_first_size * TORELANCE)) {
            continue;
        }
        if(max_val + error_tolerance < current_gain) {
            max_val = current_gain;
            max_val_pos = i;
        }
    }
    if(max_val_pos < 0) {
        for(int i = 0; i < n; i++) {
            belongs[i] ^= 1;
        }
        return -2.0;
    }
    for(int i = max_val_pos + 1; i < n; i++) {
        const int p = moved[i];
        belongs[p] ^= 1;
    }
    return max_val;
}

void crhcm_fiduccia_mattheyses(
    const int retry_number,
    const int n,
    const int adjacent_n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict to_adjacent_separator,
    const int *restrict to_adjacent,
    const int *restrict from_adjacent_separator,
    const int *restrict from_adjacent,
    const double *restrict node_weight_first,
    const double *restrict node_weight_second,
    const double *restrict adjacent_node_weight_first,
    const double *restrict adjacent_node_weight_second,
    const double *restrict edge_cost_first,
    const double edge_cost_border,
    const double *restrict edge_cost_second,
    int *restrict belongs,

    heap *restrict h1,
    heap *restrict h2,
    double *restrict d_n_buffer_1,
    double *restrict d_n_buffer_2,
    double *restrict d_n_buffer_3,
    double *restrict d_n_buffer_4,
    double *restrict d_n_buffer_5,
    int *restrict _2n_buffer,
    int *restrict n_buffer_1,
    int *restrict n_buffer_2,
    int *restrict n_buffer_3,
    int *restrict n_buffer_4,
    int *restrict n_buffer_5,
    int *restrict n_buffer_6,
    int *restrict n_buffer_7,
    int *restrict n_buffer_8,
    int *restrict n_buffer_9,
    int *restrict n_buffer_10,
    int *restrict n_buffer_11,
    bool *restrict b_n_buffer
) {
    double min_cost;
    for(int i = 0; i < retry_number; i++) {
        int *current_belongs = n_buffer_11;
        crhcm_first_separation(
            i == 0,
            n,
            edge_cost_first,
            edge_cost_second,
            current_belongs,
            _2n_buffer
        );
        double cost = crhcm_calculate_separation_cost(
            n,
            adjacent_n,
            edge_separator,
            edge,
            from_adjacent_separator,
            from_adjacent,
            node_weight_first,
            node_weight_second,
            adjacent_node_weight_first,
            adjacent_node_weight_second,
            edge_cost_first,
            edge_cost_border,
            edge_cost_second,
            current_belongs
        );
        while(1) {
            double error_tolerance = cost / pow(10.0, 9.0);
            if(error_tolerance < 0.000001) {
                error_tolerance = 0.000001;
            }
            double improve = crhcm_fiduccia_mattheyses_inner(
                n,
                adjacent_n,
                edge_separator,
                edge,
                to_adjacent_separator,
                to_adjacent,
                from_adjacent_separator,
                from_adjacent,
                node_weight_first,
                node_weight_second,
                adjacent_node_weight_first,
                adjacent_node_weight_second,
                edge_cost_first,
                edge_cost_border,
                edge_cost_second,
                error_tolerance,
                current_belongs,
                h1,
                h2,
                d_n_buffer_1,
                d_n_buffer_2,
                d_n_buffer_3,
                d_n_buffer_4,
                d_n_buffer_5,
                _2n_buffer,
                n_buffer_1,
                n_buffer_2,
                n_buffer_3,
                n_buffer_4,
                n_buffer_5,
                n_buffer_6,
                n_buffer_7,
                n_buffer_8,
                n_buffer_9,
                n_buffer_10,
                b_n_buffer
            );
            if(improve < -1.0) {
                break;
            }
            cost -= improve;
        }
        cost = crhcm_calculate_separation_cost(
            n,
            adjacent_n,
            edge_separator,
            edge,
            from_adjacent_separator,
            from_adjacent,
            node_weight_first,
            node_weight_second,
            adjacent_node_weight_first,
            adjacent_node_weight_second,
            edge_cost_first,
            edge_cost_border,
            edge_cost_second,
            current_belongs
        );
        if(i == 0 || min_cost > cost) {
            min_cost = cost;
            for(int j = 0; j < n; j++) {
                belongs[j] = current_belongs[j];
            }
        }
    }
}

double crhcm_calculate_cost(
    const weighted_range_sum_data *wrs,
    int x_order_begin,
    int x_order_end,
    int y_order_begin,
    int y_order_end
) {
    return box_sum(wrs, x_order_begin, x_order_end - 1, y_order_begin, y_order_end - 1) / (x_order_end - x_order_begin) / (y_order_end - y_order_begin);
}

double crhcm_calculate_single_cost(
    const int cost_type,
    const int distance
) {
    if(distance == 0) {
        return 0.0;
    } else {
        switch(cost_type) {
        case 0:
            return distance;
        case 1:
            return sqrt(distance);
        case 2:
            return log2(distance);
        case 3:
            return sqrt(log2(distance));
        case 4:
            return pow(log2(distance), 0.25);
        }
    }
}

int crhcm_shrink(
    const int min_shrink_size,
    const int shrink_depth,
    const int n,
    int *restrict node_weight,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict edge_weight,

    int *restrict layerd_n,
    int **restrict layerd_edge_separator,
    int **restrict layerd_edge,
    int **restrict layerd_belongs,

    int *restrict edge_size_buffer,
    int *restrict n_buffer_1,
    int *restrict n_buffer_2,
    int *restrict n_buffer_3
) {
    const int COLLAPSE_THRESHOLD = 100;

    layerd_n[shrink_depth] = n;
    layerd_edge_separator[shrink_depth] = (int*)edge_separator;
    layerd_edge[shrink_depth] = (int*)edge;
    if(edge_separator[n] == 0 || shrink_depth + 1 >= MAX_SHRINK_DEPTH) {
        return shrink_depth;
    }

    int *perm = n_buffer_1;
    generate_random_permutation(n, perm);
    int *belongs = memory_pool_allocalte(n);
    for(int i = 0; i < n; i++) {
        belongs[i] = -1;
    }
    int *centers = n_buffer_2;
    int *used = n_buffer_3;
    for(int i = 0; i < n; i++) {
        used[i] = -1;
    }
    int node_count = 0;
    int current = 0;
    for(int i = 0; i < n; i++) {
        const int p = perm[i];
        const int vwp = node_weight[p];
        if(belongs[p] < 0) {
            belongs[p] = node_count;
            centers[node_count] = p;
            if(edge_separator[p + 1] - edge_separator[p] <= COLLAPSE_THRESHOLD) {
                int edge_weight_sum = 0;
                for(int j = edge_separator[p]; j < edge_separator[p + 1]; j++) {
                    used[edge[j]] = edge_weight[j];
                    edge_weight_sum += edge_weight[j];
                }
                double max_score = -1.0;
                int max_score_node = -1;
                for(int j = edge_separator[p]; j < edge_separator[p + 1]; j++) {
                    const int t = edge[j];
                    if(belongs[t] < 0 && edge_separator[t + 1] - edge_separator[t] <= COLLAPSE_THRESHOLD) {
                        int common_edge_weight = 2 * edge_weight[j];
                        int not_common_edge_weight = edge_weight_sum - 2 * edge_weight[j];
                        for(int k = edge_separator[t]; k < edge_separator[t + 1]; k++) {
                            const int q = edge[k];
                            if(used[q] < 0) {
                                not_common_edge_weight += edge_weight[k];
                            } else {
                                not_common_edge_weight -= used[q];
                                common_edge_weight += edge_weight[k] + used[q];
                            }
                        }
                        const double score = (double)common_edge_weight / (common_edge_weight + not_common_edge_weight);
                        if(max_score < score) {
                            max_score = score;
                            max_score_node = t;
                        }
                    }
                }
                if(max_score_node >= 0) {
                    belongs[max_score_node] = node_count;
                }
                for(int j = edge_separator[p]; j < edge_separator[p + 1]; j++) {
                    used[edge[j]] = -1;
                }
            }
            node_count++;
        }
    }
    if(node_count < min_shrink_size || node_count == n) {
        memory_pool_deallocalte(n, belongs);
        return shrink_depth;
    }
    layerd_belongs[shrink_depth] = belongs;

    int edge_count = 0;
    for(int i = 0; i < node_count; i++) {
        const int c = centers[i];
        for(int j = edge_separator[c]; j < edge_separator[c + 1]; j++) {
            const int p = edge[j];
            const int bp = belongs[p];
            if(bp == i) {
                for(int k = edge_separator[p]; k < edge_separator[p + 1]; k++) {
                    const int t = edge[k];
                    const int bt = belongs[t];
                    if(bt != i && used[bt] < 0) {
                        used[bt] = 1;
                        edge_count++;
                    }
                }
            } else {
                if(used[bp] < 0) {
                    used[bp] = 1;
                    edge_count++;
                }
            }
        }
        for(int j = edge_separator[c]; j < edge_separator[c + 1]; j++) {
            const int p = edge[j];
            const int bp = belongs[p];
            if(bp == i) {
                for(int k = edge_separator[p]; k < edge_separator[p + 1]; k++) {
                    const int t = edge[k];
                    used[belongs[t]] = -1;
                }
            } else {
                used[bp] = -1;
            }
        }
    }

    int *sub_node_weight = n_buffer_1;
    int *sub_edge_separator = memory_pool_allocalte(node_count + 1);
    int *sub_edge = memory_pool_allocalte(edge_count);
    int *sub_edge_weight = edge_size_buffer;
    sub_edge_separator[0] = 0;
    for(int i = 0; i < node_count; i++) {
        const int c = centers[i];
        sub_edge_separator[i + 1] = sub_edge_separator[i];
        sub_node_weight[i] = node_weight[c];
        for(int j = edge_separator[c]; j < edge_separator[c + 1]; j++) {
            const int p = edge[j];
            const int bp = belongs[p];
            if(bp == i) {
                sub_node_weight[i] += node_weight[p];
                for(int k = edge_separator[p]; k < edge_separator[p + 1]; k++) {
                    const int t = edge[k];
                    const int bt = belongs[t];
                    if(bt != i) {
                        const int w = edge_weight[k];
                        if(used[bt] < 0) {
                            used[bt] = sub_edge_separator[i + 1];
                            sub_edge[sub_edge_separator[i + 1]] = bt;
                            sub_edge_weight[sub_edge_separator[i + 1]] = w;
                            sub_edge_separator[i + 1]++;
                        } else {
                            sub_edge_weight[used[bt]] += w;
                        }
                    }
                }
            } else {
                const int w = edge_weight[j];
                if(used[bp] < 0) {
                    used[bp] = sub_edge_separator[i + 1];
                    sub_edge[sub_edge_separator[i + 1]] = bp;
                    sub_edge_weight[sub_edge_separator[i + 1]] = w;
                    sub_edge_separator[i + 1]++;
                } else {
                    sub_edge_weight[used[bp]] += w;
                }
            }
        }
        for(int j = edge_separator[c]; j < edge_separator[c + 1]; j++) {
            const int p = edge[j];
            const int bp = belongs[p];
            if(bp == i) {
                for(int k = edge_separator[p]; k < edge_separator[p + 1]; k++) {
                    const int t = edge[k];
                    used[belongs[t]] = -1;
                }
            } else {
                used[bp] = -1;
            }
        }
    }

    return crhcm_shrink(
        min_shrink_size,
        shrink_depth + 1,
        node_count,
        sub_node_weight,
        sub_edge_separator,
        sub_edge,
        sub_edge_weight,
        layerd_n,
        layerd_edge_separator,
        layerd_edge,
        layerd_belongs,
        edge_weight,
        node_weight,
        n_buffer_2,
        n_buffer_3
    );
}

double crhcm_refinement_calculate_additional_cost(
    const weighted_range_sum_data *node_wrs,
    const weighted_range_sum_data *edge_wrs,
    const double node_cost_weight,
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict order_begin,
    const int *restrict order_end,
    const int *restrict belongs,
    const int additional_node,
    const int additional_meta
) {
    double node_additional_cost = 0.0;
    double edge_additional_cost = 0.0;
    for(int i = edge_separator[additional_node]; i < edge_separator[additional_node + 1]; i++) {
        const int p = edge[i];
        const int mp = belongs[p];
        if(mp == additional_meta) {
            continue;
        }
        int before = -1;
        int after = -1;
        if(mp < additional_meta) {
            before = mp;
        } else {
            after = mp;
        }
        bool affect = true;
        for(int j = edge_separator[p]; j < edge_separator[p + 1]; j++) {
            const int t = edge[j];
            if(t != additional_node) {
                const int mt = belongs[t];
                if(mt == additional_meta) {
                    affect = false;
                    break;
                }
                if(mt < additional_meta) {
                    if(before < mt) {
                        before = mt;
                    }
                } else {
                    if(after < 0 || after > mt) {
                        after = mt;
                    }
                }
            }
        }
        if(affect) {
            if(before >= 0 && after >= 0) {
                node_additional_cost +=
                    crhcm_calculate_cost(node_wrs, order_begin[before], order_end[before], order_begin[additional_meta], order_end[additional_meta])
                    + crhcm_calculate_cost(node_wrs, order_begin[additional_meta], order_end[additional_meta], order_begin[after], order_end[after])
                    - crhcm_calculate_cost(node_wrs, order_begin[before], order_end[before], order_begin[after], order_end[after]);
            } else if(before >= 0) {
                node_additional_cost += crhcm_calculate_cost(node_wrs, order_begin[before], order_end[before], order_begin[additional_meta], order_end[additional_meta]);
            } else if(after >= 0) {
                node_additional_cost += crhcm_calculate_cost(node_wrs, order_begin[additional_meta], order_end[additional_meta], order_begin[after], order_end[after]);
            }
        }
        if(additional_meta < mp) {
            edge_additional_cost += crhcm_calculate_cost(edge_wrs, order_begin[additional_meta], order_end[additional_meta], order_begin[mp], order_end[mp]);
        } else {
            edge_additional_cost += crhcm_calculate_cost(edge_wrs, order_begin[mp], order_end[mp], order_begin[additional_meta], order_end[additional_meta]);
        }
    }
    return node_cost_weight * node_additional_cost + (1.0 - node_cost_weight) * edge_additional_cost;
}

double crhcm_refinement_calculate_improve(
    const weighted_range_sum_data *node_wrs,
    const weighted_range_sum_data *edge_wrs,
    const double node_cost_weight,
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int *restrict order_begin,
    const int *restrict order_end,
    const int *restrict belongs,
    const int change_node,
    const int new_meta
) {
    double current_cost = crhcm_refinement_calculate_additional_cost(
        node_wrs,
        edge_wrs,
        node_cost_weight,
        n,
        edge_separator,
        edge,
        order_begin,
        order_end,
        belongs,
        change_node,
        belongs[change_node]
    );
    double new_cost = crhcm_refinement_calculate_additional_cost(
        node_wrs,
        edge_wrs,
        node_cost_weight,
        n,
        edge_separator,
        edge,
        order_begin,
        order_end,
        belongs,
        change_node,
        new_meta
    );
    return current_cost - new_cost;
}

int crhcm_refinement(
    const weighted_range_sum_data *node_wrs,
    const weighted_range_sum_data *edge_wrs,
    const double node_cost_weight,
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const int pairs_n,
    const int *restrict pairs,
    const int meta_size,
    int *restrict meta_separator,
    int *restrict meta,

    int *restrict n_buffer_1,
    int *restrict n_buffer_2,
    int *restrict n_buffer_3,
    int *restrict n_buffer_4,
    bool *restrict b_n_buffer
) {
    int *order_begin = n_buffer_1;
    int *order_end = n_buffer_2;
    int *belongs = n_buffer_3;
    int sum = 0;
    for(int i = 0; i < meta_size; i++) {
        order_begin[i] = sum;
        sum += meta_separator[i + 1] - meta_separator[i];
        order_end[i] = sum;
        for(int j = meta_separator[i]; j < meta_separator[i + 1]; j++) {
            belongs[meta[j]] = i;
        }
    }
    int *perm = n_buffer_4;
    generate_random_permutation(pairs_n, perm);
    bool *calculated = b_n_buffer;
    for(int i = 0; i < meta_size; i++) {
        calculated[i] = false;
    }
    for(int i = 0; i < pairs_n; i++) {
        const int p = perm[i];
        const int a = pairs[p * 2];
        const int b = pairs[p * 2 + 1];
        const int mp = belongs[a];
        double max_improve = 0.0;
        int max_improve_node = -1;
        int max_improve_meta = -1;
        calculated[mp] = true;
        for(int j = edge_separator[a]; j < edge_separator[a + 1]; j++) {
            const int t = belongs[edge[j]];
            if(mp != t && !calculated[t]) {
                calculated[t] = true;
                const double improve = crhcm_refinement_calculate_improve(
                    node_wrs,
                    edge_wrs,
                    node_cost_weight,
                    n,
                    edge_separator,
                    edge,
                    order_begin,
                    order_end,
                    belongs,
                    a,
                    t
                );
                if(max_improve < improve) {
                    max_improve = improve;
                    max_improve_node = a;
                    max_improve_meta = t;
                }
            }
        }
        for(int j = edge_separator[a]; j < edge_separator[a + 1]; j++) {
            calculated[belongs[edge[j]]] = false;
        }
        calculated[mp] = true;
        for(int j = edge_separator[b]; j < edge_separator[b + 1]; j++) {
            const int t = belongs[edge[j]];
            if(mp != t && !calculated[t]) {
                calculated[t] = true;
                const double improve = crhcm_refinement_calculate_improve(
                    node_wrs,
                    edge_wrs,
                    node_cost_weight,
                    n,
                    edge_separator,
                    edge,
                    order_begin,
                    order_end,
                    belongs,
                    b,
                    t
                );
                if(max_improve < improve) {
                    max_improve = improve;
                    max_improve_node = b;
                    max_improve_meta = t;
                }
            }
        }
        for(int j = edge_separator[b]; j < edge_separator[b + 1]; j++) {
            calculated[belongs[edge[j]]] = false;
        }
        calculated[mp] = false;
        if(max_improve_node >= 0) {
            belongs[max_improve_node] = max_improve_meta;
        }
    }
    int *count = n_buffer_1;
    for(int i = 0; i < meta_size; i++) {
        count[i] = 0;
    }
    for(int i = 0; i < n; i++) {
        count[belongs[i]]++;
    }
    meta_separator[0] = 0;
    meta_separator[1] = 0;
    for(int i = 0; i < meta_size - 1; i++) {
        meta_separator[i + 2] = meta_separator[i + 1] + count[i];
    }
    for(int i = 0; i < n; i++) {
        const int p = belongs[i];
        meta[meta_separator[p + 1]] = i;
        meta_separator[p + 1]++;
    }
    int new_meta_size = 0;
    for(int i = 0; i < meta_size; i++) {
        if(count[i] > 0) {
            new_meta_size++;
            meta_separator[new_meta_size] = meta_separator[i + 1];
        }
    }
    return new_meta_size;
}

void coarsen_refine_hybrid_cost_minimization_inner(
    const int min_shrink_size,
    const int node_cost_type,
    const int edge_cost_type,
    const double node_cost_weight,
    const double node_border_magnification,
    const double edge_border_magnification,
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
) {
    memory_pool_init();

    double *costs = (double*)malloc(sizeof(double) * (n + 10));
    for(int i = 0; i < n + 10; i++) {
        costs[i] = crhcm_calculate_single_cost(node_cost_type, i);
    }
    weighted_range_sum_data node_wrs;
    weighted_range_sum_init(&node_wrs, n + 10, costs);
    for(int i = 0; i < n + 10; i++) {
        costs[i] = crhcm_calculate_single_cost(edge_cost_type, i);
    }
    weighted_range_sum_data edge_wrs;
    weighted_range_sum_init(&edge_wrs, n + 10, costs);
    free(costs);

    int *order_begin = (int*)malloc(sizeof(int) * n);
    int *order_end = (int*)malloc(sizeof(int) * n);
    int meta_size = 1;
    int *meta_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *meta = (int*)malloc(sizeof(int) * n);

    int *new_meta_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *new_meta = (int*)malloc(sizeof(int) * n);

    int *inverse = (int*)malloc(sizeof(int) * n);
    int *adjacents = (int*)malloc(sizeof(int) * n);
    int *adjacent_inverse = (int*)malloc(sizeof(int) * n);
    for(int i = 0; i < n; i++) {
        adjacent_inverse[i] = -1;
    }
    int *sub_edge_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *sub_edge = (int*)malloc(sizeof(int) * edge_separator[n]);
    int *to_adjacent_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *to_adjacent = (int*)malloc(sizeof(int) * edge_separator[n]);
    int *from_adjacent_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *from_adjacent = (int*)malloc(sizeof(int) * edge_separator[n]);
    double *node_weight_first = (double*)malloc(sizeof(double) * n);
    double *node_weight_second = (double*)malloc(sizeof(double) * n);
    double *adjacent_node_weight_first = (double*)malloc(sizeof(double) * n);
    double *adjacent_node_weight_second = (double*)malloc(sizeof(double) * n);
    double *edge_cost_first = (double*)malloc(sizeof(double) * n);
    double *edge_cost_second = (double*)malloc(sizeof(double) * n);
    int *belongs = (int*)malloc(sizeof(int) * n);

    heap h1;
    heap_allocate(n, &h1);
    heap h2;
    heap_allocate(n, &h2);
    double *d_n_buffer_1 = (double*)malloc(sizeof(double) * n);
    double *d_n_buffer_2 = (double*)malloc(sizeof(double) * n);
    double *d_n_buffer_3 = (double*)malloc(sizeof(double) * n);
    double *d_n_buffer_4 = (double*)malloc(sizeof(double) * n);
    double *d_n_buffer_5 = (double*)malloc(sizeof(double) * n);
    int *_2n_buffer = (int*)malloc(sizeof(int) * 2 * n);
    int *n_buffer_1 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_2 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_3 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_4 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_5 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_6 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_7 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_8 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_9 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_10 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_11 = (int*)malloc(sizeof(int) * n);
    bool *b_n_buffer = (bool*)malloc(sizeof(bool) * n);

    int *node_weight = n_buffer_1;
    for(int i = 0; i < n; i++) {
        node_weight[i] = 1;
    }
    int *edge_weight_copy = (int*)malloc(sizeof(int) * edge_separator[n]);
    for(int i = 0; i < edge_separator[n]; i++) {
        edge_weight_copy[i] = 1;
    }
    int *edge_weight_buffer = (int*)malloc(sizeof(int) * edge_separator[n]);
    int *layerd_n = (int*)malloc(sizeof(int) * MAX_SHRINK_DEPTH);
    int **layerd_edge_separator = (int**)malloc(sizeof(int*) * MAX_SHRINK_DEPTH);
    int **layerd_edge = (int**)malloc(sizeof(int*) * MAX_SHRINK_DEPTH);
    int **layerd_belongs = (int**)malloc(sizeof(int*) * MAX_SHRINK_DEPTH);
    int current_layer = crhcm_shrink(
        min_shrink_size,
        0,
        n,
        node_weight,
        edge_separator,
        edge,
        edge_weight_copy,
        layerd_n,
        layerd_edge_separator,
        layerd_edge,
        layerd_belongs,
        edge_weight_buffer,
        n_buffer_2,
        n_buffer_3,
        n_buffer_4
    );
    free(edge_weight_copy);
    free(edge_weight_buffer);
    int *restrict current_edge_separator = layerd_edge_separator[current_layer];
    int *restrict current_edge = layerd_edge[current_layer];
    for(int i = 0; i < layerd_n[current_layer]; i++) {
        order_begin[i] = 0;
        order_end[i] = layerd_n[current_layer];
    }
    meta_separator[0] = 0;
    meta_separator[1] = layerd_n[current_layer];
    for(int i = 0; i < layerd_n[current_layer]; i++) {
        meta[i] = i;
    }

    int *tmp = (int*)malloc(sizeof(int) * 2 * n);

    while(meta_size > 0) {
        while(current_layer > 0 && layerd_n[current_layer] / meta_size < min_shrink_size) {
            const int pre_n = layerd_n[current_layer];
            memory_pool_deallocalte(current_edge_separator[pre_n], current_edge);
            memory_pool_deallocalte(pre_n + 1, current_edge_separator);
            current_layer--;

            const int current_n = layerd_n[current_layer];
            int *current_belongs = layerd_belongs[current_layer];

            int *count = n_buffer_1;
            for(int i = 0; i < pre_n; i++) {
                count[i] = 0;
            }
            for(int i = 0; i < current_n; i++) {
                count[current_belongs[i]]++;
            }

            int *separator = _2n_buffer;
            separator[0] = 0;
            separator[1] = 0;
            for(int i = 0; i < pre_n - 1; i++) {
                separator[i + 2] = separator[i + 1] + count[i];
            }
            int *belongs_inverse = n_buffer_2;
            for(int i = 0; i < current_n; i++) {
                const int b = current_belongs[i];
                belongs_inverse[separator[b + 1]] = i;
                separator[b + 1]++;
            }

            new_meta_separator[0] = 0;
            for(int i = 0; i < meta_size; i++) {
                new_meta_separator[i + 1] = new_meta_separator[i];
                for(int j = meta_separator[i]; j < meta_separator[i + 1]; j++) {
                    const int p = meta[j];
                    for(int k = separator[p]; k < separator[p + 1]; k++) {
                        new_meta[new_meta_separator[i + 1]] = belongs_inverse[k];
                        new_meta_separator[i + 1]++;
                    }
                }
            }

            int *tmp = meta_separator;
            meta_separator = new_meta_separator;
            new_meta_separator = tmp;
            tmp = meta;
            meta = new_meta;
            new_meta = tmp;

            memory_pool_deallocalte(current_n, current_belongs);

            current_edge_separator = layerd_edge_separator[current_layer];
            current_edge = layerd_edge[current_layer];

            int pairs_n = 0;
            int *pairs = n_buffer_3;
            for(int i = 0; i < pre_n; i++) {
                if(count[i] == 2) {
                    pairs[pairs_n * 2] = belongs_inverse[separator[i]];
                    pairs[pairs_n * 2 + 1] = belongs_inverse[separator[i] + 1];
                    pairs_n++;
                }
            }
            meta_size = crhcm_refinement(
                &node_wrs,
                &edge_wrs,
                node_cost_weight,
                current_n,
                current_edge_separator,
                current_edge,
                pairs_n,
                pairs,
                meta_size,
                meta_separator,
                meta,
                n_buffer_4,
                n_buffer_5,
                n_buffer_6,
                n_buffer_7,
                b_n_buffer
            );

            int sum = 0;
            for(int i = 0; i < meta_size; i++) {
                const int e = sum + meta_separator[i + 1] - meta_separator[i];
                for(int j = meta_separator[i]; j < meta_separator[i + 1]; j++) {
                    const int p = meta[j];
                    order_begin[p] = sum;
                    order_end[p] = e;
                }
                sum = e;
            }
        }

        int new_meta_size = 0;
        new_meta_separator[0] = 0;
        for(int i = 0; i < meta_size; i++) {
            const int b = meta_separator[i];
            const int e = meta_separator[i + 1];
            const int ob = order_begin[meta[b]];
            if(e - b == 1) {
                if(current_layer == 0) {
                    order[ob] = meta[b];
                } else {
                    new_meta_separator[new_meta_size + 1] = new_meta_separator[new_meta_size] + 1;
                    new_meta[new_meta_separator[new_meta_size]] = meta[b];
                    new_meta_size++;
                }
                continue;
            }
            const int oe = order_end[meta[b]];
            const int predicted_om = (oe + ob) / 2;

            const double node_cost_to_border = node_border_magnification * crhcm_calculate_cost(&node_wrs, ob, predicted_om, predicted_om, oe);
            const double edge_cost_to_border = edge_border_magnification * crhcm_calculate_cost(&edge_wrs, ob, predicted_om, predicted_om, oe);

            for(int j = 0; j < e - b; j++) {
                const int p = meta[b + j];
                inverse[p] = j;
            }
            int adjacent_size = 0;
            sub_edge_separator[0] = 0;
            to_adjacent_separator[0] = 0;
            for(int j = 0; j < e - b; j++) {
                const int p = meta[b + j];
                sub_edge_separator[j + 1] = sub_edge_separator[j];
                to_adjacent_separator[j + 1] = to_adjacent_separator[j];
                int adjacent_before_b = -1;
                int adjacent_before_e = -1;
                int adjacent_after_b = -1;
                int adjacent_after_e = -1;
                edge_cost_first[j] = 0.0;
                edge_cost_second[j] = 0.0;
                for(int k = current_edge_separator[p]; k < current_edge_separator[p + 1]; k++) {
                    const int t = current_edge[k];
                    const int tob = order_begin[t];
                    const int toe = order_end[t];
                    if(tob == ob) {
                        sub_edge[sub_edge_separator[j + 1]] = inverse[t];
                        sub_edge_separator[j + 1]++;
                    } else {
                        if(adjacent_inverse[t] < 0) {
                            adjacents[adjacent_size] = t;
                            adjacent_inverse[t] = adjacent_size;
                            adjacent_size++;
                        }
                        to_adjacent[to_adjacent_separator[j + 1]] = adjacent_inverse[t];
                        to_adjacent_separator[j + 1]++;

                        if(tob <= ob) {
                            if(adjacent_before_b < tob) {
                                adjacent_before_b = tob;
                                adjacent_before_e = toe;
                            }

                            edge_cost_first[j] += crhcm_calculate_cost(&edge_wrs, tob, toe, ob, predicted_om);
                            edge_cost_second[j] += crhcm_calculate_cost(&edge_wrs, tob, toe, predicted_om, oe);
                        } else {
                            if(adjacent_after_b < 0 || adjacent_after_b > tob) {
                                adjacent_after_b = tob;
                                adjacent_after_e = toe;
                            }

                            edge_cost_first[j] += crhcm_calculate_cost(&edge_wrs, ob, predicted_om, tob, toe);
                            edge_cost_second[j] += crhcm_calculate_cost(&edge_wrs, predicted_om, oe, tob, toe);
                        }
                    }
                }
                double b = node_cost_to_border;
                node_weight_first[j] = 0.0;
                node_weight_second[j] = 0.0;
                if(adjacent_before_b >= 0) {
                    b += crhcm_calculate_cost(&node_wrs, adjacent_before_b, adjacent_before_e, ob, predicted_om);
                    node_weight_first[j] -= crhcm_calculate_cost(&node_wrs, adjacent_before_b, adjacent_before_e, ob, predicted_om);
                    node_weight_second[j] -= crhcm_calculate_cost(&node_wrs, adjacent_before_b, adjacent_before_e, predicted_om, oe);
                }
                if(adjacent_after_b >= 0) {
                    b += crhcm_calculate_cost(&node_wrs, predicted_om, oe, adjacent_after_b, adjacent_after_e);
                    node_weight_first[j] -= crhcm_calculate_cost(&node_wrs, ob, predicted_om, adjacent_after_b, adjacent_after_e);
                    node_weight_second[j] -= crhcm_calculate_cost(&node_wrs, predicted_om, oe, adjacent_after_b, adjacent_after_e);
                }
                node_weight_first[j] += b;
                node_weight_second[j] += b;
            }
            from_adjacent_separator[0] = 0;
            for(int j = 0; j < adjacent_size; j++) {
                const int p = adjacents[j];
                adjacent_inverse[p] = -1;
                from_adjacent_separator[j + 1] = from_adjacent_separator[j];
                const int pob = order_begin[p];
                const int poe = order_end[p];
                int adjacent_before_b = -1;
                int adjacent_before_e = -1;
                int adjacent_after_b = -1;
                int adjacent_after_e = -1;
                if(pob <= ob) {
                    adjacent_before_b = pob;
                    adjacent_before_e = poe;
                } else {
                    adjacent_after_b = pob;
                    adjacent_after_e = poe;
                }
                for(int k = current_edge_separator[p]; k < current_edge_separator[p + 1]; k++) {
                    const int t = current_edge[k];
                    const int tob = order_begin[t];
                    const int toe = order_end[t];
                    const double expected_t = (toe - 1 + tob) / 2.0;
                    if(tob == ob) {
                        from_adjacent[from_adjacent_separator[j + 1]] = inverse[t];
                        from_adjacent_separator[j + 1]++;
                    } else {
                        if(tob <= ob) {
                            if(adjacent_before_b < tob) {
                                adjacent_before_b = tob;
                                adjacent_before_e = toe;
                            }
                        } else {
                            if(adjacent_after_b < 0 || adjacent_after_b > tob) {
                                adjacent_after_b = tob;
                                adjacent_after_e = toe;
                            }
                        }
                    }
                }
                double b = node_cost_to_border;
                adjacent_node_weight_first[j] = 0.0;
                adjacent_node_weight_second[j] = 0.0;
                if(adjacent_before_b >= 0) {
                    b += crhcm_calculate_cost(&node_wrs, adjacent_before_b, adjacent_before_e, ob, predicted_om);
                    adjacent_node_weight_first[j] -= crhcm_calculate_cost(&node_wrs, adjacent_before_b, adjacent_before_e, ob, predicted_om);
                    adjacent_node_weight_second[j] -= crhcm_calculate_cost(&node_wrs, adjacent_before_b, adjacent_before_e, predicted_om, oe);
                }
                if(adjacent_after_b >= 0) {
                    b += crhcm_calculate_cost(&node_wrs, predicted_om, oe, adjacent_after_b, adjacent_after_e);
                    adjacent_node_weight_first[j] -= crhcm_calculate_cost(&node_wrs, ob, predicted_om, adjacent_after_b, adjacent_after_e);
                    adjacent_node_weight_second[j] -= crhcm_calculate_cost(&node_wrs, predicted_om, oe, adjacent_after_b, adjacent_after_e);
                }
                adjacent_node_weight_first[j] += b;
                adjacent_node_weight_second[j] += b;
            }

            for(int j = 0; j < e - b; j++) {
                if(node_weight_first[j] < 0.0 || node_weight_second[j] < 0.0) {
                    if(node_weight_first[j] < node_weight_second[j]) {
                        node_weight_second[j] -= node_weight_first[j];
                        node_weight_first[j] = 0.0;
                    } else {
                        node_weight_first[j] -= node_weight_second[j];
                        node_weight_second[j] = 0.0;
                    }
                }

                node_weight_first[j] *= node_cost_weight;
                node_weight_second[j] *= node_cost_weight;

                edge_cost_first[j] *= (1.0 - node_cost_weight);
                edge_cost_second[j] *= (1.0 - node_cost_weight);
            }
            for(int j = 0; j < adjacent_size; j++) {
                if(adjacent_node_weight_first[j] < 0.0 || adjacent_node_weight_second[j] < 0.0) {
                    if(adjacent_node_weight_first[j] < adjacent_node_weight_second[j]) {
                        adjacent_node_weight_second[j] -= adjacent_node_weight_first[j];
                        adjacent_node_weight_first[j] = 0.0;
                    } else {
                        adjacent_node_weight_first[j] -= adjacent_node_weight_second[j];
                        adjacent_node_weight_second[j] = 0.0;
                    }
                }

                adjacent_node_weight_first[j] *= node_cost_weight;
                adjacent_node_weight_second[j] *= node_cost_weight;
            }

            int retry_number = (n - 1) / layerd_n[current_layer] + 1;
            if(retry_number > 512) {
                retry_number = 512;
            }
            crhcm_fiduccia_mattheyses(
                retry_number,
                e - b,
                adjacent_size,
                sub_edge_separator,
                sub_edge,
                to_adjacent_separator,
                to_adjacent,
                from_adjacent_separator,
                from_adjacent,
                node_weight_first,
                node_weight_second,
                adjacent_node_weight_first,
                adjacent_node_weight_second,
                edge_cost_first,
                edge_cost_to_border * (1.0 - node_cost_weight),
                edge_cost_second,
                belongs,
                &h1,
                &h2,
                d_n_buffer_1,
                d_n_buffer_2,
                d_n_buffer_3,
                d_n_buffer_4,
                d_n_buffer_5,
                _2n_buffer,
                n_buffer_1,
                n_buffer_2,
                n_buffer_3,
                n_buffer_4,
                n_buffer_5,
                n_buffer_6,
                n_buffer_7,
                n_buffer_8,
                n_buffer_9,
                n_buffer_10,
                n_buffer_11,
                b_n_buffer
            );

            int count_first = 0;
            int count_second = 0;
            for(int j = 0; j < e - b; j++) {
                if(belongs[j] == 0) {
                    count_first++;
                } else {
                    count_second++;
                }
            }
            assert(count_first > 0 && count_second > 0);
            const int om = ob + count_first;
            new_meta_separator[new_meta_size + 1] = new_meta_separator[new_meta_size];
            new_meta_separator[new_meta_size + 2] = new_meta_separator[new_meta_size] + count_first;
            for(int j = 0; j < e - b; j++) {
                const int p = meta[b + j];
                if(belongs[j] == 0) {
                    order_end[p] = om;
                    new_meta[new_meta_separator[new_meta_size + 1]] = p;
                    new_meta_separator[new_meta_size + 1]++;
                } else {
                    order_begin[p] = om;
                    new_meta[new_meta_separator[new_meta_size + 2]] = p;
                    new_meta_separator[new_meta_size + 2]++;
                }
            }
            new_meta_size += 2;
        }

        meta_size = new_meta_size;
        int *tmp = meta_separator;
        meta_separator = new_meta_separator;
        new_meta_separator = tmp;
        tmp = meta;
        meta = new_meta;
        new_meta = tmp;
    }

    weighted_range_sum_free(&node_wrs);
    weighted_range_sum_free(&edge_wrs);
    free(order_begin);
    free(order_end);
    free(meta_separator);
    free(meta);
    free(new_meta_separator);
    free(new_meta);
    free(inverse);
    free(adjacents);
    free(adjacent_inverse);
    free(sub_edge_separator);
    free(sub_edge);
    free(to_adjacent_separator);
    free(to_adjacent);
    free(from_adjacent_separator);
    free(from_adjacent);
    free(node_weight_first);
    free(node_weight_second);
    free(adjacent_node_weight_first);
    free(adjacent_node_weight_second);
    free(edge_cost_first);
    free(edge_cost_second);
    free(belongs);
    heap_free(&h1);
    heap_free(&h2);
    free(d_n_buffer_1);
    free(d_n_buffer_2);
    free(d_n_buffer_3);
    free(d_n_buffer_4);
    free(d_n_buffer_5);
    free(_2n_buffer);
    free(n_buffer_1);
    free(n_buffer_2);
    free(n_buffer_3);
    free(n_buffer_4);
    free(n_buffer_5);
    free(n_buffer_6);
    free(n_buffer_7);
    free(n_buffer_8);
    free(n_buffer_9);
    free(n_buffer_10);
    free(n_buffer_11);
    free(b_n_buffer);
    free(layerd_n);
    free(layerd_edge_separator);
    free(layerd_edge);
    free(layerd_belongs);
    memory_pool_free();
}

void coarsen_refine_hybrid_cost_minimization(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
) {
    coarsen_refine_hybrid_cost_minimization_inner(512, 2, 2, 0.5, 1.0, 1.0, n, edge_separator, edge, order);
}
