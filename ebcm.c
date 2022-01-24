#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"
#include "reordering.h"

double ebcm_calculate_separation_cost(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const double *restrict cost_at_first,
    const double *restrict cost_at_second,
    const double cost_of_internal_edge,
    const int *restrict belongs
) {
    double sum = 0.0;
    for(int i = 0; i < n; i++) {
        if(belongs[i] == 0) {
            sum += cost_at_first[i];
        } else {
            sum += cost_at_second[i];
        }
        for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
            const int p = edge[j];
            if(i < p && belongs[i] != belongs[p]) {
                sum += cost_of_internal_edge;
            }
        }
    }
    return sum;
}

void ebcm_first_separation(
    const int n,
    const double *restrict cost_at_first,
    const double *restrict cost_at_second,
    int *restrict belongs,

    int *restrict _2n_buffer
) {
    const int MAX_WEIGHT = 10000;

    double min_val = cost_at_first[0] - cost_at_second[0];
    double max_val = min_val;
    for(int i = 1; i < n; i++) {
        const double d = cost_at_first[i] - cost_at_second[i];
        if(min_val > d) {
            min_val = d;
        }
        if(max_val < d) {
            max_val = d;
        }
    }
    if(max_val - min_val < 0.001) {
        for(int i = 0; i < n / 2; i++) {
            belongs[i] = 0;
        }
        for(int i = n / 2; i < n; i++) {
            belongs[i] = 1;
        }
        return;
    }
    int *tmp = _2n_buffer;
    for(int i = 0; i < n; i++) {
        const int d = (cost_at_first[i] - cost_at_second[i] - min_val) * MAX_WEIGHT / (max_val - min_val);
        tmp[i * 2] = d;
        tmp[i * 2 + 1] = i;
    }

    first_element_sort(n, tmp, MAX_WEIGHT);

    for(int i = 0; i < n / 2; i++) {
        belongs[tmp[i * 2 + 1]] = 0;
    }
    for(int i = n / 2; i < n; i++) {
        belongs[tmp[i * 2 + 1]] = 1;
    }
}

double ebcm_fiduccia_mattheyses_inner(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const double *restrict cost_at_first,
    const double *restrict cost_at_second,
    const double cost_of_internal_edge,
    const double error_tolerance,
    int *restrict belongs,

    heap *restrict h1,
    heap *restrict h2,
    double *restrict d_n_buffer_1,
    double *restrict d_n_buffer_2,
    double *restrict d_n_buffer_3,
    int *restrict n_buffer_1,
    int *restrict n_buffer_2,
    int *restrict n_buffer_3,
    int *restrict n_buffer_4,
    bool *restrict b_n_buffer
) {
    const double TORELANCE = 1.3;

    heap *h_first = h1;
    heap *h_second = h2;
    int *bucket_first = n_buffer_1;
    int *bucket_second = n_buffer_2;
    double *bucket_gain_first = d_n_buffer_1;
    double *bucket_gain_second = d_n_buffer_2;
    int *inverse = n_buffer_3;
    int count_first = 0;
    int count_second = 0;
    for(int i = 0; i < n; i++) {
        if(belongs[i] == 0) {
            inverse[i] = count_first;
            bucket_first[count_first] = i;
            bucket_gain_first[count_first] = cost_at_first[i] - cost_at_second[i];
            for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
                if(belongs[edge[j]] == 0) {
                    bucket_gain_first[count_first] -= cost_of_internal_edge;
                } else {
                    bucket_gain_first[count_first] += cost_of_internal_edge;
                }
            }
            count_first++;
        } else {
            inverse[i] = count_second;
            bucket_second[count_second] = i;
            bucket_gain_second[count_second] = cost_at_second[i] - cost_at_first[i];
            for(int j = edge_separator[i]; j < edge_separator[i + 1]; j++) {
                if(belongs[edge[j]] == 0) {
                    bucket_gain_second[count_second] += cost_of_internal_edge;
                } else {
                    bucket_gain_second[count_second] -= cost_of_internal_edge;
                }
            }
            count_second++;
        }
    }
    heap_initialize(count_first, bucket_gain_first, h_first);
    heap_initialize(count_second, bucket_gain_second, h_second);

    bool *used = b_n_buffer;
    for(int i = 0; i < n; i++) {
        used[i] = false;
    }
    int moved_size = 0;
    int *moved = n_buffer_4;
    double *moved_gain = d_n_buffer_3;
    int current_count_first = count_first;
    int current_count_second = count_second;
    while(1) {
        const bool c_first = !heap_empty(h_first);
        const bool c_second = !heap_empty(h_second);
        int m;
        if(c_first && c_second) {
            if(current_count_first > current_count_second && (current_count_first + 1) > (current_count_second - 1) * TORELANCE) {
                m = 0;
            } else if(current_count_second > current_count_first && (current_count_second + 1) > (current_count_first - 1) * TORELANCE) {
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

        if(m == 0) {
            const int p = bucket_first[heap_front(h_first).pos];
            moved[moved_size] = p;
            moved_gain[moved_size] = heap_front(h_first).data;
            heap_pop(h_first);
            used[p] = true;
            for(int i = edge_separator[p]; i < edge_separator[p + 1]; i++) {
                const int t = edge[i];
                if(!used[t]) {
                    if(belongs[t] == 0) {
                        heap_update(inverse[t], 2.0 * cost_of_internal_edge, h_first);
                    } else {
                        heap_update(inverse[t], -2.0 * cost_of_internal_edge, h_second);
                    }
                }
            }
            current_count_first--;
            current_count_second++;
        } else {
            const int p = bucket_second[heap_front(h_second).pos];
            moved[moved_size] = p;
            moved_gain[moved_size] = heap_front(h_second).data;
            heap_pop(h_second);
            used[p] = true;
            for(int i = edge_separator[p]; i < edge_separator[p + 1]; i++) {
                const int t = edge[i];
                if(!used[t]) {
                    if(belongs[t] == 0) {
                        heap_update(inverse[t], -2.0 * cost_of_internal_edge, h_first);
                    } else {
                        heap_update(inverse[t], 2.0 * cost_of_internal_edge, h_second);
                    }
                }
            }
            current_count_first++;
            current_count_second--;
        }
        moved_size++;
    }

    double max_val = 0.0;
    int max_val_pos = -1;
    double current_gain = 0.0;
    current_count_first = count_first;
    current_count_second = count_second;
    for(int i = 0; i < moved_size; i++) {
        if(belongs[moved[i]] == 0) {
            current_count_first--;
            current_count_second++;
        } else {
            current_count_first++;
            current_count_second--;
        }
        current_gain += moved_gain[i];
        if((current_count_first - current_count_second > 1 && current_count_first > current_count_second * TORELANCE) || (current_count_second - current_count_first > 1 && current_count_second > current_count_first * TORELANCE)) {
            continue;
        }
        if(max_val + error_tolerance < current_gain) {
            max_val = current_gain;
            max_val_pos = i;
        }
    }
    if(max_val_pos < 0) {
        return -2.0;
    }
    for(int i = 0; i <= max_val_pos; i++) {
        const int p = moved[i];
        if(belongs[p] == 0) {
            belongs[p] = 1;
        } else {
            belongs[p] = 0;
        }
    }
    return max_val;
}

void ebcm_fiduccia_mattheyses(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    const double *restrict cost_at_first,
    const double *restrict cost_at_second,
    const double cost_of_internal_edge,
    int *restrict belongs,

    heap *restrict h1,
    heap *restrict h2,
    double *restrict d_n_buffer_1,
    double *restrict d_n_buffer_2,
    double *restrict d_n_buffer_3,
    int *restrict _2n_buffer,
    int *restrict n_buffer_1,
    int *restrict n_buffer_2,
    int *restrict n_buffer_3,
    bool *restrict b_n_buffer
) {
    ebcm_first_separation(
        n,
        cost_at_first,
        cost_at_second,
        belongs,
        _2n_buffer
    );
    double cost = ebcm_calculate_separation_cost(
        n,
        edge_separator,
        edge,
        cost_at_first,
        cost_at_second,
        cost_of_internal_edge,
        belongs
    );
    while(1) {
        double error_tolerance = cost / pow(10.0, 9.0);
        if(error_tolerance < 0.000001) {
            error_tolerance = 0.000001;
        }
        double improve = ebcm_fiduccia_mattheyses_inner(
            n,
            edge_separator,
            edge,
            cost_at_first,
            cost_at_second,
            cost_of_internal_edge,
            error_tolerance,
            belongs,
            h1,
            h2,
            d_n_buffer_1,
            d_n_buffer_2,
            d_n_buffer_3,
            _2n_buffer,
            n_buffer_1,
            n_buffer_2,
            n_buffer_3,
            b_n_buffer
        );
        if(improve < -1.0) {
            break;
        }
        cost -= improve;
    }
}

double ebcm_calculate_cost(
    const weighted_range_sum_data *wrs,
    int x_order_begin,
    int x_order_end,
    int y_order_begin,
    int y_order_end
) {
    if(x_order_begin > y_order_begin) {
        int tmp = x_order_begin;
        x_order_begin = y_order_begin;
        y_order_begin = tmp;
        tmp = x_order_end;
        x_order_end = y_order_end;
        y_order_end = tmp;
    }
    return box_sum(wrs, x_order_begin, x_order_end - 1, y_order_begin, y_order_end - 1) / (x_order_end - x_order_begin) / (y_order_end - y_order_begin);
}

void edge_based_cost_minimization_inner(
    const int cost_type,
    const double border_magnification,
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
) {
    double *costs = (double*)malloc(sizeof(double) * (n + 10));
    for(int i = 0; i < n + 10; i++) {
        if(i == 0) {
            costs[i] = 0;
        } else {
            switch (cost_type) {
            case 0:
                costs[i] = i;
                break;
            case 1:
                costs[i] = sqrt(i);
                break;
            case 2:
                costs[i] = log2(i);
                break;
            case 3:
                costs[i] = sqrt(log2(i));
                break;
            case 4:
                costs[i] = pow(log2(i), 0.25);
                break;
            }
        }
    }
    weighted_range_sum_data wrs;
    weighted_range_sum_init(&wrs, n + 10, costs);
    free(costs);

    int *order_begin = (int*)malloc(sizeof(int) * n);
    int *order_end = (int*)malloc(sizeof(int) * n);
    for(int i = 0; i < n; i++) {
        order_begin[i] = 0;
        order_end[i] = n;
    }

    int meta_size = 1;
    int *meta_separator = (int*)malloc(sizeof(int) * (n + 1));
    meta_separator[0] = 0;
    meta_separator[1] = n;
    int *meta = (int*)malloc(sizeof(int) * n);
    for(int i = 0; i < n; i++) {
        meta[i] = i;
    }

    int *new_meta_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *new_meta = (int*)malloc(sizeof(int) * n);

    int *inverse = (int*)malloc(sizeof(int) * n);
    int *sub_edge_separator = (int*)malloc(sizeof(int) * (n + 1));
    int *sub_edge = (int*)malloc(sizeof(int) * edge_separator[n]);
    double *cost_at_first = (double*)malloc(sizeof(double) * n);
    double *cost_at_second = (double*)malloc(sizeof(double) * n);
    int *belongs = (int*)malloc(sizeof(int) * n);

    heap h1;
    heap_allocate(n, &h1);
    heap h2;
    heap_allocate(n, &h2);
    double *d_n_buffer_1 = (double*)malloc(sizeof(double) * n);
    double *d_n_buffer_2 = (double*)malloc(sizeof(double) * n);
    double *d_n_buffer_3 = (double*)malloc(sizeof(double) * n);
    int *_2n_buffer = (int*)malloc(sizeof(int) * 2 * n);
    int *n_buffer_1 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_2 = (int*)malloc(sizeof(int) * n);
    int *n_buffer_3 = (int*)malloc(sizeof(int) * n);
    bool *b_n_buffer = (bool*)malloc(sizeof(bool) * n);
    while(meta_size > 0) {
        int new_meta_size = 0;
        new_meta_separator[0] = 0;
        for(int i = 0; i < meta_size; i++) {
            const int b = meta_separator[i];
            const int e = meta_separator[i + 1];
            const int ob = order_begin[meta[b]];
            if(e - b == 1) {
                order[ob] = meta[b];
                continue;
            }
            const int oe = order_end[meta[b]];
            const int predicted_om = (oe + ob) / 2;

            for(int j = 0; j < e - b; j++) {
                const int p = meta[b + j];
                inverse[p] = j;
            }
            sub_edge_separator[0] = 0;
            for(int j = 0; j < e - b; j++) {
                const int p = meta[b + j];
                sub_edge_separator[j + 1] = sub_edge_separator[j];
                cost_at_first[j] = 0.0;
                cost_at_second[j] = 0.0;
                for(int k = edge_separator[p]; k < edge_separator[p + 1]; k++) {
                    const int t = edge[k];
                    const int tob = order_begin[t];
                    if(tob == ob) {
                        sub_edge[sub_edge_separator[j + 1]] = inverse[t];
                        sub_edge_separator[j + 1]++;
                    } else {
                        const int toe = order_end[t];
                        cost_at_first[j] += ebcm_calculate_cost(&wrs, ob, predicted_om, tob, toe);
                        cost_at_second[j] += ebcm_calculate_cost(&wrs, predicted_om, oe, tob, toe);
                    }
                }
            }

            const double cost_of_internal_edge = border_magnification * ebcm_calculate_cost(&wrs, ob, predicted_om, predicted_om, oe);
            ebcm_fiduccia_mattheyses(
                e - b,
                sub_edge_separator,
                sub_edge,
                cost_at_first,
                cost_at_second,
                cost_of_internal_edge,
                belongs,
                &h1,
                &h2,
                d_n_buffer_1,
                d_n_buffer_2,
                d_n_buffer_3,
                _2n_buffer,
                n_buffer_1,
                n_buffer_2,
                n_buffer_3,
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

    weighted_range_sum_free(&wrs);
    free(order_begin);
    free(order_end);
    free(meta_separator);
    free(meta);
    free(new_meta_separator);
    free(new_meta);
    free(inverse);
    free(sub_edge_separator);
    free(sub_edge);
    free(cost_at_first);
    free(cost_at_second);
    free(belongs);
    heap_free(&h1);
    heap_free(&h2);
    free(d_n_buffer_1);
    free(d_n_buffer_2);
    free(d_n_buffer_3);
    free(_2n_buffer);
    free(n_buffer_1);
    free(n_buffer_2);
    free(n_buffer_3);
    free(b_n_buffer);
}

void edge_based_cost_minimization(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
) {
    edge_based_cost_minimization_inner(2, 1.0, n, edge_separator, edge, order);
}
