#include <math.h>
#include <stdlib.h>
#include "lib.h"

#define SERIAL_THRESHOLD 5

int calculate_index(const int size, const int b, const int l) {
    const int d = l - b;
    return d * size - (d - 1) * d / 2 + b;
}

void range_sum_init(range_sum_data *restrict rs, const int n, const double *restrict data) {
    rs->raw = (double*)malloc(sizeof(double) * n);
    for(int i = 0; i < n; i++) {
        rs->raw[i] = data[i];
    }

    int pre_interval = n;
    for(int li = 0; li < LAYER_DEPTH; li++) {
        if(pre_interval < SERIAL_THRESHOLD) {
            rs->layers[li].interval = -1;
        } else {
            rs->layers[li].interval = sqrt(pre_interval);
            rs->layers[li].size = (n - 1) / rs->layers[li].interval;
            int maxd = (pre_interval - 2) / rs->layers[li].interval + 1;
            if(maxd > rs->layers[li].size - 1) {
                maxd = rs->layers[li].size - 1;
            }
            rs->layers[li].baby_up = (double*)malloc(sizeof(double) * n);
            rs->layers[li].baby_down = (double*)malloc(sizeof(double) * n);
            rs->layers[li].giant = (double*)malloc(sizeof(double) * calculate_index(rs->layers[li].size, 0, maxd));
            double sum = 0.0;
            for(int i = 0; i < n; i++) {
                if(i % rs->layers[li].interval == 0) {
                    sum = 0.0;
                }
                sum += data[i];
                rs->layers[li].baby_up[i] = sum;
            }
            sum = 0.0;
            for(int i = n - 1; i >= 0; i--) {
                sum += data[i];
                rs->layers[li].baby_down[i] = sum;
                if(i % rs->layers[li].interval == 0) {
                    sum = 0.0;
                }
            }
            for(int i = 0; i < rs->layers[li].size; i++) {
                rs->layers[li].giant[i] = rs->layers[li].baby_down[i * rs->layers[li].interval];
            }
            for(int d = 1; d < maxd; d++) {
                for(int i = 0; i + d < rs->layers[li].size; i++) {
                    const int m = i + d / 2;
                    const int p = calculate_index(rs->layers[li].size, i, i + d);
                    const int p1 = calculate_index(rs->layers[li].size, i, m);
                    const int p2 = calculate_index(rs->layers[li].size, m + 1, i + d);
                    rs->layers[li].giant[p] = rs->layers[li].giant[p1] + rs->layers[li].giant[p2];
                }
            }
            pre_interval = rs->layers[li].interval;
        }
    }
}

double range_sum(const range_sum_data *restrict rs, const int b, const int l) {
    if(b > l) {
        return 0;
    }
    for(int li = 0; li < LAYER_DEPTH; li++) {
        if(rs->layers[li].interval < 0) {
            break;
        } else {
            int bb = b / rs->layers[li].interval;
            int lb = l / rs->layers[li].interval;
            if(bb != lb) {
                double m = 0.0;
                if(bb + 1 <= lb - 1) {
                    m = rs->layers[li].giant[calculate_index(rs->layers[li].size, bb + 1, lb - 1)];
                }
                return (rs->layers[li].baby_down[b] + rs->layers[li].baby_up[l]) + m;
            }
        }
    }
    double sum = 0.0;
    for(int i = b; i <= l; i++) {
        sum += rs->raw[i];
    }
    return sum;
}

void range_sum_free(range_sum_data *restrict rs) {
    free(rs->raw);
    for(int li = 0; li < LAYER_DEPTH; li++) {
        if(rs->layers[li].interval >= 0) {
            free(rs->layers[li].baby_up);
            free(rs->layers[li].baby_down);
            free(rs->layers[li].giant);
        }
    }
}

void weighted_range_sum_init(weighted_range_sum_data *restrict wrs, const int n, const double *restrict data) {
    range_sum_init(&wrs->rs, n, data);

    int pre_interval = n;
    for(int li = 0; li < LAYER_DEPTH; li++) {
        if(pre_interval < SERIAL_THRESHOLD) {
            wrs->layers[li].interval = -1;
        } else {
            wrs->layers[li].interval = sqrt(pre_interval);
            wrs->layers[li].size = (n - 1) / wrs->layers[li].interval;
            int maxd = (pre_interval - 2) / wrs->layers[li].interval + 1;
            if(maxd > wrs->layers[li].size - 1) {
                maxd = wrs->layers[li].size - 1;
            }
            wrs->layers[li].baby_up = (double*)malloc(sizeof(double) * n);
            wrs->layers[li].baby_down = (double*)malloc(sizeof(double) * n);
            wrs->layers[li].giant = (double*)malloc(sizeof(double) * calculate_index(wrs->layers[li].size, 0, maxd));
            double sum = 0.0;
            int begin = 0;
            for(int i = 0; i < n; i++) {
                if(i % wrs->layers[li].interval == 0) {
                    sum = 0.0;
                    begin = i;
                }
                sum += (i - begin + 1) * data[i];
                wrs->layers[li].baby_up[i] = sum;
            }
            sum = 0.0;
            begin = n - 1;
            for(int i = n - 1; i >= 0; i--) {
                sum += range_sum(&wrs->rs, i, begin);
                wrs->layers[li].baby_down[i] = sum;
                if(i % wrs->layers[li].interval == 0) {
                    sum = 0.0;
                    begin = i - 1;
                }
            }
            for(int i = 0; i < wrs->layers[li].size; i++) {
                wrs->layers[li].giant[i] = wrs->layers[li].baby_down[i * wrs->layers[li].interval];
            }
            for(int d = 1; d < maxd; d++) {
                for(int i = 0; i + d < wrs->layers[li].size; i++) {
                    const int m = i + d / 2;
                    const int p = calculate_index(wrs->layers[li].size, i, i + d);
                    const int p1 = calculate_index(wrs->layers[li].size, i, m);
                    const int p2 = calculate_index(wrs->layers[li].size, m + 1, i + d);
                    wrs->layers[li].giant[p] =
                        range_sum(
                            &wrs->rs,
                            (m + 1) * wrs->layers[li].interval,
                            (i + d + 1) * wrs->layers[li].interval - 1
                        ) * (m - i + 1) * wrs->layers[li].interval
                        + wrs->layers[li].giant[p1]
                        + wrs->layers[li].giant[p2];
                }
            }
            pre_interval = wrs->layers[li].interval;
        }
    }
}

double weighted_range_sum(const weighted_range_sum_data *restrict wrs, const int b, const int l) {
    if(b > l) {
        return 0;
    }
    for(int li = 0; li < LAYER_DEPTH; li++) {
        if(wrs->layers[li].interval < 0) {
            break;
        } else {
            int bb = b / wrs->layers[li].interval;
            int lb = l / wrs->layers[li].interval;
            if(bb != lb) {
                double res =
                    (wrs->layers[li].baby_down[b]
                    + wrs->layers[li].baby_up[l])
                    + range_sum(
                        &wrs->rs,
                        (bb + 1) * wrs->layers[li].interval,
                        l
                    ) * ((bb + 1) * wrs->layers[li].interval - b);
                if(bb + 1 <= lb - 1) {
                    res =
                        (res
                        + range_sum(
                            &wrs->rs,
                            lb * wrs->layers[li].interval,
                            l
                        ) * (lb - bb - 1) * wrs->layers[li].interval)
                        + wrs->layers[li].giant[calculate_index(wrs->layers[li].size, bb + 1, lb - 1)];
                }
                return res;
            }
        }
    }
    double sum = 0.0;
    for(int i = b; i <= l; i++) {
        sum += (i - b + 1) * wrs->rs.raw[i];
    }
    return sum;
}

void weighted_range_sum_free(weighted_range_sum_data *restrict wrs) {
    range_sum_free(&wrs->rs);
    for(int li = 0; li < LAYER_DEPTH; li++) {
        if(wrs->layers[li].interval >= 0) {
            free(wrs->layers[li].baby_up);
            free(wrs->layers[li].baby_down);
            free(wrs->layers[li].giant);
        }
    }
}

double box_sum(const weighted_range_sum_data *restrict wrs, const int small_b, const int small_l, const int large_b, const int large_l) {
    if(small_l - small_b < large_l - large_b) {
        return
            weighted_range_sum(wrs, large_b - small_l, large_b - small_b)
            + (small_l - small_b + 1) * range_sum(&wrs->rs, large_b - small_b + 1, large_l - small_b)
            - weighted_range_sum(wrs, large_l - small_l + 1, large_l - small_b);
    } else {
        return
            weighted_range_sum(wrs, large_b - small_l, large_l - small_l)
            + (large_l - large_b + 1) * range_sum(&wrs->rs, large_l - small_l + 1, large_l - small_b)
            - weighted_range_sum(wrs, large_b - small_b + 1, large_l - small_b);
    }
}

double triangle_sum(const weighted_range_sum_data *restrict wrs, const int top) {
    return top * range_sum(&wrs->rs, 1, top) - weighted_range_sum(wrs, 2, top);
}
