#include <string.h>
#include "lib.h"

int fast_radix_sort_tmp[(1 << RADIX_BIT_LENGTH) * SMALL_ARRAY_THRESHOLD * 2];
int *restrict fast_radix_sort_ref[1 << RADIX_BIT_LENGTH];

void first_element_insertion_sort(
    int *restrict x,
    const int n
) {
    for(int i = 2; i < n; i += 2) {
        int t = x[i];
        int v = x[i + 1];
        for(int j = i - 2; j >= 0; j -= 2) {
            if(x[j] > t) {
                x[j + 2] = x[j];
                x[j + 3] = x[j + 1];
            } else {
                x[j + 2] = t;
                x[j + 3] = v;
                goto aa;
            }
        }
        x[0] = t;
        x[1] = v;
aa:;
    }
}

void first_element_fast_radix_sort(
    int *restrict x,
    const int n,
    const int max_val
) {
    const int mask = (1 << RADIX_BIT_LENGTH) - 1;
    for(int i = 0; i < INT_BIT_LENGTH; i += RADIX_BIT_LENGTH) {
        if(max_val < (1 << i)) {
            break;
        }
        for(int j = 0; j < (1 << RADIX_BIT_LENGTH); j++) {
            fast_radix_sort_ref[j] = fast_radix_sort_tmp + j * n;
        }
        for(int j = 0; j < n; j += 2) {
            const int t = (x[j] >> i) & mask;
            fast_radix_sort_ref[t][0] = x[j];
            fast_radix_sort_ref[t][1] = x[j + 1];
            fast_radix_sort_ref[t] += 2;
        }
        int count = 0;
        for(int j = 0; j < (1 << RADIX_BIT_LENGTH); j++) {
            const int *b = fast_radix_sort_tmp + j * n;
            const int s = fast_radix_sort_ref[j] - b;
            memcpy(x + count, b, sizeof(int) * s);
            count += s;
        }
    }
}

void first_element_slow_radix_sort(
    int *restrict x,
    const int n,
    int mask
) {
    int *p1 = x;
    int *p2 = x + n - 2;
    const int *e = x + n;
    while(1) {
        while((*p1 & mask) == 0 && p1 < e) {
            p1 += 2;
        }
        while((*p2 & mask) != 0 && p2 > x) {
            p2 -= 2;
        }
        if(p1 >= p2) {
            break;
        }
        int t = p1[0];
        p1[0] = p2[0];
        p2[0] = t;
        t = p1[1];
        p1[1] = p2[1];
        p2[1] = t;
        p1 += 2;
        p2 -= 2;
    }
    if(mask == 1) {
        return;
    }
    const int m = p1 - x;
    if(m < (INT_BIT_LENGTH - 1) * RADIX_BIT_LENGTH * 2 && mask > (1 << (m / (RADIX_BIT_LENGTH * 2)))) {
        first_element_insertion_sort(x, m);
    } else if(m <= SMALL_ARRAY_THRESHOLD * 2) {
        first_element_fast_radix_sort(x, m, mask - 1);
    } else {
        first_element_slow_radix_sort(x, m, mask / 2);
    }
    const int k = n - m;
    if(k < (INT_BIT_LENGTH - 1) * RADIX_BIT_LENGTH * 2 && mask > (1 << (k / (RADIX_BIT_LENGTH * 2)))) {
        first_element_insertion_sort(p1, k);
    } else if(k <= SMALL_ARRAY_THRESHOLD * 2) {
        first_element_fast_radix_sort(p1, k, mask - 1);
    } else {
        first_element_slow_radix_sort(p1, k, mask / 2);
    }
}

void first_element_sort(
    const int l,
    int *restrict x,
    const int max_val
) {
    if(l < (INT_BIT_LENGTH - 1) * RADIX_BIT_LENGTH && max_val >= (1 << (l / RADIX_BIT_LENGTH))) {
        first_element_insertion_sort(x, l * 2);
    } else if(l <= SMALL_ARRAY_THRESHOLD) {
        first_element_fast_radix_sort(x, l * 2, max_val);
    } else {
        if(max_val == 0) {
            return;
        }
        int mask = 1;
        while(mask < max_val) {
            mask *= 2;
        }
        first_element_slow_radix_sort(x, l * 2, mask);
    }
}
