#include <string.h>
#include "lib.h"

void single_insertion_sort(
    int *restrict x,
    const int n
) {
    for(int i = 1; i < n; i++) {
        int t = x[i];
        for(int j = i - 1; j >= 0; j--) {
            if(x[j] > t) {
                x[j + 1] = x[j];
            } else {
                x[j + 1] = t;
                goto aa;
            }
        }
        x[0] = t;
aa:;
    }
}

void single_fast_radix_sort(
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
        for(int j = 0; j < n; j++) {
            const int t = (x[j] >> i) & mask;
            fast_radix_sort_ref[t][0] = x[j];
            fast_radix_sort_ref[t]++;
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

void single_slow_radix_sort(
    int *restrict x,
    const int n,
    int mask
) {
    int *p1 = x;
    int *p2 = x + n - 1;
    const int *e = x + n;
    while(1) {
        while((*p1 & mask) == 0 && p1 < e) {
            p1++;
        }
        while((*p2 & mask) != 0 && p2 > x) {
            p2--;
        }
        if(p1 >= p2) {
            break;
        }
        const int t = *p1;
        *p1 = *p2;
        *p2 = t;
        p1++;
        p2--;
    }
    if(mask == 1) {
        return;
    }
    const int m = p1 - x;
    if(m < (INT_BIT_LENGTH - 1) * RADIX_BIT_LENGTH && mask > (1 << (m / RADIX_BIT_LENGTH))) {
        single_insertion_sort(x, m);
    } else if(m <= SMALL_ARRAY_THRESHOLD) {
        single_fast_radix_sort(x, m, mask - 1);
    } else {
        single_slow_radix_sort(x, m, mask / 2);
    }
    const int k = n - m;
    if(k < (INT_BIT_LENGTH - 1) * RADIX_BIT_LENGTH && mask > (1 << (k / RADIX_BIT_LENGTH))) {
        single_insertion_sort(p1, k);
    } else if(k <= SMALL_ARRAY_THRESHOLD) {
        single_fast_radix_sort(p1, k, mask - 1);
    } else {
        single_slow_radix_sort(p1, k, mask / 2);
    }
}

void single_sort(
    const int l,
    int *restrict x,
    const int max_val
) {
    if(l < (INT_BIT_LENGTH - 1) * RADIX_BIT_LENGTH && max_val >= (1 << (l / RADIX_BIT_LENGTH))) {
        single_insertion_sort(x, l);
    } else if(l <= SMALL_ARRAY_THRESHOLD) {
        single_fast_radix_sort(x, l, max_val);
    } else {
        if(max_val == 0) {
            return;
        }
        int mask = 1;
        while(mask < max_val) {
            mask *= 2;
        }
        single_slow_radix_sort(x, l, mask);
    }
}
