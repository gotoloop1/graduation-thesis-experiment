#include <stdint.h>
#include "lib.h"

uint32_t xor128_x = 123456789;
uint32_t xor128_y = 362436069;
uint32_t xor128_z = 521288629;
uint32_t xor128_w = 88675123;

void xor128_set_seed(unsigned int w) {
    xor128_w = w;
    for(int i = 0; i < 1000; i++) {
        xor128();
    }
}

unsigned int xor128() {
    uint32_t t = xor128_x ^ (xor128_x << 11);
    xor128_x = xor128_y; xor128_y = xor128_z; xor128_z = xor128_w;
    return xor128_w = (xor128_w ^ (xor128_w >> 19)) ^ (t ^ (t >> 8)); 
}

void generate_random_permutation(
    const int n,
    int *restrict perm
) {
    for(int i = 0; i < n; i++) {
        perm[i] = i;
    }
    for(int i = n - 1; i >= 1; i--) {
        const int j = xor128() % (i + 1);
        const int t = perm[i];
        perm[i] = perm[j];
        perm[j] = t;
    }
}
