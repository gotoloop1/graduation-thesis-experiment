#include <stdlib.h>
#include "lib.h"

int *memory_pool[50][1000];
int memory_pool_rem[50];

void memory_pool_init() {
    for(int i = 0; i < 50; i++) {
        memory_pool_rem[i] = 0;
    }
}

int *memory_pool_allocalte(int size) {
    int cls = 0;
    while(size > (1 << cls)) {
        cls++;
    }
    if(memory_pool_rem[cls] == 0) {
        return (int*)malloc(sizeof(int) * (1 << cls));
    } else {
        int *ret = memory_pool[cls][memory_pool_rem[cls] - 1];
        memory_pool_rem[cls]--;
        return ret;
    }
}

void memory_pool_deallocalte(int size, int *p) {
    int cls = 0;
    while(size > (1 << cls)) {
        cls++;
    }
    memory_pool[cls][memory_pool_rem[cls]] = p;
    memory_pool_rem[cls]++;
}

void memory_pool_free() {
    for(int i = 0; i < 50; i++) {
        for(int j = 0; j < memory_pool_rem[i]; j++) {
            free(memory_pool[i][j]);
        }
        memory_pool_rem[i] = 0;
    }
}

double *double_memory_pool[50][1000];
int double_memory_pool_rem[50];

void double_memory_pool_init() {
    for(int i = 0; i < 50; i++) {
        double_memory_pool_rem[i] = 0;
    }
}

double *double_memory_pool_allocalte(int size) {
    int cls = 0;
    while(size > (1 << cls)) {
        cls++;
    }
    if(double_memory_pool_rem[cls] == 0) {
        return (double*)malloc(sizeof(double) * (1 << cls));
    } else {
        double *ret = double_memory_pool[cls][double_memory_pool_rem[cls] - 1];
        double_memory_pool_rem[cls]--;
        return ret;
    }
}

void double_memory_pool_deallocalte(int size, double *p) {
    int cls = 0;
    while(size > (1 << cls)) {
        cls++;
    }
    double_memory_pool[cls][double_memory_pool_rem[cls]] = p;
    double_memory_pool_rem[cls]++;
}

void double_memory_pool_free() {
    for(int i = 0; i < 50; i++) {
        for(int j = 0; j < double_memory_pool_rem[i]; j++) {
            free(double_memory_pool[i][j]);
        }
        double_memory_pool_rem[i] = 0;
    }
}
