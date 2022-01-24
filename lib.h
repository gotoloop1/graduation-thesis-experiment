#ifndef LIB_INCLUDED
#define LIB_INCLUDED

#include <stdbool.h>

#define INT_BIT_LENGTH 32
#define RADIX_BIT_LENGTH 4
#define SMALL_ARRAY_THRESHOLD (1 << 16)
extern int fast_radix_sort_tmp[(1 << RADIX_BIT_LENGTH) * SMALL_ARRAY_THRESHOLD * 2];
extern int *restrict fast_radix_sort_ref[1 << RADIX_BIT_LENGTH];
void single_sort(
    const int l,
    int *restrict x,
    const int max_val
);
void first_element_sort(
    const int l,
    int *restrict x,
    const int max_val
);

typedef struct {
    double data;
    int pos;
} heap_element;
typedef struct {
    int n;
    heap_element *elements;
    int *pos;
    bool *used;
} heap;
void heap_allocate(const int size, heap *restrict h);
void heap_initialize(const int n, const double *restrict data, heap *restrict h);
bool heap_empty(const heap *restrict h);
heap_element heap_front(const heap *restrict h);
void heap_pop(heap *restrict h);
void heap_update(const int pos, const double add, heap *restrict h);
void heap_free(heap *restrict h);
void dump_heap(const heap *restrict h, const int *restrict map);

#define LAYER_DEPTH 3
typedef struct {
    double *restrict raw;
    struct {
        int interval;
        int size;
        double *restrict baby_up;
        double *restrict baby_down;
        double *restrict giant;
    } layers[LAYER_DEPTH];
} range_sum_data;
typedef struct {
    range_sum_data rs;
    struct {
        int interval;
        int size;
        double *restrict baby_up;
        double *restrict baby_down;
        double *restrict giant;
    } layers[LAYER_DEPTH];
} weighted_range_sum_data;
void range_sum_init(range_sum_data *restrict rs, const int n, const double *restrict d);
double range_sum(const range_sum_data *restrict rs, const int b, const int l);
void range_sum_free(range_sum_data *restrict rs);
void weighted_range_sum_init(weighted_range_sum_data *restrict wrs, const int n, const double *restrict data);
double weighted_range_sum(const weighted_range_sum_data *restrict wrs, const int b, const int l);
void weighted_range_sum_free(weighted_range_sum_data *restrict wrs);
double box_sum(const weighted_range_sum_data *restrict wrs, const int small_b, const int small_l, const int large_b, const int large_l);
double triangle_sum(const weighted_range_sum_data *restrict wrs, const int top);

void memory_pool_init();
int *memory_pool_allocalte(int size);
void memory_pool_deallocalte(int size, int *p);
void memory_pool_free();
void double_memory_pool_init();
double *double_memory_pool_allocalte(int size);
void double_memory_pool_deallocalte(int size, double *p);
void double_memory_pool_free();

void xor128_set_seed(unsigned int w);
unsigned int xor128();
void generate_random_permutation(const int n, int *restrict perm);

#endif
