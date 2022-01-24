#include <stdlib.h>
#include "lib.h"

#include <stdio.h>

void heap_allocate(const int size, heap *restrict h) {
    h->elements = (heap_element*)malloc(sizeof(heap_element) * size);
    h->pos = (int*)malloc(sizeof(int) * size);
    h->used = (bool*)malloc(sizeof(bool) * size);
}

void heap_update_up(int i, heap *restrict h) {
    while(i > 0) {
        const int p = (i - 1) / 2;
        if(h->elements[p].data > h->elements[i].data) {
            break;
        }
        heap_element tmp = h->elements[p];
        h->elements[p] = h->elements[i];
        h->elements[i] = tmp;
        h->pos[h->elements[i].pos] = i;
        i = p;
    }
    h->pos[h->elements[i].pos] = i;
}

void heap_update_down(int i, heap *restrict h) {
    while(1) {
        const int p1 = i * 2 + 1;
        const int p2 = i * 2 + 2;
        const bool c1 = p1 < h->n && h->used[p1];
        const bool c2 = p2 < h->n && h->used[p2];
        int m;
        if(c1 && c2) {
            if(h->elements[p1].data > h->elements[p2].data) {
                m = p1;
            } else {
                m = p2;
            }
        } else if(c1) {
            m = p1;
        } else if(c2) {
            m = p2;
        } else {
            break;
        }
        if(h->elements[i].data > h->elements[m].data) {
            break;
        }
        heap_element tmp = h->elements[m];
        h->elements[m] = h->elements[i];
        h->elements[i] = tmp;
        h->pos[h->elements[i].pos] = i;
        i = m;
    }
    h->pos[h->elements[i].pos] = i;
}

void heap_update_vacuum(int i, heap *restrict h) {
    while(1) {
        const int p1 = i * 2 + 1;
        const int p2 = i * 2 + 2;
        const bool c1 = p1 < h->n && h->used[p1];
        const bool c2 = p2 < h->n && h->used[p2];
        int m;
        if(c1 && c2) {
            if(h->elements[p1].data > h->elements[p2].data) {
                m = p1;
            } else {
                m = p2;
            }
        } else if(c1) {
            m = p1;
        } else if(c2) {
            m = p2;
        } else {
            break;
        }
        h->elements[i] = h->elements[m];
        h->pos[h->elements[i].pos] = i;
        i = m;
    }
    h->used[i] = false;
}

void heap_initialize(const int n, const double *restrict data, heap *restrict h) {
    h->n = n;
    for(int i = 0; i < n; i++) {
        heap_element he = {data[i], i};
        h->elements[i] = he;
        h->pos[i] = i;
        h->used[i] = true;
        heap_update_up(i, h);
    }
}

bool heap_empty(const heap *restrict h) {
    return !h->used[0];
}

heap_element heap_front(const heap *restrict h) {
    return h->elements[0];
}

void heap_pop(heap *restrict h) {
    heap_update_vacuum(0, h);
}

void heap_update(const int pos, const double add, heap *restrict h) {
    const int p = h->pos[pos];
    h->elements[p].data += add;
    if(add > 0) {
        heap_update_up(p, h);
    } else {
        heap_update_down(p, h);
    }
}

void heap_free(heap *restrict h) {
    free(h->elements);
    free(h->pos);
    free(h->used);
}

void dump_heap(const heap *restrict h, const int *restrict map) {
    for(int i = 0; i < h->n; i++) {
        if(h->used[i]) {
            printf("%d %lf\n", map[h->elements[i].pos], h->elements[i].data);
        }
    }
}
