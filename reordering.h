#ifndef REORDERING_INCLUDED
#define REORDERING_INCLUDED

void reverse_cuthill_mckee(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
);

void edge_based_cost_minimization(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
);

void node_based_cost_minimization(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
);

void hybrid_cost_minimization(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
);

void coarsen_refine_hybrid_cost_minimization(
    const int n,
    const int *restrict edge_separator,
    const int *restrict edge,
    int *restrict order
);

#endif
