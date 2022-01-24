CC = icc

ifeq ($(CC),icc)
C_COMPILER_OPTION := -std=c99 -lm -lpthread -fopenmp -O3 -qopt-subscript-in-range -xHost
else
C_COMPILER_OPTION := -std=c99 -lm -lpthread -fopenmp -O3
endif

reordering: Makefile lib.h reordering.h matrix.h single_sort.c first_element_sort.c heap.c range_sum.c memory.c random.c matrix.c rcm.c ebcm.c nbcm.c hcm.c crhcm.c reordering.c
	$(CC) single_sort.c first_element_sort.c heap.c range_sum.c memory.c random.c matrix.c rcm.c ebcm.c nbcm.c hcm.c crhcm.c reordering.c $(C_COMPILER_OPTION) -o reordering

mspmv: Makefile lib.h matrix.h first_element_sort.c matrix.c mspmv.c
	$(CC) first_element_sort.c matrix.c mspmv.c $(C_COMPILER_OPTION) -o mspmv

mspmv_16: Makefile lib.h matrix.h first_element_sort.c matrix.c mspmv.c
	$(CC) first_element_sort.c matrix.c mspmv.c -DMSPMV_NUM=16 $(C_COMPILER_OPTION) -o mspmv_16

mspmv_8: Makefile lib.h matrix.h first_element_sort.c matrix.c mspmv.c
	$(CC) first_element_sort.c matrix.c mspmv.c -DMSPMV_NUM=8 $(C_COMPILER_OPTION) -o mspmv_8

mspmv_4: Makefile lib.h matrix.h first_element_sort.c matrix.c mspmv.c
	$(CC) first_element_sort.c matrix.c mspmv.c -DMSPMV_NUM=4 $(C_COMPILER_OPTION) -o mspmv_4

mspmv_2: Makefile lib.h matrix.h first_element_sort.c matrix.c mspmv.c
	$(CC) first_element_sort.c matrix.c mspmv.c -DMSPMV_NUM=2 $(C_COMPILER_OPTION) -o mspmv_2

mspmv_1: Makefile lib.h matrix.h first_element_sort.c matrix.c mspmv.c
	$(CC) first_element_sort.c matrix.c mspmv.c -DMSPMV_NUM=1 $(C_COMPILER_OPTION) -o mspmv_1

not_mspmv: Makefile lib.h matrix.h first_element_sort.c matrix.c not_mspmv.c
	$(CC) first_element_sort.c matrix.c not_mspmv.c $(C_COMPILER_OPTION) -o not_mspmv
