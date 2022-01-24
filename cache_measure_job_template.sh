#!/bin/bash

#PJM -L rscgrp=regular
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=56
#PJM -L elapse=07:00:00
#PJM -g go34
#PJM -j

events='mem_load_retired.fb_hit,mem_load_retired.l1_hit,mem_load_retired.l1_miss,mem_load_retired.l2_hit,mem_load_retired.l2_miss,mem_load_retired.l3_hit,mem_load_retired.l3_miss'

command_file='command_file'
matrix_file='matrix_file'
order_file='order_file'
resulf_file='resulf_file'
repeat='2048'

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_random.txt" "${repeat}" 2> "${resulf_file}_random.txt"

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_original.txt" "${repeat}" 2> "${resulf_file}_original.txt"

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_RCM.txt" "${repeat}" 2> "${resulf_file}_RCM.txt"

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_EBCM.txt" "${repeat}" 2> "${resulf_file}_EBCM.txt"

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_NBCM.txt" "${repeat}" 2> "${resulf_file}_NBCM.txt"

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_HCM.txt" "${repeat}" 2> "${resulf_file}_HCM.txt"

perf stat -e "${events}" "${command_file}" "${matrix_file}" "${order_file}_CRHCM.txt" "${repeat}" 2> "${resulf_file}_CRHCM.txt"
