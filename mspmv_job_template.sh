#!/bin/bash

#PJM -L rscgrp=regular
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=omp_thread_num
#PJM -L elapse=07:00:00
#PJM -g go34
#PJM -j

command_file='command_file'
matrix_file='matrix_file'
order_file='order_file'
resulf_file='resulf_file'
repeat='repeat'

echo -n "" > "${resulf_file}"

echo -n 'random ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_random.txt" "$((${repeat} / 16))" >> "${resulf_file}"

echo -n 'original ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_original.txt" "${repeat}" >> "${resulf_file}"

echo -n 'RCM ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_RCM.txt" "${repeat}" >> "${resulf_file}"

echo -n 'EBCM ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_EBCM.txt" "${repeat}" >> "${resulf_file}"

echo -n 'NBCM ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_NBCM.txt" "${repeat}" >> "${resulf_file}"

echo -n 'HCM ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_HCM.txt" "${repeat}" >> "${resulf_file}"

echo -n 'CRHCM ' >> "${resulf_file}"
"${command_file}" "${matrix_file}" "${order_file}_CRHCM.txt" "${repeat}" >> "${resulf_file}"
