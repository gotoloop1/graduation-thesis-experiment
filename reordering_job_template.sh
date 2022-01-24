#!/bin/bash

#PJM -L rscgrp=regular
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=1
#PJM -L elapse=04:00:00
#PJM -g go34
#PJM -j

matrix_file='matrix_file'
order_file='order_file'
resulf_file='resulf_file'
seed='seed'

echo -n "" > "${resulf_file}"

echo -n 'random ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_random.txt" 0 "${seed}" >> "${resulf_file}"

echo -n 'original ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_original.txt" 1 "${seed}" >> "${resulf_file}"

echo -n 'RCM ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_RCM.txt" 2 "${seed}" >> "${resulf_file}"

echo -n 'EBCM ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_EBCM.txt" 3 "${seed}" >> "${resulf_file}"

echo -n 'NBCM ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_NBCM.txt" 4 "${seed}" >> "${resulf_file}"

echo -n 'HCM ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_HCM.txt" 5 "${seed}" >> "${resulf_file}"

echo -n 'CRHCM ' >> "${resulf_file}"
../../../reordering "${matrix_file}" "${order_file}_CRHCM.txt" 6 "${seed}" >> "${resulf_file}"
