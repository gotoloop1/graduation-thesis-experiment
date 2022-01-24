#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
mspmv_numbers=("1" "2" "4" "8" "16")

mkdir -p working/mspmv_mnum
mkdir -p results/mspmv_mnum

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        for n in ${mspmv_numbers[@]}; do
            name="working/mspmv_mnum/${t}_${s}_${n}"
            mkdir -p "${name}"
            cp mspmv_job_template.sh "${name}/mspmv_job.sh"
            sed -i -e "s|#PJM --omp thread=omp_thread_num|#PJM --omp thread=56|g" "${name}/mspmv_job.sh"
            sed -i -e "s|command_file='command_file'|command_file='../../../mspmv_${n}'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|matrix_file='matrix_file'|matrix_file='../../../data/${t}.mtx'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|order_file='order_file'|order_file='../../../orders/${t}_${s}'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|resulf_file='resulf_file'|resulf_file='../../../results/mspmv_mnum/${t}_${s}_${n}.txt'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|repeat='repeat'|repeat='2048'|g" "${name}/mspmv_job.sh"
        done
    done
done
