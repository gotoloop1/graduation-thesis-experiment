#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
programs=("not_mspmv" "mspmv")

mkdir -p working/cache_measure_thread
mkdir -p results/cache_measure_thread

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        for p in ${programs[@]}; do
            name="working/cache_measure_thread/${t}_${s}_${p}"
            mkdir -p "${name}"
            cp cache_measure_job_template.sh "${name}/cache_measure_job.sh"
            sed -i -e "s|command_file='command_file'|command_file='../../../${p}'|g" "${name}/cache_measure_job.sh"
            sed -i -e "s|matrix_file='matrix_file'|matrix_file='../../../data/${t}.mtx'|g" "${name}/cache_measure_job.sh"
            sed -i -e "s|order_file='order_file'|order_file='../../../orders/${t}_${s}'|g" "${name}/cache_measure_job.sh"
            sed -i -e "s|resulf_file='resulf_file'|resulf_file='../../../results/cache_measure_thread/${t}_${s}_${p}'|g" "${name}/cache_measure_job.sh"
        done
    done
done
