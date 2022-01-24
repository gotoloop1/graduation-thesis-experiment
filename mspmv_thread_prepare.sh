#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
thread_numbers=("1" "2" "4" "8" "14" "28" "56")

mkdir -p working/mspmv_thread
mkdir -p results/mspmv_thread

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        for n in ${thread_numbers[@]}; do
            case "${n}" in
                "56" | "28" ) repeat="2048" ;;
                "14" | "8" ) repeat="1024" ;;
                "4" | "2" ) repeat="512" ;;
                "1" ) repeat="256" ;;
            esac
            name="working/mspmv_thread/${t}_${s}_${n}"
            mkdir -p "${name}"
            cp mspmv_job_template.sh "${name}/mspmv_job.sh"
            sed -i -e "s|#PJM --omp thread=omp_thread_num|#PJM --omp thread=${n}|g" "${name}/mspmv_job.sh"
            sed -i -e "s|command_file='command_file'|command_file='../../../mspmv'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|matrix_file='matrix_file'|matrix_file='../../../data/${t}.mtx'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|order_file='order_file'|order_file='../../../orders/${t}_${s}'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|resulf_file='resulf_file'|resulf_file='../../../results/mspmv_thread/${t}_${s}_${n}.txt'|g" "${name}/mspmv_job.sh"
            sed -i -e "s|repeat='repeat'|repeat='${repeat}'|g" "${name}/mspmv_job.sh"
        done
    done
done
