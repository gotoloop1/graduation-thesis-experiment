#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")

mkdir -p orders
mkdir -p working/reordering
mkdir -p results/reordering

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        name="working/reordering/${t}_${s}"
        mkdir -p "${name}"
        cp reordering_job_template.sh "${name}/reordering_job.sh"
        sed -i -e "s|matrix_file='matrix_file'|matrix_file='../../../data/${t}.mtx'|g" "${name}/reordering_job.sh"
        sed -i -e "s|order_file='order_file'|order_file='../../../orders/${t}_${s}'|g" "${name}/reordering_job.sh"
        sed -i -e "s|resulf_file='resulf_file'|resulf_file='../../../results/reordering/${t}_${s}.txt'|g" "${name}/reordering_job.sh"
        sed -i -e "s|seed='seed'|seed='${s}'|g" "${name}/reordering_job.sh"
    done
done
