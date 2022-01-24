#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
thread_numbers=("1" "2" "4" "8" "14" "28" "56")

submit_begin="0"
submit_end="770"
count="0"

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        for n in ${thread_numbers[@]}; do
            if [ "${count}" -ge "${submit_begin}" -a "${count}" -lt "${submit_end}" ]; then
                name="working/mspmv_thread/${t}_${s}_${n}"
                cd "${name}"
                pjsub mspmv_job.sh
                cd ../../../
            fi
            count="$((${count} + 1))"
        done
    done
done
