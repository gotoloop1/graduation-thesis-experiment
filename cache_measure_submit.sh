#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
programs=("not_mspmv" "mspmv")

submit_begin="0"
submit_end="220"
count="0"

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        for p in ${programs[@]}; do
            if [ "${count}" -ge "${submit_begin}" -a "${count}" -lt "${submit_end}" ]; then
                name="working/cache_measure_thread/${t}_${s}_${p}"
                cd "${name}"
                pjsub cache_measure_job.sh
                cd ../../../
            fi
            count="$((${count} + 1))"
        done
    done
done
