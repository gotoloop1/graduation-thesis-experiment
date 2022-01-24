#!/bin/bash -xe

testcases=("bcsstk36" "bmw7st_1" "ct20stif" "Fault_639" "Flan_1565" "Hook_1498" "ldoor" "msdoor" "nasasrb" "pwtk" "smt")
seeds=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")

for t in ${testcases[@]}; do
    for s in ${seeds[@]}; do
        name="working/reordering/${t}_${s}"
        cd "${name}"
        pjsub reordering_job.sh
        cd ../../../
    done
done
