#!/bin/sh

JOBS=(
    minimal-pipeline.bash
)

BASEDIR=$(dirname "$0")

echo ==========================
echo $i
echo ==========================
bash -e build.bash 

# All other jobs
for i in "${JOBS[@]}"
do

    echo ==========================
    echo $i
    echo ==========================

    bash -e $BASEDIR/jobs/$i
    exitcode=$?

    if [ $exitcode -ne 0 ]; then
        echo ***JOB FAILED*** 
        echo $i
        exit 1
    else
        echo --- $i PASSED ---
    fi
done

echo ==========================
echo ALL TEST PASSED
echo ==========================
