#!/bin/sh
#
#================================================================
# run_tests.sh
#================================================================
#
#% DESCRIPTION
#%    Runs integration tests
#%    Each test should be defined as a bash script and appended in the JOBS variable below

JOBS=(
    minimal-pipeline.bash
)

IMAGE=bov-tb:test
BASEDIR=$(dirname "$0")

# Build docker image
echo ==========================
echo Build
echo ==========================

bash -e -c "docker build -t $IMAGE $BASEDIR/../"

# Integration Tests
for i in "${JOBS[@]}"
do
    echo ==========================
    echo $i
    echo ==========================

    # Run job in docker
    sudo docker run -v $BASEDIR/jobs:/jobs/ --rm $IMAGE /bin/bash/ -e /jobs/$i
    exitcode=$?

    # Handle exit code
    if [ $exitcode -ne 0 ]; then
        echo *** JOB FAILED: $i ***
        exit 1
    else
        echo --- JOB PASSED: $i ---
    fi
done

echo ==========================
echo ALL TEST PASSED
echo ==========================
