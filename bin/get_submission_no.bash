#!/bin/bash

# Error handling
set -eo pipefail

#=======
# Inputs
pair_id=$1

set +e; submission_number=$(echo $pair_id | grep -Eo '[0-9]{2,2}\-[0-9]{4,5}\-[0-9]{2,2}'); set -e
if [ -n "$submission_number" ]; then
    echo $submission_number
else
    echo $pair_id
fi