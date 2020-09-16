#!/bin/sh

READS=$1
RESULTS=$2

sudo docker run --rm -it -v $READS:/reads/ -v $RESULTS:/results/ bov-tb