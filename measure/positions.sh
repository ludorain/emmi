#! /usr/bin/env bash

if [ $# -ne 5 ]; then
    echo " usage: $0 [database] [board] [sensor] [tag] [axis] "
    exit 1
fi
DATABASE=$1
BOARD=$2
SENSOR=$3
TAG=$4
AXIS=$5

while read -r board sensor tag axis start step count; do
    [[ ${board} != ${BOARD} || ${sensor} != ${SENSOR} || ${tag} != ${TAG} || ${axis} != ${AXIS} ]] && continue
    end=$(( start + (count - 1) * step ))
    [[ $step = 0 ]] && values=$start || values=$(seq $start $step $end | xargs)
    echo ${values}
done < ${DATABASE}
