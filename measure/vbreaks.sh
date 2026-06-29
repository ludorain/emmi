#! /usr/bin/env bash

if [ $# -ne 4 ]; then
    echo " usage: $0 [database] [board] [sensor] [what] "
    exit 1
fi
DATABASE=$1
BOARD=$2
SENSOR=$3
WHAT=$4

while read -r board sensor at20c coeff; do
    [[ ${board} != ${BOARD} || ${sensor} != ${SENSOR} ]] && continue
    [[ ${WHAT} == "at20c" ]] && echo ${at20c}
    [[ ${WHAT} == "coeff" ]] && echo ${coeff}
done < ${DATABASE}
