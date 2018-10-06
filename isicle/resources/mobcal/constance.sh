#!/bin/bash

source /etc/bashrc
module purge
module load intel/18.0.0

PARAMS=$1
ATOMS=$2
FILE=$3

isicle/resources/mobcal/mobcal_constance "$PARAMS" "$ATOMS" "$FILE" "${FILE%.*}".out
