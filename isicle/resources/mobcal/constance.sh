#!/bin/bash

. /etc/bashrc
module purge
module load intel/18.0.0

PARAMS=$1
ATOMS=$2
FILE=$3

resources/mobcal/mobcal_constance "$PARAMS" "$ATOMS" "$FILE" "${FILE%.*}".out
