#!/bin/bash

source /etc/bashrc
module purge
module load intel/ips_18

PARAMS=$1
ATOMS=$2
FILE=$3

isicle/resources/mobcal/mobcal_cascade "$PARAMS" "$ATOMS" "$FILE" "${FILE%.*}".out