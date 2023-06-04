#!/bin/bash

#ZW: I got the following script from https://www.biostars.org/p/13452/ (written by Frederic Mahe). It cleanes up NAMES.dmp so that only scientific names are used. Other names and abbreviations shouldn't be used anymore.
WORKINGDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/"

NAMES="${WORKINGDIRECTORY}""names.dmp"

## Limit search space to scientific names
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"

