#!/bin/bash

WORKINGDIR="${1}"
FILENAME="${2}"

#Extract virus name only, export the names to virusNames.txt, then replace _ with spaces

#Use this line if only one virus name is contained between '[...]' per line in the bed file (as in Aag_Contigs_EVEs_NEW.bed)
cat "${WORKINGDIR}""$FILENAME" | grep -oP '\[\K[^]]+' > "${WORKINGDIR}"virusNamesByEVEentry_"$FILENAME".txt

#Use this line if there are multiple instances of '[virusName]' per line. This only looks at the first column to get the appropriate name
#cat "${WORKINGDIR}""$FILENAME" | cut -f 1 | grep -oP '\[\K[^]]+' > "${WORKINGDIR}"virusNames.txt

#Remove the underscores in the virus names
sed -i 's/_/ /g' "${WORKINGDIR}"virusNamesByEVEentry_"$FILENAME".txt
