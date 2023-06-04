#!/bin/bash

#Directory with extractVirusNamesForEVEgrouping.sh, IdentifyRankForEVEgrouping.sh, and ClassifyEVEtaxonomy.py
SCRIPTDIRECTORY="/path/to/directory/containing/necessary/scripts/"

#Directory containing bed file to be analyzed. Include ending "/"
#Should also contain names.dmp, nodes.dmp
WORKINGDIRECTORY="/path/to/directory/containing/bedFileWithColumnContainigVirusNamesToBeClassified"

#Name of above bed file to analyze
FILETOVISUALIZE="closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt"

#For how I set this up (probablly needs to be fixed later), each line of input NEEDS to have a corresponding line in output (virusNamesByEVEentry.txt). Need to ensure number of lines is the SAME!!!
${SCRIPTDIRECTORY}/extractVirusNamesForEVEgrouping.sh "${WORKINGDIRECTORY}" "${FILETOVISUALIZE}"

#Output taxonomy hierarchy for each virus in input
cat "${WORKINGDIRECTORY}"virusNamesByEVEentry_"${FILETOVISUALIZE}".txt | while read SPECIESNAME ; do bash ${SCRIPTDIRECTORY}/IdentifyRankForEVEgrouping.sh "$SPECIESNAME" "${WORKINGDIRECTORY}" "${FILETOVISUALIZE}"; done

wc -l "${WORKINGDIRECTORY}""${FILETOVISUALIZE}"
wc -l "${WORKINGDIRECTORY}"virusNamesByEVEentry_"${FILETOVISUALIZE}".txt
wc -l "${WORKINGDIRECTORY}"TaxonomyHierarchyByEVE_"${FILETOVISUALIZE}".txt
wc -l "${WORKINGDIRECTORY}"RankingHierarchyByEVE_"${FILETOVISUALIZE}".txt

#Create and output txt file with added EVE taxonomy info
python ${SCRIPTDIRECTORY}/ClassifyEVEtaxonomy.py "${WORKINGDIRECTORY}" "${FILETOVISUALIZE}"
