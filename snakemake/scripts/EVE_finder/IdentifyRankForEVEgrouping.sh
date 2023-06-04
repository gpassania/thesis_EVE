#!/bin/bash

#ZW: I based the below script on code from https://www.biostars.org/p/13452/ (written by Frederic Mahe)
SPECIESNAME="${1}"
WORKINGDIR="${2}"
FILENAME="${3}"
NAMES="${WORKINGDIR}""names.dmp"
NODES="${WORKINGDIR}""nodes.dmp"
TAXONOMY=""
RANKCATEGORIES=""

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
   grep --max-count=1 "^${1}"$'\t' "${2}" | cut --fields="${3}"
}

#Obtain taxid from species name
get_name_or_taxidFromSpecies()
{
   grep --max-count=1 "${1}"$'\t' "${2}" | cut --fields="${3}"
}

# Get the taxid corresponding to the species name
TAXID=$(get_name_or_taxidFromSpecies "${SPECIESNAME}" "${NAMES}" "1")

#ZW:Obtain the rank at which current TAXID is at. Need to check for Flaviviridae family. Odd numbered collumns because all even columns are "|"
#RANK=$(get_name_or_taxid "${TAXID}" "${NODES}" "5")

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
while [[ "${TAXID}" -gt 1 ]] ; do
    # Obtain the scientific name corresponding to a taxid
    NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    # Obtain the parent taxa taxid
    PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    RANK=$(get_name_or_taxid "${TAXID}" "${NODES}" "5")
    # Build the taxonomy path
    TAXONOMY="${NAME}\t${TAXONOMY}"
    RANKCATEGORIES="${RANK}\t${RANKCATEGORIES}"
    TAXID="${PARENT}"
done
#Create two output files. Each line of each output file corresponds to a virus entry from initial bed file.
#TaxonomyHierarchyByEVE is a tab seperated output file of each taxonomy assignment for a given virus (according to NCBI records)
#RankingHierarchyByEVE is a tab seperated output file of the rank (ie Kingdom, order, family, genus...) corresponding to each entry in TaxonomyHierarchyByEVE
echo -e ${TAXONOMY} >> "${WORKINGDIR}"TaxonomyHierarchyByEVE_"$FILENAME".txt
echo -e ${RANKCATEGORIES} >> "${WORKINGDIR}"RankingHierarchyByEVE_"$FILENAME".txt
exit 0
