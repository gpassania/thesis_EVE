#!/usr/bin/env python

## Imports ##
import csv
import sys

Virus_TAXID_List =  "RNA_and_ssDNA_viral_TAXIDs"


NR_file = sys.argv[1]
genus = sys.argv[1].split("_")[0].lower()

### Define all of the TAXIDs that are viral ###


with open(Virus_TAXID_List) as qF:
	VIRAL_TAXIDs = [space.split(" ")[0] for space in [part[6] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]]

	
### Find the Virus name from the Bedfile ###

Hits=[]

with open(NR_file) as bF:
	entries = bF.readlines()[0:]
	
for entry in entries:
	tax = entry.split('\t')[11]
	
	if tax.rstrip() in VIRAL_TAXIDs:	
		Hits.append(entry)


for entry in entries:	
	name = entry.split('\t')[3]
	if 'hypothetical' not in name.lower() and 'predicted' not in name.lower() and genus not in name.lower():
		Hits.append(entry)
		
with open(NR_file + "_filtered", "wb") as f:
	for Test in Hits:
		f.write(Test)