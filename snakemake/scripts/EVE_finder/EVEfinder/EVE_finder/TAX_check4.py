#!/usr/bin/env python

######			######			######			######			######			######			######

#####		Imports		#####
import csv
import sys

######			######			######			######			######			######			######

#####		Globals		#####

### Define TAXIDs for all viruses, Blast_NR viruses and TAXIDs, and Fasta viruses	###

Virus_TAXID_List =  "RNA_and_ssDNA_viral_TAXIDs"
Blast_NR =  sys.argv[1]
Fasta_file = sys.argv[2]
Bedfile = sys.argv[3]

######			######			######			######			######			######			######

#####		Functions		#####

### Open Fasta file and define the Header and the Sequence. Creates a dictionary	###

with open(Fasta_file) as file:
    names=[part[0] for part in [entry.partition('\n') for entry in file.read().split('>')[1:]]]
with open(Fasta_file) as file:
    seqs=[part[2].replace("\n","") for part in [entry.partition('\n') for entry in file.read().split('>')[1:]]]

Potentials = dict(zip(names,seqs)) #return dictionary of name and sequence

### Open Bed file and define the Name and the entire entry. Creates a dictionary	###

with open(Bedfile) as bF:
    bed_entries = bF.readlines()[0:]
with open(Bedfile) as bF:
    bed_names=[part[3] for part in [entry.split('\t') for entry in bF.readlines()[0:]]]

bed_potentials = dict(zip(bed_names,bed_entries)) #return dictionary of name and entry

### Open Viral_TAXID file and define the Taxonomy ID and the Virus name	###

Virus_list = []

with open(Virus_TAXID_List) as qF:
	TAXIDs = [part[6] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]

for IDs in TAXIDs:
	IDs = IDs.rstrip()
	#print IDs
	for ID in IDs.split(" "):
		#print ID
		Virus_list.append(ID)

### Open Blast_NR file and define the Taxonomy ID and the Virus name. Creates a dictionary	###

with open(Blast_NR) as Ef:
	entries = Ef.readlines()[0:]
	
IDs=[]	
	
for entry in entries:
#	print entry.split("\t")[0]
	if entry.split("\t")[0] not in IDs:
		IDs.append(entry.split("\t")[0])

### Define Scores. Score for each potential EVE is the number of Blast_NR results that are viral over the total	###	

Confidence_scores=[]

for ID in IDs:
	Plus = 0
	Minus = 0
	for entry in entries:
#		print ID.split("|")[3]
#		print entry.split("\t")[0].split("|")[3]
		if ID.split("|")[3] == entry.split("\t")[0].split("|")[3]:
#			print entry.split("\t")[11].rstrip()
			if entry.split("\t")[11].rstrip() in Virus_list:
				Plus = Plus + 1
			else:
				Minus = Minus + 1
	Confidence_scores.append(float(Plus * 100 / (Plus + Minus)))
	
Confidence_dict = dict(zip(IDs, Confidence_scores))

### Define top result from Blast_NR file	###	

Scores=[]
Sorted_Scores=[]
Threshold=[]
Final=[]
for ID in IDs:
	for entry in entries:
		if entry.split("\t")[0] in ID:
#			print entry.split("\t")[0], float(entry.split("\t")[4])
			Scores.append(float(entry.split("\t")[4]))
			Sorted_Scores.append(float(entry.split("\t")[4]))
			Threshold.append(entry)
	Score_dict=dict(zip(Scores, Threshold))
	Sorted_Scores = sorted(Sorted_Scores)
#	print Sorted_Scores
#	print Score_dict[Sorted_Scores[0]]
	Final.append(Score_dict[Sorted_Scores[0]])
	Scores=[]
	Sorted_Scores=[]
	Threshold=[]	
		

### Define Hits. Hits are when the top result from the Blast_NR file is viral. Unknown is when there is no Blast_NR hit	###	

Found=[]
Hits = []
Misses = []

for unknown in names:
	for f in Final:
#		print f.split("\t")[0]
		if unknown == f.split("\t")[0]:
			Found.append(unknown)

for pair in Final:
#	print int(pair.split("\t")[11].rstrip())
	if pair.split("\t")[11].rstrip() in Virus_list:
		if pair not in Hits:
			Hits.append(pair)
	
for pair in Final:
	if pair not in Hits:
		if pair not in Misses:
			Misses.append(pair)
				
				
######			######			######			######			######			######			######
#####		Bedfile Output		#####

Bed_Hits=[]
Bed_Misses=[]
Bed_UNKNOWN=[]

for Hit in Hits:
#	print Hit.split("\t")[0]
	if Hit.split("\t")[0] in bed_potentials:
		Bed_Hits.append((bed_potentials[Hit.split("\t")[0]].rstrip() + "\t" + str(Confidence_dict[Hit.split("\t")[0]]) + "\n"))

for Miss in Misses:
	if Miss.split("\t")[0] in bed_potentials:
		Bed_Misses.append((bed_potentials[Miss.split("\t")[0]].rstrip() + "\t" + str(Confidence_dict[Miss.split("\t")[0]]) + "\n"))

Bed_UNKNOWN=[]
for unknown in names:
	if unknown not in Found:
		Bed_UNKNOWN.append(bed_potentials[unknown].rstrip() + "\tN/A\n")

with open(Bedfile.split(".")[0] + '_Hits.bed', "wb") as f:
	for Test in Bed_Hits:
		f.write(Test)

with open(Bedfile.split(".")[0] + '_Misses.bed', "wb") as f:
	for Test in Bed_Misses:
		f.write(Test)

with open(Bedfile.split(".")[0] + '_UNKNOWN.bed', "wb") as f:
	for Test in Bed_UNKNOWN:
		f.write(Test)

#####		Fasta Output		#####

Hit_output = []
Miss_output = []

for Hit in Hits:
	if Hit.split("\t")[0] in Potentials:
		Hit_output.append(">" + Hit.split("\t")[0])	
		Hit_output.append(Potentials[Hit.split("\t")[0]])

for Miss in Misses:
	if Miss.split("\t")[0] in Potentials:
		Miss_output.append(">" + Miss.split("\t")[0])	
		Miss_output.append(Potentials[Miss.split("\t")[0]])
		
UNKNOWN=[]
for unknown in names:
	if unknown not in Found:
			UNKNOWN.append(">" + unknown)
			UNKNOWN.append(Potentials[unknown])
	
with open(Fasta_file.split(".")[0] + '_Hits.fasta', "wb") as f:
	for Test in Hit_output:
		f.write(Test + "\n")

with open(Fasta_file.split(".")[0] + '_Misses.fasta', "wb") as f:
	for Test in Miss_output:
		f.write(Test + "\n")

with open(Fasta_file.split(".")[0] + '_UNKNOWN.fasta', "wb") as f:
	for Test in UNKNOWN:
		f.write(Test + "\n")

#####		Output Check		#####

print str(Fasta_file.split("_")[0]) + ' ' + str(Fasta_file.split("_")[1]) + "\t" + str(len(Hits)) + " HITS" 
print str(Fasta_file.split("_")[0]) + ' ' + str(Fasta_file.split("_")[1]) + "\t" + str(len(Misses)) + " MISSES"
print str(Fasta_file.split("_")[0]) + ' ' + str(Fasta_file.split("_")[1]) + "\t" + str(len(UNKNOWN)/2) + " UNKNOWN"  
print str(Fasta_file.split("_")[0]) + ' ' + str(Fasta_file.split("_")[1]) + "\t" + str(len(names)) + " TOTAL" 
if len(Hits) + len(Misses) + len(UNKNOWN)/2 == len(names) == len(Bed_UNKNOWN)+len(Bed_Misses)+len(Bed_Hits):
	print "CHECK"
else:
	print "ERROR"

######			######			######			######			######			######			######
		
