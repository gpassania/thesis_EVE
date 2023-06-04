#!/usr/bin/env python

## Imports ##
import csv
import sys
## DEFINITIONS ##

infile = sys.argv[1]
outfile = sys.argv[2]

a = open(infile)
csvfile = open(outfile, 'wb')
spamwriter = csv.writer(csvfile, delimiter='\t')
Host_names = []
evalues = []
sorted_evalues = []
strands = []
Virus_starts = []
Virus_ends = []
Virus_lengths = []
pidents = []
qcovs = []
	
EVEs = 0

for line in a:
	line = line.strip()
	line = line.split("\t")
	
	
###		Definitions		###	
	
	Host_name = line [3]
	evalue = line[4]
	strand = line[5]
	Virus_start = line[6]
	Virus_end = line[7]
	Virus_length = line[8]
	pident = line[9]
	qcov = line[10]
	
###				###

	if "," in Host_name:
		Host_name = Host_name.split(",")
		for name in Host_name:
			Host_names.append(name)
	else:
		Host_names.append(Host_name)

	if "," in evalue:
		evalue = evalue.split(",")
		for value in evalue:
			evalues.append(str(float(value)))
			sorted_evalues.append(str(float(value)))
	else:
		evalues.append(str(float(evalue)))
		sorted_evalues.append(str(float(evalue)))

	if "," in strand:
		strand = strand.split(",")[0]


	if "," in Virus_start:
		Virus_start = Virus_start.split(",")
		for start in Virus_start:
			Virus_starts.append(start)
	else:
		Virus_starts.append(Virus_start)

	if "," in Virus_end:
		Virus_end = Virus_end.split(",")
		for end in Virus_end:
			Virus_ends.append(end)
	else:
		Virus_ends.append(Virus_end)

	if "," in Virus_length:
		Virus_length = Virus_length.split(",")
		for length in Virus_length:
			Virus_lengths.append(length)
	else:
		Virus_lengths.append(Virus_length)

	if "," in pident:
		pident = pident.split(",")
		for id in pident:
			pidents.append(id)
	else:
		pidents.append(pident)
		
	if "," in qcov:
		qcov = qcov.split(",")
		for cov in qcov:
			qcovs.append(cov)
	else:
		qcovs.append(qcov)		
	
###				###

	sorted_evalues = sorted(sorted_evalues)		
	
	Host_names_dict = dict(zip(evalues, Host_names))
	Virus_starts_dict = dict(zip(evalues, Virus_starts))
	Virus_ends_dict = dict(zip(evalues, Virus_ends))
	Virus_lengths_dict = dict(zip(evalues, Virus_lengths))
	pidents_dict = dict(zip(evalues, pidents))
	qcovs_dict = dict(zip(evalues, qcovs))

###				###
	line[3] = Host_names_dict[str(sorted_evalues[0])]
	line[4] = sorted_evalues[0]
	line[5] = strand
	line[6] = Virus_starts_dict[str(sorted_evalues[0])]
	line[7] = Virus_ends_dict[str(sorted_evalues[0])]
	line[8] = Virus_lengths_dict[str(sorted_evalues[0])]
	line[9] = pidents_dict[str(sorted_evalues[0])]
	line[10] = int(float(qcovs_dict[str(sorted_evalues[0])]))
	
	Host_names = []
	evalues = []
	sorted_evalues = []
	strands = []
	Virus_starts = []
	Virus_ends = []
	Virus_lengths = []
	pidents = []
	qcovs = []

	spamwriter.writerow(line)
	EVEs = EVEs + 1


