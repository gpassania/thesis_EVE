#!/usr/bin/env python

import os
import sys
import csv

queryFile = sys.argv[1]
outfile = queryFile.split(".")[:-1]
outfile = ".".join(outfile)



with open(queryFile) as qF:
	Host_chrom = [part[0] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	Host_start = [int(part[1]) for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	Host_end = [int(part[2]) for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	named = [part[3] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	score = [part[4] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]	
with open(queryFile) as qF:
	framed = [part[5] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	pident = [float(part[6]) for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	qcovs = [part[7] for part in [entry.split('\t') for entry in qF.readlines()[0:]]]	
with open(queryFile) as qF:
	Virus_start = [int(part[8]) for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	Virus_end = [int(part[9]) for part in [entry.split('\t') for entry in qF.readlines()[0:]]]
with open(queryFile) as qF:
	Virus_length = [int(part[10]) for part in [entry.split('\t') for entry in qF.readlines()[0:]]]

name = []
for names in named:
	names = names.strip("\n")
	names = names.lstrip(" ")
	names = names.replace(" ", "_")
	name.append(names)	

c = ','.join('%s_%s_%s_%s' % t for t in zip(name, Host_chrom, Host_start, Host_end))
list = c.split(',')

start_end = []
unsorted = zip(Host_start, Host_end)
strand = []

for se in unsorted:
	if se[0] > se[1]:
		strand.append('-')
	elif se[0] < se[1]:
		strand.append('+')
	se = sorted(se)
	start_end.append(se)
#start_end = sorted(start_end)
#print start_end

def unzip(iterable):
    return zip(*iterable)

if len(start_end):
	start, end = unzip(start_end)    

queries = zip(Host_chrom, start, end, list, score, strand, Virus_start, Virus_end, Virus_length, pident, qcovs)
print ("Putative EVEs for " + queryFile.split('_')[0] + " " + queryFile.split('_')[1] + ":"), len(queries)
#print queries

with open(outfile + '_EVEs.bed', "wb") as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerows(queries)
	

    
    
    
    
    
    
    
    
    
    
