#!/usr/env/python3

"""
Author: Gian Passania

Script to combine results from extracts files and add assembly levels

usage: 
python3 combine_results_detect.py [extractfile1] [extractfile2] ... [last_extractfile]

extractfile: full path to extracts file resulting from detect_final_filter.py
"""

#import statements

from sys import argv
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "gian.passania@wur.nl"
import time
import subprocess as sp


#Functions
	
def append_lines(extractfile):
	"""
	for each extractfile, get assembly level and write to file
	
	extractsfile: str, full path to extractfile
	"""
	#set variables
	outpath = extractfile.split("/")[:-3]
	outname = "/".join(outpath) + "/combined_extracts.txt"
	ref = extractfile.split("/")[-3]
	#if file does not exist: create file and write headers
	if not os.path.exists(outname):
		x = open(outname, "w")
		x.write("Scaffold\tStart\tEnd\tLength\tID\tFraction\tVirusName\tvirus_eval\tPercentage_ID\tSubject_eval\tSampleID\tAssemblyLevel\n")
		x.close()
	#run ncbi datasets to get assembly level
	acc = ref.split("_")[0:2]
	ID = "_".join(acc)
	cmd = "datasets summary genome accession {} > tempfile".format(ID)
	sp.check_call(cmd, shell=True)
	with open(os.getcwd()+"/tempfile", "r") as reader:
		line = reader.readlines()
		line = line[0]
		idx = line.find("assembly_level")
		idx += 17
		level = line[idx:idx+10]
		if "Scaffold" in level:
			level = level[:-2]
		if "Contig" in level:
			level = level[:-4]
	f = open(outname, "a")
	#append lines to combined extracts file
	with open(extractfile, "r") as reader:
		for line in reader:
			if line.startswith("Scaffold"):
				continue
			line = line.strip()
			parts = line.split("\t")
			parts.append(ref)
			parts.append(level+"\n")
			newline = "\t".join(parts)
			f.write(newline)
	f.close()



#main function
def main():
	for files in argv[1:]:
		append_lines(files)

if __name__ == "__main__":
	main()
