#!/usr/env/python3

"""
Author: Gian Passania

Script to create a summary of vyper result tables of all analyzed fastq files

usage:
python3 summartize_top.py [top10files]
"""

#import statements

from sys import argv
import os

#Functions

def append_lines(top10file):
	"""
	for each top10file, append lines to new summary file
	
	top10file: str, full path to top10file resulting from vyper
	"""
	outpath = top10file.split("/")[:-2]
	outname = "/".join(outpath) + "/top10summary.txt"
	if not os.path.exists(outname):
		x = open(outname, "w")
		x.write("WeightedHits\tUniqueHits\tVirus\tAccession\tSampleID\n")
		x.close()
	f = open(outname, "a")
	with open(top10file, "r") as reader:
		sampleID = top10file.split(".")[0]
		sampleID = sampleID.split("/")[-1] + "\n"
		for line in reader:
			if line.startswith("Weighted"):
				continue
			line = line.strip()
			parts = line.split()
			parts.append(sampleID)
			newline = "\t".join(parts)
			f.write(newline)
	f.close()


#main function
def main():
	for files in argv[1:]:
		append_lines(files)

if __name__ == "__main__":
	main()
