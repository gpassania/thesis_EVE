#!/usr/bin/python3
"""
Author: Gian Passania

Script to summarize overlapping regions to one liners and add taxonomy names

usage:

python3 sum_overlaps.py [overlapping_regions.txt]
"""
#import statements
from sys import argv
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "gian.passania@wur.nl"
import time
#functions

def get_oneliners(regionsfile):
	"""
	for each overlapping region, calc length and write oneliner for each region
	
	regionsfile: str, full path to overlapping regionsfile
	"""
	#set lastline
	with open(regionsfile, "r") as reader:
		lastline = reader.readlines()[-1]
	with open(regionsfile, "r") as reader:
		newline_list = []
		init = 1
		starters = []
		enders = []
		for line in reader:
			if line.startswith("Chr"):
				continue
			line = line.strip()
			parts = line.split("\t")
			clusternum = int(parts[-1])
			if clusternum != init:
				start = min(starters)
				end = max(enders)
				oldline[1] = str(start)
				oldline[2] = str(end)
				length = end - start
				oldline.append(str(length))
				newline = "\t".join(oldline)
				newline_list.append(newline)
				starters = []
				enders = []
				starters.append(int(parts[1]))
				enders.append(int(parts[2]))
				init +=1
			else:
				starters.append(int(parts[1]))
				enders.append(int(parts[2]))
				oldline = parts
				
			if line in lastline:
				start = min(starters)
				end = max(enders)
				oldline[1] = str(start)
				oldline[2] = str(end)
				length = end - start
				oldline.append(str(length))
				newline = "\t".join(oldline)

	
	return newline_list		

def entrez(newline_list):
	"""
	for every oneliner, get full organism name and append to line
	
	newline_list: list, contains oneliners for each overlapping region
	"""
	final = []
	for line in newline_list:
		parts = line.split("\t")
		acc = parts[5]
		handle = Entrez.efetch(db="nucleotide", id = acc, rettype= "gb", retmode="text")
		x = SeqIO.read(handle, "genbank")
		print(x.annotations)
		parts.append(x.annotations["organism"])
		finalline = "\t".join(parts)
		final.append(finalline)
		time.sleep(2)
	return final
		
def write_output(final, outname):
	"""
	write output to overlapping_summary.txt
	
	final: list, contains oneliners for each overlapping region  + organisms names
	outname: str, full path for directory where output is placed
	"""
	f = open(outname, "w")
	f.write("Chr\tCluster_start\tCluster_end\thits1\tVirus1\tID1\thits2\tVirus2\tID2\thits3\tVirus3\tID3\tSampleID\tclusternumber\tLength\tfullVirusName\n")
	for line in final:
		f.write(line)
		f.write("\n")
	f.close()
	
def main():
	"""
	main function of this script
	"""
	regionsfile = argv[1]
	outname = regionsfile.split("/")[:-1]
	outname = "/".join(outname)
	outname += "/overlapping_summary.txt"
	newline_list = get_oneliners(regionsfile)
	final = entrez(newline_list)
	write_output(final, outname)

	
if __name__ == "__main__":
	main()
