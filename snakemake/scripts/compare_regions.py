#!/usr/bin/python3
"""
Author: Gian Passania

Script to compare regions found across fastq sets

usage:

python3 compare_regions.py [regionfiles]

regionfiles: full paths to regionfiles 
"""
#import statements
from sys import argv
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "gian.passania@wur.nl"
import time
#functions

def add_fileID(regionfile):
	"""
	add fileID to each line and write to combined file
	
	regionfile: str, full path to regionfile
	"""
	outpath = regionfile.split("/")[:-2]
	outname = "/".join(outpath) + "/combined_clusters.txt"
	fileID = regionfile.split("/")[-1]
	fileID = fileID.split(".")[0]
	if not os.path.exists(outname):
		x = open(outname, "w")
		x.write("Chr\tCluster_start\tCluster_end\thits1\tVirus1\tID1\thits2\tVirus2\tID2\thits3\tVirus3\tID3\n")
		x.close()
	f = open(outname, "a")
	with open(regionfile, "r") as reader:
		for line in reader:
			if line.startswith("Chr"):
				continue
			line = line.strip()
			parts = line.split("\t")
			parts.append(fileID)
			newline = "\t".join(parts)
			f.write(newline)
			f.write("\n")
		f.close()
	return outname
	
def get_line_list(outname):
	"""
	get a list of all lines in combined file
	
	outname: str, full path to combined file
	"""
	line_list = []
	with open(outname, "r") as reader:
		for line in reader:
			if line.startswith("Chr"):
				continue
			line_list.append(line)
	return line_list
	
def get_same_regions(outname, line_list):
	"""
	get similar regions across all samples into a list of lists
	
	outname: str, full path to created combined file
	line_list: str, list of all lines in combined file
	"""
	same_regions = {}
	same_chroms = {}
	used_lines = []
	with open(outname, "r") as reader:
		for line in reader:
			if line in used_lines:
				continue
			for line2 in line_list:
				if line2 in used_lines:
					continue
				parts1 = line.split("\t")
				parts2 = line2.split("\t")
				if parts1[-1] == parts2[-1]:
					continue
				if parts1[4] != parts2[4]:
					continue
				if parts1[1] >= parts2[1] and parts1[1] <= parts2[2]:
					if line not in same_regions.keys():
						same_regions[line] = []
					same_regions[line].append(line2)
					used_lines.append(line)
					used_lines.append(line2)
				elif parts1[1] <= parts2[1] and parts1[1] >= parts2[2]:
					if line not in same_regions.keys():
						same_regions[line] = []
					same_regions[line].append(line2)
					used_lines.append(line)
					used_lines.append(line2)
				else:
					continue
					
		for key in same_regions.keys():
			chrom = key.split("\t")[0]
			for line in same_regions[key]:
				chrom2 = line.split("\t")[0]
				if chrom == chrom2:
					if key not in same_chroms.keys():
						same_chroms[key] = []
					same_chroms[key].append(line)
	return same_chroms


	
def write_output1(same_chroms, outname):
	"""
	write overlapping regions to file
	
	same_chroms: dict, where key = line and values = list of lines in same regions
	outname: str, full path combined file
	"""
	outfile = "/".join(outname.split("/")[:-1]) + "/overlapping_regions.txt"
	f = open(outfile, "w")
	f.write("Chr\tCluster_start\tCluster_end\thits1\tVirus1\tID1\thits2\tVirus2\tID2\thits3\tVirus3\tID3\tSampleID\tclusternumber\n")
	count = 1
	for key in same_chroms.keys():
		keys = key.strip()
		keyparts = keys.split("\t")
		keyparts.append(str(count)+"\n")
		f.write("\t".join(keyparts))
		for line in same_chroms[key]:
			line = line.strip()
			parts = line.split("\t")
			parts.append(str(count)+"\n")
			f.write("\t".join(parts))
		count +=1
	f.close()

def write_output2(same_chroms, line_list,  outname):
	"""
	write unique regions to file
	
	same_chroms: dict, where key = line and values = list of lines in same regions
	line_list: list, list of all lines in combined file
	outname: str, full path combined file
	"""
	outfile = "/".join(outname.split("/")[:-1]) + "/unique_regions.txt"
	f = open(outfile, "w")
	f.write("Chr\tCluster_start\tCluster_end\thits1\tVirus1\tID1\thits2\tVirus2\tID2\thits3\tVirus3\tID3\tSampleID\tFullVirusName\n")
	for lines in same_chroms.values():
		for line in lines:
			if line in line_list:
				line_list.remove(line)
	for line in line_list:
		line = line.strip()
		parts = line.split("\t")
		acc = parts[5]
		handle = Entrez.efetch(db="nucleotide", id = acc, rettype= "gb", retmode="text")
		x = SeqIO.read(handle, "genbank")
		parts.append(x.annotations["organism"])
		finalline = "\t".join(parts)
		time.sleep(1)
		f.write(finalline)
		f.write("\n")
	f.close()


def main():
	"""
	main function of this script
	"""
	for regionfile in argv[1:]:
		outname = add_fileID(regionfile)
	line_list = get_line_list(outname)
	same_chroms  = get_same_regions(outname, line_list)
	write_output1(same_chroms, outname)
	write_output2(same_chroms, line_list, outname)
if __name__ == "__main__":
	main()
