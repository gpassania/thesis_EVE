#!/usr/bin/env/python3
"""
split reference genome into seperate file for each chromosome
"""

#import statements
from sys import argv

#functions
def get_chrom(reference):
	with open(reference) as reader:
		chrom_dict = {}
		for line in reader:
			line = line.strip()
			if line.startswith(">") and "contig:" in line:
				header = line
				header = header.split()[-5]
				header = header[:-1]
				print(header)
				chrom_dict[header] = []
			elif line.startswith(">") and "chromosome:" in line:
				header = line
				header = header.split()[-1]
				header = "chr" + header
				print(header)
				chrom_dict[header] = []
			else:
				chrom_dict[header].append(line)
		for key in chrom_dict:
			chrom_dict[key] = "".join(chrom_dict[key])
	return chrom_dict

def write_to_files(chrom_dict):
		for key in chrom_dict:
			f = open(key+".fasta", "w")
			f.write(">"+key+"\n")
			f.write(chrom_dict[key])
			f.close()

def main():
	reference = argv[1]
	chrom_dict = get_chrom(reference)
	write_to_files(chrom_dict)
	
if __name__ == "__main__":
	main()
