#!/bin/usr/python3

"""
Author: Gian Passania
script for filtering results from vyper blatsam, removing low coverage hits

usage:
python3 vyper_pre_filter.py [blatsam_file] [new_file]
"""

#import statements
from sys import argv

#functions

def filter_blatsam(blatsam_file, new_file):
	"""
	filter blatsam file based on mapQV and length
	
	blatsam: str, full path to blatsam file
	new_file: str, full path for new filtered file
	"""
	f = open(new_file, "w")
	with open(blatsam_file) as reader:
		for line in reader:
			if "RunID" in line:
				f.write(line)
			if line.startswith("#"):
				continue
			parts = line.split("\t")
			if int(parts[5]) >= 25 and int(parts[15]) >= 125:
				newline = "\t".join(parts)
				f.write(newline)
	f.close()
	
def main():
	blatsam_file = argv[1]
	new_file = argv[2]
	filter_blatsam(blatsam_file, new_file)


if __name__ == "__main__":
	main()

