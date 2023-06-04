#!/bin/usr/python3

"""
Author: Gian Passania
script for seperating the "salltitles" line in blastx output and selecting the viral accession

usage:
python3 fix_titles.py [blastx_file] [new_file]
"""

#import statements
from sys import argv

#functions

def change_salltitles(blastx_file, new_file):
	"""
	check every salltitles and select the viral title
	
	blastx_file: str, full path to blastx output file
	new_file: str, full path to new output file
	"""
	f = open(new_file, "w")
	with open(blastx_file) as reader:
		for line in reader:
			line = line.strip()
			parts = line.split("\t")
			alltitles = parts[3]
			if "<>" in alltitles:
				if "virus" in alltitles:
					titles = alltitles.split("<>")
					for title in titles:
						if "virus" in title:
							parts[3] = title
			newline = "\t".join(parts)
			f.write(newline)
			f.write("\n")
	
def main():
	blastx_file = argv[1]
	new_file = argv[2]
	change_salltitles(blastx_file, new_file)


if __name__ == "__main__":
	main()

