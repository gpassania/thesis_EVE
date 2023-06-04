#!/usr/bin/python3
"""
Author: Gian Passania

Script to select data from detection pipeline results

usage:

python3 detect_select.py {complete table}

where complete_table is the path to the output of the pipeline results named "completetable.txt"
"""
#import statements
from sys import argv

#functions

def select(complete_table):
	"""
	selection of results based on evalues and Fraction
	
	complete_table: txt file
	"""
	outputfilename = complete_table.split("/")[-3]
	outpath = complete_table.split("/")[:-1]
	outpath.append(outputfilename)
	f = open("/".join(outpath)+"_extracts.txt", "w")
	f.write("Scaffold\tStart\tEnd\tLength\tID\tFraction\tVirusName\tvirus_eval\tPercentage_ID\tSubject_eval\n")
	with open(complete_table) as reader:
		for line in reader:
			infolist=[]
			line = line.strip()
			if line.startswith("ID"):
				continue
			parts = line.split("\t")
			if parts[8] == "0.0" or parts[8] == "0":
				continue
			if float(parts[11]) <= float(parts[22]):
				infolist.append(parts[1])
				infolist.append(parts[2])
				infolist.append(parts[3])
				infolist.append(parts[4])
				infolist.append(parts[5])
				infolist.append(parts[8])
				infolist.append(parts[10])
				infolist.append(parts[11])
				infolist.append(parts[14])
				infolist.append(parts[22])
				infolist.append("\n")
				info = "\t".join(infolist)
				f.write(info)
		f.close()


def main():
	complete_table = argv[1]
	select(complete_table)

if __name__ == "__main__":
	main()
