#!/usr/bin/python3
"""
Author: Gian Passania

Script to further filter selected data obtained by detect_select.py

usage:

python3 detect_filter_final.py {extracts file}
"""
#import statements
from sys import argv

#functions

def get_IDs(extracts_file):
	"""
	for selected results where viral and subject evalues == 0, select best hit
	
	extracts_file: tab seperated txt file containing selected results
	"""
	equale_list = []
	normal_list = []
	with open(extracts_file, "r") as reader:
		for line in reader:
			parts = line.split("\t")
			if parts[7] == parts[9]:
				info = [parts[4], parts[6]]
				equale_list.append(info)
			else:
				normal_list.append(line)
	return equale_list, normal_list

def check_blastx(equale_list, blastx_file):
	"""
	for results with equal evalues, check blastx output to find best hit
	
	equale_list: list of IDs with equal evalues
	"""
	line_list = []
	IDs = []
	keep = []
	with open(blastx_file, "r") as reader:
		for line in reader:
			for ID, name in equale_list:
				if ID not in IDs:
					if ID in line:
						IDs.append(ID)
						info = [ID, line]
						line_list.append(info)
						continue
					else:
						continue
	for ID, name in equale_list:
		for line in line_list:
			if ID in line[0]:
				if name in line[1]:
					keep.append(ID)
	print(keep)
	return keep
	
def write_out(keep, extracts_file, normal_list):
	"""
	write a new file of filtered results
	
	keep: list of IDs to keep in results
	extracts_file: tab seperated txt file containing results
	normal_list: list of lines passing origional filter
	"""
	f = open(extracts_file[:-4]+"_final.txt", "w")
	for line in normal_list:
		f.write(line)
	with open(extracts_file, "r") as reader:
		for ID in keep:
			for line  in reader:
				if ID in line:
					f.write(line)
					break
	f.close()
		
def main():
	"""
	main function
	"""
	extracts_file = argv[1]
	blastx_file = argv[2]
	equale_list, normal_list = get_IDs(extracts_file)
	keep = check_blastx(equale_list, blastx_file)
	write_out(keep, extracts_file, normal_list)

if __name__ == "__main__":
	main()
