#!/usr/envs/bin/python3
"""
script to get taxids from the salltitles line in blastx output

usage:
python3 acc_to_taxid.py [blastx_file]
"""
#import statements
from sys import argv


def get_acc(blastx_file):
	"""
	get accessions from blastx output file
	
	blastx_file: file, blastx output in tabular format
	format :'6 qseqid qstart qend salltitles evalue qframe pident qcovs sstart send slen'
	"""
	with open(blastx_file) as reader:
		acc_line_dict = {}
		count = 0
		for line in reader:
			line.strip()
			acc = line.split()[3]
			if "||" in acc:
				acc = acc.split("|")[2]
			elif "|" in acc:
				acc = acc.split("|")[1]
			else:
				acc = acc
			if acc_line_dict.get(acc) == None:
				acc_line_dict[acc] = [line]
			else:
				acc_line_dict[acc].append(line)
	return acc_line_dict

def get_taxid(acc_line_dict, ncbi_file):
	"""
	loop through protein2taxid file and match the accessions with taxids
	
	acc_line_dict: dict;
						keys: 0 to len(blastx_file)
						values: accessions
	ncbi_file: file, protein2accessions file from NCBI
	"""
	taxid_line_dict = {}
	with open(ncbi_file, "r") as reader:
		count = 0
		for line in reader:
			line = line.strip()
			acc = line.split()[0]
			taxid = line.split()[1]
			if acc_line_dict.get(acc) == None:
				continue
			else:
				taxid_line_dict[acc] = taxid
				print(count)
				count+=1
				print(acc, taxid)
	return taxid_line_dict

def write_out(taxid_line_dict, outputfile, acc_line_dict):
	"""
	write output to file
	"""
	f = open(outputfile, "w")
	lineslist = []
	for key in taxid_line_dict.keys():
		for i in range(0,len(acc_line_dict[key])):
			lines = acc_line_dict[key]
			lineslist.append(str(lines[i].strip())+"\t"+taxid_line_dict[key]+"\n")
	for line in lineslist:
		f.write(line)
	f.close()

def main():
	"""
	main function of this script
	"""
	blastx_file = argv[1]
	ncbi_file = argv[2]
	outputfile = argv[3]
	acc_line_dict = get_acc(blastx_file)
	taxid_line_dict = get_taxid(acc_line_dict, ncbi_file)
	write_out(taxid_line_dict, outputfile, acc_line_dict)


if __name__ == "__main__":
	main()
