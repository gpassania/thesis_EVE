#!/usr/bin/env/python3

"""
author: Gian Passania
script to write paths to snakemake configuration file

usage:
python3 add_configs.py [reference] [database] [fastq] [script path] [protein database]

reference: full path to representative genome, example: /home/reference/species/representative_genome.fasta
database: full path to nucleotide database, example: /home/database/RVDB.fasta
fastq: full path to fastq files, example, example: /home/fastqfiles/species/
script_path: full path to script folder, example: /home/snakemake/scripts
protein_database: full path to protein database, example: /home/database/RVDB_prot.fasta
"""

#import statements
from sys import argv
import os

#functions
def write_to_config(reference, database, database_prot, fastq, script_path) :
	"""
	writes paths given to config.yaml file
	
	reference: str, path to reference genome fasta directory
	database: str, path to database fasta file
	fastq_1: str, path to first fastq file
	fastq_2: str, path to second fastq file 
	script_path: str, path to scripts directory
	"""
	mainpath = reference.split("/")[:-3]
	species = reference.split("/")[-2]
	reference_path = "/".join(reference.split("/")[:-1])
	f = open("config.yaml", "w")
	f.write("species:\n    {}\n".format(species))
	f.write("mainpath:\n    {}/\n".format("/".join(mainpath)))
	f.write("database:\n    {}\n".format(database))
	f.write("scripts_path:\n    {}\n".format(script_path))
	f.write("database_prot:\n    {}\n".format(database_prot))
	f.write("representative_ref:\n    {}\n".format(reference))
	f.write("reference_path:\n    {}\n".format(reference_path))
	f.write("fastq_files_path:\n    {}\n".format(fastq))
	f.write("\n")
	f.write("refgenomes:\n    ")
	for dirs in os.listdir(reference_path):
		if dirs.endswith(".fna"):
			f.write("{}: {}{}\n    ".format(dirs, reference_path+"/", dirs))
	f.write("\n")
	f.write("fastq_names:\n    ")
	for files in os.listdir(fastq):
		f.write("{}: {}{}/{}\n    ".format(files, fastq, files, files))
	f.close()
	
def main():
	"""
	main function of script
	"""
	reference = argv[1]
	database = argv[2]
	fastq = argv[3]
	script_path = argv[4]
	database_prot = argv[5]
	write_to_config(reference, database, database_prot, fastq, script_path)
if __name__ == "__main__":
	main()
