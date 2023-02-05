#!/usr/bin/env/python3

"""
author: Gian Passania
script to write paths to snakemake configuration file
"""

#import statements
from sys import argv
import os

#functions
def write_to_config(reference, database, fastq_1 = "test", fastq_2 = "test2"):
	"""
	writes paths given to config.yaml file
	
	reference: str, path to reference genome fasta file
	database: str, path to database fasta file
	fastq_1: str, path to first fastq file
	fastq_2: str, path to second fastq file 
	"""
	mainpath = reference.split("/")[:-1]
	refpath = "/".join(mainpath) + "/split_ref"
	split_files = os.listdir(refpath)
	f = open("config.yaml", "w")
	f.write("mainpath: {}/\n".format("/".join(mainpath)))
	f.write("reference: {}\n".format(reference))
	f.write("database: {}\n".format(database))
	f.write("fastq_1: {}\n".format(fastq_1))
	f.write("fastq_2: {}\n".format(fastq_2))
	f.write("samples:\n")
	for i in range(len(split_files)):
		f.write(" {}: {}\n".format(split_files[i],split_files[i]))
	f.close()
	
def main():
	reference = argv[1]
	database = argv[2]
	write_to_config(reference, database)
if __name__ == "__main__":
	main()
