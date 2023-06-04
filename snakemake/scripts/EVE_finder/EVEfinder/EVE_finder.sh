#!/bin/bash

#$ -l mem_free=5G
#$ -l h_rt=0:10:0
#$ -cwd
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -pe smp 8
#$ -R yes

######			######			######			######			######			######			######

######			Definitions				######

path=$1
file_list=$2
Database=$3
input=$(sed -n "${SGE_TASK_ID}p" "${path}${file_list}")
input=${input/.fa}

######			######			######			######			######			######			######

######			EVE Identification			######

### [Blast #1] Using Blastx to search a database of viral proteins using the genome of interest as a query ###

blastx \
-query ${path}${input}.fa \
-out ${path}${input/.fa}.result \
-evalue 1e-6 -num_threads $NSLOTS \
-db /andino/hpse1/home/mkunitomi/Refseq_Libraries/${Database} \
-outfmt '6 qseqid qstart qend salltitles evalue qframe pident qcovs sstart send slen'

######			######			######			######			######			######			######

######			Top Hit Finder			######

### [Blast to Bed] Convert blast output to BED format ###

EVE_finder/Blast_to_Bed3.py ${input}.result

### [Bed Sort] sort the bed file (required for downstream steps) ###

sortBed -i ${input}_EVEs.bed > ${input}_sorted_EVEs.bed

### [Index Genome Bedtools] Create a genome index for Bedtools ###

EVE_finder/Index_Genome_Bedtools.sh &>/dev/null

### [Bed Merger] Use Bedtools merge to merge overlapping hits ###

mergeBed -i ${input}_sorted_EVEs.bed -c 4,5,6,7,8,9,10,11 -o collapse,collapse,distinct,collapse,collapse,collapse,collapse,collapse > ${input}_merge_EVEs.bed

### [Top Hit MergedBed] Select the top hit from the merged bed file. ###

EVE_finder/Top_score_BED2.py ${input}_merge_EVEs.bed ${input}_merged_EVEs.bed

### [Bed_to_Fasta] Pull hit sequence from the genome and create a multifasta file ###

bedtools getfasta -s -name -fi ${input}.fa -bed ${input}_merged_EVEs.bed -fo ${input}_merged_EVEs.fasta

###	Cleanup	###
rm ${input}_sorted_EVEs.bed
rm ${input}_merge_EVEs.bed
######			######			######			######			######			######			######

######			False Postive Blast Search			######

### [Blast #1] Using Blastx to search the RefSeq NR protein database using the hit fasta file as a query ###

blastx \
-query ${path}${input}_merged_EVEs.fasta \
-out ${path}${input}_BLAST_NR \
-db /andino/hpse1/home/mkunitomi/Refseq_Libraries/refseq_protein \
-evalue 1e-4 -num_threads $NSLOTS \
-num_alignments 100 \
-outfmt '6 qseqid qstart qend salltitles evalue qframe pident qcovs sstart send slen staxids'

######			######			######			######			######			######			######

######			False Postive Filter			######

### [BLAST_NR_filter] Removes all hypothetical, predicted, and genus specific results	###

EVE_finder/BLAST_NR_filter.py ${input}_BLAST_NR

### [TAX_ID] For each ID for each hit, checks to see if it is viral. If so, add the hit to list of validated hits in Fasta format ###

EVE_finder/TAX_check4.py ${input}_BLAST_NR_filtered ${input}_merged_EVEs.fasta ${input}_merged_EVEs.bed

######			######			######			######			######			######			######

###		END OF SCRIPT		###















