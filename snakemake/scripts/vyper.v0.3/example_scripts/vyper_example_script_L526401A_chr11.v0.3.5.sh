#!/bin/bash
#PBS -r y
#PBS -j oe
#PBS -q exc_b2
#PBS -M m.forster@ikmb.uni-kiel.de
##PBS -l select=1:exc_s=true:ncpus=8:mem=31gb
##PBS -l select=1:ncpus=16:mem=120gb
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N vyper_chr11_L52
#PBS -o 'vyper_example_script_L526401A_chr11.v0.3.5.sh.output'

# Get Job Number (because PBS often loses its PBS_JOBID)
pbs_number=$(cut -d '.' -f 1 <<< ${PBS_JOBID})

# echo date and Job ID (via TMPDIR)
date
echo $TMPDIR

############################################################
## Example script, 1-lane of whole genome
## Run-time on 1-core node (v0.3.3): < 4 minutes
############################################################

############################################################
## Dear Vy-PER test User, please modify:
## - the obove PBS-Pro-header lines to suit your cluster
##   and job scheduler (you don't need 120gb or 16 cpus),
## - (A), (B), and (C) to suit your folder structure.
############################################################
## Prerequisites and dependencies:
##
## (1) Download the chromosome-11-only example file:
## wget http://www.ikmbtmp.uni-kiel.de/pibase/vy-per/L526401A.chr11_and_250k_unmapped.fq.tar.gz
## and extract:
## tar -xvzf L526401A.chr11_and_250k_unmapped.fq.tar.gz
##
## (We extracted the chr11 reads from the liver cancer example 
##  data on the VirusSeq homepage:
##  http://odin.mdacc.tmc.edu/~xsu1/VirusSeq.html,
##  published by Chen et al, PMID: 23162058.
##  VirusSeq is NOT needed for this Vy-PER example script, but
##  if you don't know VirusSeq, why not try it too?)
##
## (2) Python 2.7 with pysam 0.6 
## http://code.google.com/p/pysam/downloads/list
##
## (3) R 3.0.1
## http://www.r-project.org
## In R, you need to install two packages:
## a) the Bioconductor package quantsmooth:
## source("http://bioconductor.org/biocLite.R")
## biocLite("quantsmooth")
## (see http://www.bioconductor.org/packages/release/bioc/html/quantsmooth.html)
## b) the CRAN package RColorBrewer:
## install.packages("RColorBrewer")
## (see http://cran.r-project.org/web/packages/RColorBrewer/index.html)
##
## (4) BLAT
## http://genome.ucsc.edu/FAQ/FAQblat.html#blat3
##
## (5) BWA and SAMtools
## http://sourceforge.net/projects/bio-bwa/files/
## http://sourceforge.net/projects/samtools/files/
##
## (6) Phobos 3.3.12
##  http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm
## 
## (7) the reference sequences from the Vy-PER download page
## You will need to rebuild the BWA reference sequence files
## if your version of BWA is not compatible with 0.6.2-r126
##
## (8) If you have access to a RIVYERA, see the instructions
## in section (D). If not, you need to go into (D) and set
## riyvera=0.
############################################################
## For Vy-PER questions contact Michael Forster
## m.forster@ikmb.uni-kiel.de
############################################################

############################################################
## (A) Max number of threads allowed in this PBS job
############################################################

threads=1

############################################################
# (B) Base folders and reference sequences
############################################################

# FASTQ sequences
c=/ifs/data/nfs_share/sukmb175/sequences/publication_testdata/vyper/virusseq_testdata_chr11

# results
d=/ifs/data/nfs_share/sukmb175/projects/publication_testdata/vyper/virusseq_testdata_chr11.v0.3.5
mkdir -p ${d}

# host genome reference sequences for BWA:
ref=/ifs/data/nfs_share/sukmb175/references/hg19/bwa/hg19.complete.fa
reffai=${ref}.fai

# host genome reference sequences for BLAT:
blatref=/ifs/data/nfs_share/sukmb175/references/hg19/blat/hg19.complete.2bit

# viruses' genome reference for BLAT:

# The following genomes can be used for detecting only HIV integration loci:
# refvb=/ifs/data/nfs_share/sukmb175/references/viruses/2015_01_23/blat/hiv.2bit

# The following genomes can be used to detect 5801 different virus species or strains:
refvb=/ifs/data/nfs_share/sukmb175/references/viruses/2015_01_23/blat/all_viruses.2bit

# The following genomes were used for the published manuscript and do not contain HIV genomes:
# refvb=/ifs/data/nfs_share/sukmb175/references/viruses/2012_11_07/blat/all_viruses.2bit
refvbfa=${refvb}.fa
refvbfai=${refvbfa}.fai

############################################################
## (C) FASTQ input files and title of summary PDF
############################################################

# WARNING: The BWA-alignment is the prerequisite both for
# classical variant-calling and for Vy-PER. 
# If you unnecessarily duplicate this BWA-alignment
# before running Vy-PER, it will cost several hours of 16-core
# run time or several DAYS of SINGLE core run time.

# set variable NOT to perform the BWA-alignment but to 
# proceed directly to the Vy-PER steps:
nobwa=1

lib=L526401A
fastq1=${c}/${lib}_1.fq.gz
fastq2=${c}/${lib}_2.fq.gz

readsize=101
center=ICMB_Kiel_Germany

# title of summary PDF
plottitle='Liver cancer L526401A chr11 (RNA)'
echo Plot title:
echo ${plottitle}

############################################################
## (D) RIVYERA specific things
############################################################

# set variable to use the RIVYERA:
rivyera=1

##### RSA keys need to be generated manually, once #####
# Create two pairs of rsa keys:
# 
# Create a public ssh key on the cluster headnode:
# ssh-keygen
#   Generating public/private rsa key pair.
#   Enter file in which to save the key (/home/sukmb175/.ssh/id_rsa): /home/sukmb175/.ssh/my_new_key.rsa
#   Enter passphrase (empty for no passphrase): 
#   Enter same passphrase again: 
#   Your identification has been saved in my_new_key.rsa.
#   Your public key has been saved in my_new_key.rsa.pub.
#
# Append this public ssh key on the cluster headnode to /home/myuser/.ssh/authorized_keys:
#   cat my_new_key.rsa.pub >> .ssh/authorized_keys
#
# Create another public ssh key on the cluster headnode and copy the public ssh key 
# to the_rivyera_server:/home/myuser/.ssh/your_cluster_public_key (you can also
# append it to the_rivyera_server:/home/myuser/.ssh/authorized_keys, but you MUST 
# specify the /home/myuser/.ssh/your_cluster_public_key in the shellscript
# smith_waterman_resultscopy.sh, in line "rsafile=/home/mforster/.ssh/nukey"
#####

# PBS headnode:
headnode=rzcl00g

# cluster headnode:
clusthead=ikmbhead.rz.uni-kiel.de

# user account on cluster headnode:
clustuser=sukmb175

# RIVYERA user account:
rivuser=mforster

# RIVYERA server address:
rivserver=rivyera.ikmb.uni-kiel.de

# project and subproject:
proj=publication_testdata
subp=vyper

# smith-waterman masterscript:
swm=/home/${rivuser}/scripts/generic/smith_waterman.sh

# hg19 reference on RIVYERA:
rivref=/home/${rivuser}/data/refs/hg19/hg19.complete.fa

############################################################
## END OF DEFINITIONS. THE REST IS THE SAME FOR ALL SAMPLES
############################################################


if [  $nobwa -ne 1 ]
then

	############################################################
	## BWA alignment as 16-thread job on on 16-core node 
	############################################################
	
	# output files:
	out1=${d}/${lib}_R1.sai
	out2=${d}/${lib}_R2.sai
	
	echo "1) BWA step 1: alignment of single ends"
	
	bwa aln \
	       -n 2 -q 15 -l 5000 -t ${threads} \
	       ${ref} ${fastq1} > ${out1}
	
	echo ${out1}
	date
	
	bwa aln \
	       -n 2 -q 15 -l 5000 -t ${threads} \
	       ${ref} ${fastq2} > ${out2}
	
	echo ${out2}
	date
	
	############################################################
	## BWA Pairing as single-core job 
	## (15 cores idle, but that's just to keep this example simple
	## and easy-to-run. In real life, you will split up the
	## real jobs to make optimal use of your linux cluster.)
	############################################################
	
	# output files:
	out3pe=${d}/${lib}_PE.sai
	out4=${d}/${lib}_PE
	s=${lib}_PE
	
	echo "2) BWA step 2: pairing"
	
	bwa sampe \
	-a 500 \
	-r "@RG\tID:Illumina\tSM:${lib}\tLB:lib_2x${readsize}\tDS:${ref}\tCN:${center}" \
	${ref} ${out1} ${out2} ${fastq1} ${fastq2} > ${out3pe}
	
	echo ${out3}
	date
	
	echo "3) sam2bam"
	
	samtools view -uhSt ${ref}.fai ${out3pe} | samtools sort - ${out4}
	echo SORTED BAM
	date
	
	samtools index ${out4}.bam
	echo INDEXED BAM
	date
else

	# define some file names referenced later:
	out3pe=${d}/${lib}_PE.sai
	out4=${d}/${lib}_PE

fi

############################################################
## Vy-Per part 1: extract reads and filter low-complexity reads
############################################################

# output folder:
e=${d}/oneEndUnmapped
mkdir $e

# output files:
out=${e}/${lib}_PE_oneEndUnmapped.sam
fas=${e}/${lib}_PE.fa
pho=${e}/${lib}_PE.phobos.3.txt
fasq=${e}/${lib}_PE_oneEndUnmapped.fastq

# extract unmapped read end into separate sam file:

echo "4) Vy-PER part 1: samtools extract unmapped end"
samtools view -f 4 -F 264 ${out4}.bam > ${out}
date

# convert into fasta using simple single-end python script
# because bam2fastq only converts read pairs into a pair of fastq files

echo "5) Vy-PER part 1: vyper_sam2fas_se"
Vy-PER_sam2fas_se ${out} ${fas} 
date

echo "6) Vy-PER part 1: phobos STR detection"
phobos --outputFormat 3 ${fas} ${pho}
date

echo "7) Vy-PER part 1: vyper_sam2fas_se filtering of reads with more than 30bp non-STR region"

# output folders:
f=${e}/non_str
vblat=${f}/virusMappedBLAT
mkdir -p $vblat

fas4=${f}/${lib}_oneEndUnmapped.non_str.fa
Vy-PER_sam2fas_se -fp3 ${pho} 30 ${out} ${fas4}
date

############################################################
## Vy-Per part 2: BLAT against virus genomes
############################################################

echo "8) Vy-PER part 2: BLAT unmapped reads to all viruses genomes"
outblat=${vblat}/${lib}_oneEndUnmapped.hg19sloppy.unmapped.blat_to_viruses.txt
blat -out=blast8 -noTrimA -t=dna -q=dna -maxGap=0 -fastMap ${refvb} ${fas4}  ${outblat} 
date

############################################################
## Vy-Per part 3: summary table 
############################################################

# output folder for per-hit fasta & reference file pairs 
# (e.g. for inspecting with an interactive aligner/viewer/editor, 
# such as DNA Alignment, www.fluxus-engineering.com/align.htm)

outfa=${vblat}/vyper_fasta_files
mkdir ${outfa}

cd ${vblat}

sam=${out3pe}
out5=${vblat}/${lib}_oneEndUnmapped.hg19sloppy.unmapped.blat_to_viruses.vyper_blatsam.txt

echo "9) Vy-PER part 3: vyper_blatsam summary file and folder of fasta files" 
Vy-PER_blatsam ${lib} ${refvbfa} ${outblat} ${fas4} ${ref} ${sam} ${out5} ${outfa}
date

############################################################
## Vy-PER part 4: final filtering
## (can merge multiple-lane lane results from Vy-Per parts 1-3)
## DEFAULT FINAL FILTERING
############################################################

outfilt=${vblat}/vyper_final_filtering_10_chimeras
workdir=${outfilt}/work
mkdir -p ${workdir}

# leave the tmp__* files here, in case you to inspect them:
cd ${workdir}

#### default command line parameters: 
# virus cluster region size in host genome: 1000 bps (as the DNA fragments
# for Illumina paired-end libraries are typically 300-650bps long),
# min 10 virus hits within cluster,
# min 10 hits from same virus within cluster,
# 1% Illumina-phiX174-library spike-in (is the usual fraction in WGS, WXS, WTS),
# max 50% STR or homopolymer fraction in a virus hit sequence
# if alternative virus hits were detected in the same host genome position: 
#   min length of 2nd (3rd) sequence = 95% of length of best hit
#   max 3bps length difference of 2nd (3rd) sequence to best hit
#   length_host_hit/length_virus_hit > 90% --> host hit

params="-p 1000 10 10 0.01 0.5 0.95 3 0.90"

# No FPGA-based Smith-Waterman option (is off by default)
sw=0
# output file from FPGA-based Smith-Waterman (just a placeholder here):
swout=${outfilt}/${lib}.out_detailed.txt
# min read-mapping quality to hg19
mq=0

# Final output files from Vy-Per_final_filtering:
outbase=${outfilt}/${lib}_out

echo "10) Vy-PER part 4: final filtering (can be repeated at this stage, or explored with different parameters"
echo " "
echo 'Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"'

Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
date

# Note that the hg19 blatref reference is not used anymore in this version of
# Vy-PER_final_filtering, because it takes DAYS to BLAT 1000 sequences to the entire human genome,
# at the required sensitive BLAT-settings. 
# So we now use the FPGA-based Smith-Waterman aligner which takes 2 minutes for the same file.

# However, as only we and a few others have a RIVYERA, you can re-activate the Blat-to-hg19 run
# in Vy-PER_final_filtering, using the following command line parameter for Vy-PER_final_filtering
# (this is only required for filtering singletons, see below HIGHLY SENSITIVE filtering):
sw=3

if [ $rivyera -ne 1 ]
then
	############################################################
	## Vy-PER part 4: Alternative to above default filtering:
	## HIGHLY SENSITIVE FINAL FILTERING WITHOUT SMITH-WATERMAN
	## (Only BLAT is used to align reads to a local genome window)
	############################################################
	
	outfilt=${vblat}/vyper_final_filtering_1_chimera
	workdir=${outfilt}/work
	mkdir -p ${workdir}
	
	# leave the tmp__* files here, in case you want to inspect them:
	cd ${workdir}
	
	#### altered command line parameters: 
	# min 1 virus hit within cluster,
	# min 1 hits from same virus within cluster,
	
	params="-p 1000 1 1 0.01 0.5 0.95 3 0.90"
	
	# No FPGA-based Smith-Waterman option (is off by default)
	sw=0
	# output file from FPGA-based Smith-Waterman (just a placeholder here):
	swout=${outfilt}/${lib}.out_detailed.txt
	# min read-mapping quality to hg19
	mq=20
	
	# Final output files from Vy-Per_final_filtering:
	outbase=${outfilt}/${lib}_out
	
	echo "11) Vy-PER part 4: final filtering (can be repeated at this stage, or explored with different parameters"
	echo " "
	echo Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
	
	Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
	date

else

	############################################################
	## Vy-PER part 4: Alternative to above default filtering:
	## HIGHLY-SENSITIVE FINAL FILTERING USING SMITH-WATERMAN FILE
	## (The Smith-Waterman file is a special formatted
	## output file from the RIVYERA-FPGA-Smith-Waterman aligner,
	## so most users will not be able to use this part 4)
	############################################################
	
	outfilt=${vblat}/vyper_final_filtering_smithwaterman
	workdir=${outfilt}/work
	mkdir -p ${workdir}
	
	### leave the tmp__* files here, in case you want to inspect them:
	cd ${workdir}
	
	### change command line parameters:
	### pick up each single virus candidate (including technical noise from phiX-chimera)
	params="-p 1000 1 1 0.01 0.5 0.95 3 0.90"
	
	############################################################
	# (a) export the fasta file ${outfilt}/${lib}_out.fpga.fa --> for FPGA-Smith-Waterman
	############################################################
	sw=1
	swout=${outfilt}/${lib}_smith-waterman_detailed.txt
	# min read-mapping quality to hg19
	mq=20
	### Final output files from Vy-Per_final_filtering:
	outbase=${outfilt}/${lib}_out
	echo "17) Vy-PER part 4: Alternative highly-sensitive final filtering STEP 1 - export ${outfilt}/${lib}_out.fpga.fa"
	echo " "
	echo Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
	Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
	date
	
	############################################################
	# (b) copy the exported fasta file to the RIVYERA:
	############################################################
	
	rivfolder=/home/${rivuser}/data/${proj}/${subp}/${lib}
	echo $rivfolder
	
	echo "ssh -T ${headnode} ssh -T ${rivuser}@${rivserver} mkdir -p ${rivfolder}"
	ssh -T ${headnode} ssh -T ${rivuser}@${rivserver} mkdir -p ${rivfolder}
	
	echo "ssh -T ${headnode} rsync -auv ${outfilt}/${lib}_out.fpga.fa ${rivuser}@${rivserver}:${rivfolder}/"
	ssh -T ${headnode} rsync -auv ${outfilt}/${lib}_out.fpga.fa ${rivuser}@${rivserver}:${rivfolder}/
	date
	
	############################################################
	# (c) submit the Smith-Waterman job on the RIVYERA
	############################################################
	
	# delete any existing flag-file from a previous run:
	rm -f ${outfilt}/${lib}_out.fpga.fa.done
	rm -f ${outfilt}/${lib}_smith-waterman_detailed.txt
	
	# build the target folder on this cluster server, from an external server's viewpoint:
	rivtarget=${clustuser}@${clusthead}:${outfilt}
	
	sleep 5
	
	ls -lh ${outfilt}/
	
	echo "ssh -T ${headnode} ssh -T ${rivuser}@${rivserver} ${swm} ${rivuser} ${proj} ${subp} ${lib} ${rivref} ${rivtarget}"
	ssh -T ${headnode} ssh -T ${rivuser}@${rivserver} ${swm} ${rivuser} ${proj} ${subp} ${lib} ${rivref} ${rivtarget}
	
	############################################################
	# (d) wait until Smith-Waterman job has completed
	############################################################
	
	# poll the results every second for the flag-file that
	# the job has finished and the results have completed copying 
	
	jobrunning=1
	while [  $jobrunning -eq 1 ]
	do
	  sleep 1
	  # If I copy files within a PBS job from the Rivyera, they are empty, so
	  # I now copy the data from the Rivyera to the PBS cluster.
	  # ssh -T ${headnode} rsync -auv ${rivuser}@${rivserver}:${rivfolder}/${lib}_out.fpga.fa.done ${outfilt}/
	  # ssh -t ${headnode} scp ${rivuser}@${rivserver}:${rivfolder}/${lib}_out.fpga.fa.done ${outfilt}/
	  if [ -a ${outfilt}/${lib}_out.fpga.fa.done ]
	  then
	     jobrunning=0
	  fi
	done
	
	# wait 5 seconds in case of operating system synchronisation delays:
	sleep 5
	
	ls -lh ${outfilt}/
		
	############################################################
	# (e) import the smith-waterman results file 
	# ${lib}_smith-waterman_detailed.txt and perform final filtering.
	############################################################
	sw=2
	swout=${outfilt}/${lib}_smith-waterman_detailed.txt
	# min read-mapping quality to hg19
	mq=20
	### Final output files from Vy-Per_final_filtering:
	outbase=${outfilt}/${lib}_out
	
	echo "18) Vy-Per part 4: Alternative highly-sensitive final filtering STEP 2 - with Smith-Waterman results file ${swout}"
	echo " "
	echo Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
	Vy-PER_final_filtering ${params} ${sw} ${swout} ${mq} ${out5} ${blatref} ${reffai} ${outbase} "${plottitle}"
	date

fi

############################################################
## PBS PRO summary stats (modify for your own job scheduler)
############################################################

echo " "
echo "-----------------------"
echo " "

qstat -f $pbs_number

