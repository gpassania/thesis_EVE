#!/bin/bash

############################################################################
# Generic Rivyera shell script to submit the smith-waterman run
# and create a "signal" file "*.done" when the run has completed 
############################################################################
# Michael Forster                              m.forster@ikmb.uni-kiel.de 
# v0.3.3   2 September 2014
############################################################################

if [ $# -ne 6 ]
then
echo "You must give six parameters: <rivyera_user> <project> <subproject> <library> <ref> <target>"
echo " e.g. smith_waterman.sh mforster ibd_ped test A0503 /home/mforster/data/refs/hg19/hg19.complete.fa user@server.uni-kiel.de:/targetfolder"
exit
fi

############################################################################
###  project, subproject, library ID, reference sequence 
###  from command line arguments
############################################################################

rivuser=${1}
proj=${2}
subp=${3}
lib=${4}
ref=${5}
targetlocation=${6}

############################################################################
###  put together the standardized input/output file names
############################################################################

infile=/home/${rivuser}/data/${proj}/${subp}/${lib}/${lib}_out.fpga.fa
outfile=/home/${rivuser}/data/${proj}/${subp}/${lib}/${lib}_smith-waterman
samfile=/home/${rivuser}/data/${proj}/${subp}/${lib}/${lib}_smith-waterman.sam
thecommand="/opt/se_software/se_sw/se_sw -df ${ref} -bit /opt/se_software/se_sw/bitfiles/ -d 1 -rc -qf ${infile} -o ${outfile} -so ${samfile} -sm 2"
logfile=/home/${rivuser}/data/${proj}/${subp}/${lib}/${lib}_smith_waterman.sh.log
donefile=/home/${rivuser}/data/${proj}/${subp}/${lib}/${lib}_out.fpga.fa.done
detailedfile=/home/${rivuser}/data/${proj}/${subp}/${lib}/${lib}_smith-waterman_detailed.txt

postcommand="/home/${rivuser}/scripts/generic/smith_waterman_resultscopy.sh ${detailedfile} ${donefile} ${targetlocation}"

echo $infile
echo $outfile
echo $samfile
echo $logfile
echo $thecommand

############################################################################
###  submit the job
############################################################################

jobctl submit --prio 10 --log ${logfile} --pre date --post "${postcommand}" "${thecommand}"


#Usage:  /opt/se_software/se_sw/se_sw [OPTIONS]
#
#   Options:   -m   MACHINE      Index of the machine to use for alignment.
#              -t   THRESHOLD    The Threshold of the alignment. Only queries with scores
#                                higher than this threshold will be reported. Default is 0.
#              -d   BEST         Align queries only against the BEST best fitting databases.
#              -max MAX          Specify which maximum to take. Choosable values are:
#                                     RAND   - Reports one maximum by random (Default).
#                                     FIRST  - Reports the leftmost maximum.
#                                     ALL    - Reports every maximum.
#              -rc               Add all reads as their reverse complement, too
#              -mg  MATE_GAP     Specifies a maximum gap between paired-end mappings.
#              -mo  MATE_OVERLAP Specifies a maximum overlap between paired-end mappings.
#              -so  OUTPUT_FILE  Specify a filename where the SAM output should be written to
#              -o   OUTPUT_FILE  Specify a filename where the output should be written to
#              -qc  MAX_Q        Specify the maximum number of query sequences
#              -dc  MAX_D        Specify the maximum number of database sequences
#              -qf  QUERY_FILE   Specify one or more FASTA/FASTQ files containing query sequences
#              -qmf MATE_FILES   Specify a multiple of two FASTA/FASTQ files containing mate-paired query sequences
#              -qs  QUERY_SKIP   Specify a number of queries of the FASTA File to skip
#              -df  DB_FILE      Specify one or more FASTA/NCBI files containing database sequences
#              -do  OUTPUT_FILE  Specify a file where Database references should be put (in FASTA Format).
#              -bit DIRECTORY    Specify the folder to search for bitfiles. Default is "bitfiles".
#              -smf SM_FILE      Specify the scoring matrix file to use for alignment.
#              -sm  INDEX        Specify the scoring matrix to use for alignment.
#                                   0: Default Scoring Matrix  (Nucleotide)
#                                   1: NUC22 Scoring Matrix    (Nucleotide)
#                                   2: NUC44 Scoring Matrix    (Nucleotide, full IUPAC)
#                                   3: BLOSUM62 Scoring Matrix (Proteins)



