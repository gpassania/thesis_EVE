Vy-PER - Virus integration detection bY Paired End Reads 

Vy-PER helps to identify reads from paired-end sequencing, where one 
end aligns to host DNA (e.g. a human genome), and the other end aligns 
to e.g. a virus genome (or any other non-host genome).
 
Paired-end sequencing refers to Illumina HiSeq whole genome, exome or
transcriptome sequencing with typically 2x50bps to 2x101bps.
Vy-PER should also work with longer reads, if BWA is executed in
long-read mode.

Version: v0.3.5, released on 10 February 2015

Authors: 
Michael Forster m.forster@ikmb.uni-kiel.de and 
Silke Szymczak s.szymczak@ikmb.uni-kiel.de


1) Overview:
=====================

There are four scripts (Python, R) and two worked example
scripts in this download. You need to install a few dependencies: 
BLAT, BWA, Phobos, Python with Pysam, R, and SAMtools, see detailed
Dependencies section, below.  

1.1 Tool scripts for Linux cluster:

Vy-PER_sam2fas_se (can be used on a per-lane basis)
Vy-PER_blatsam (can be used on a per-lane basis)
Vy-PER_final_filtering (can merge the per-lane files)
rscript_ideogram.R (used in Vy-PER_final_filtering, also useable as
stand-alone plot script with user-defined virus names and colors)

Instead of a single script, we have split up the scripts into
several tools, because whole genomes usually run on several lanes
(e.g. 3-8 lanes, depending on coverage requirements). This allows
users to align the reads on a per-lane basis, and after alignment
of a lane is complete on a 16-core node, the first two Vy-Per scripts
can be run on this lane's alignment files on a single core. This per-lane 
approach results in better parallelization and faster run times.
Prior duplicate removal or recalibration/realignment with GATK
are not necessary. The Vy-PER jobs can run on the very first raw
per-lane BWA-alignment. In parallel to the Vy-PER jobs, the computationally 
very expensive classic realignment/recalibration/lane-merging/GATK variant
calling jobs can be submitted to other Linux nodes. 

1.2 Optional tool scripts for RIVYERA (FPGA) computer:

smith_waterman.v0.2.sh
smith_waterman_resultscopy.sh

If you do not have a RIVYERA, set rivyera=0 in the example scripts (see 1.3).
If you have a RIVYERA, you need to copy these two scripts into your user folder
on the RIVYERA:
/home/${rivuser}/scripts/generic/smith_waterman.sh
/home/${rivuser}/scripts/generic/smith_waterman_resultscopy.sh

To enable automated file transfer between your Linux cluster and the RIVYERA,
you need to create two pairs of rsa keys and put the public keys of each onto
the other server:
 
a) Create a public ssh key on the cluster headnode:
 ssh-keygen
  Generating public/private rsa key pair.
  Enter file in which to save the key (/home/sukmb175/.ssh/id_rsa): /home/sukmb175/.ssh/my_new_key.rsa
  Enter passphrase (empty for no passphrase): 
  Enter same passphrase again: 
  Your identification has been saved in my_new_key.rsa.
  Your public key has been saved in my_new_key.rsa.pub.

  Append this public ssh key on the cluster headnode to /home/myuser/.ssh/authorized_keys:
  cat my_new_key.rsa.pub >> .ssh/authorized_keys

b) Create another public ssh key on the cluster headnode and copy the public ssh key 
 to the_rivyera_server:/home/myuser/.ssh/your_cluster_public_key (you can also
 append it to the_rivyera_server:/home/myuser/.ssh/authorized_keys, but you MUST 
 specify the /home/myuser/.ssh/your_cluster_public_key in the shellscript
 smith_waterman_resultscopy.sh, in line "rsafile=/home/mforster/.ssh/nukey"


1.3 Example scripts:

vyper_example_script_L526401A.v0.3.4.sh (run time about 20 minutes)
vyper_example_script_L526401A_chr11.v0.3.4.sh (run time < 4 minutes)
vyper_example_script_198T.v0.3.4.sh (run time about 20 minutes)
vyper_example_script_268T.v0.3.4.sh (run time about 20 minutes)

You can download the example transcriptome data L526401A to make sure that
the scripts are working. The small example data (chr11) runs in
about 4 minutes on a single core, and the full L526401A RNA Seq data runs in 
about 18 minutes.

The above example shell scripts will give you randomly colored
virus integrations in the PDF plots. You can also quickly run just the
plotting script rscript_ideogram.R and optionally specify a custom annotation
file in which you can define the colors and the spelling of the virus names.
See the example annotation file in the Vy-PER download:
info_annotation_virus.txt


2) Dependencies:
=================

2.1 Python 2.7
https://www.python.org/downloads/
Example:
wget http://www.python.org/ftp/python/2.7.3/Python-2.7.3.tgz
tar -xvzf Python-2.7.3.tgz
cd Python-2.7.3
./configure
make
#pwd
#/home/michael/progs/external/Python-2.7.3
# required for the packages such pysam, which are installed in local user directory:
mkdir /home/michael/progs/external/lib/python2.7/site-packages -p
cd ~
vim .bashrc
## Python 2.7.3 ##
PATH=/usit/abel/u1/andrfra/progs/external/Python-2.7.3:$PATH
# for installing python packages:
export PYTHONPATH=/home/michael/progs/external/lib/python2.7/site-packages/

2.2 Pysam 0.6 for Python 2.7 
http://code.google.com/p/pysam/downloads/list
Example:
wget http://pysam.googlecode.com/files/pysam-0.6.tar.gz
tar -xvzf pysam-0.6.tar.gz
cd pysam-0.6
#pwd
#/home/michael/progs/external/Python-2.7.3/pysam/pysam-0.6
python setup.py build
python setup.py install --prefix /home/michael/progs/external

2.3 R 3.0.1
http://www.r-project.org
In R, you need to install the Bioconductor package quantsmooth:
source("http://bioconductor.org/biocLite.R")
biocLite("quantsmooth")
(see http://www.bioconductor.org/packages/release/bioc/html/quantsmooth.html)

2.4 BLAT
WARNING: Commercial users must obtain a licence from Jim Kent.
http://genome.ucsc.edu/FAQ/FAQblat.html#blat3

2.5 BWA
http://sourceforge.net/projects/bio-bwa/files/

2.5 SAMtools
http://sourceforge.net/projects/samtools/files/

2.6 Phobos 3.3.12
WARNING: Commercial users must obtain a licence from Christoph Mayer.
http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm
 
2.7 Reference sequences from Vy-PER homepage
wget http://www.ikmbtmp.uni-kiel.de/pibase/vy-per/vy-per_refs.tar.gz
(You will need to rebuild the BWA reference sequence files
if your version of BWA is not compatible with 0.6.2-r126)


3) Example data download:
==========================

3.1 Example 1: Full transcriptome (run time about 20 minutes on single core)
Download the pre-aligned sam/bam-files (and un-tar):
wget http://ikmbtmp.uni-kiel.de/pibase/vy-per/L526401A.tar.gz
tar -xvzf L526401A.tar.gz

The original raw data are from:
wget http://odin.mdacc.tmc.edu/~xsu1/L526401A_1.fq.gz
wget http://odin.mdacc.tmc.edu/~xsu1/L526401A_2.fq.gz
Liver cancer example data from the VirusSeq homepage (VirusSeq 
is NOT needed for the Vy-PER example script, but if you don't know
VirusSeq, why not try it too? http://odin.mdacc.tmc.edu/~xsu1/VirusSeq.html,
VirusSeq was published by Chen et al, PMID: 23162058)

3.2 Example 2: chr11 of transcriptome (run time 4 minutes)
wget http://www.ikmbtmp.uni-kiel.de/pibase/vy-per/L526401A.chr11_and_250k_unmapped.tar.gz
tar -xvzf L526401A.chr11_and_250k_unmapped.tar.gz

3.3 Example 3: half-lane of liver cancer genome (run time ca. 20 minutes)
published by Sung et al, PMID: 22634754.
wget http://ikmbtmp.uni-kiel.de/pibase/vy-per/198T.tar.gz
tar -xvzf 198T.tar.gz

3.4 Example 4: half-lane of liver cancer genome (run time ca. 20 minutes)
published by Sung et al, PMID: 22634754.
wget http://ikmbtmp.uni-kiel.de/pibase/vy-per/268T.tar.gz
tar -xvzf 268T.tar.gz


4) Brief summary of the scripts:
=================================

Full command line help is displayed when the script name
is entered without any further command line parameters.

2.1 Vy-PER_sam2fas_se

This Python script can be used on a per-lane basis. Features:
a) extract unmapped ends from a paired-end SAM file,
b) optionally filter these reads using a Phobos-3-file, and 
c) save these reads in FASTQ or FASTA format.


2.2 Vy-PER_blatsam

This Python script can be used on a per-lane basis. Features:

a) for each read that seems to map to a virus, get the paired
end that maps to the human genome (from the paired-end SAM file)
b) for each read that seems to map to a virus, get the virus
or - if there are ambiguities - the top 3 viruses
c) output a detailed tab-separated summary text file. 
d) Optionally output a folder of fasta files for interactive
QC using other alignment tools.


2.3 Vy-PER_final_filtering

This Python script can merge the per-lane files. Features:
a) Final filtering
b) Filtered summary tables (tab-separated text files)
c) PDF summary plot using rscript_ideogram.R


2.4 rscript_ideogram.R

This R script generates a summary plot (PDF) for virus integrations on
chr1-chr22, chrX, chrY. Features:
a) reads the .hit_clusters.txt file from Vy-PER_final_filtering
b) optionally reads virus name aliases and virus colors from extra file
c) generates summary plot (PDF) of virus/hg19 chimera paired-ends
d) Not summarised in the plot are any integration candidates detected on
unmapped contigs or mtDNA (chrM).


5) Disclaimer:
===============

The software, documentation and any help or support are provided "as-is"
and without any warranty of any kind. Installation of the required dependencies
must be carried out by a competent linux user or administrator at their own risk.

By downloading or using the Vy-PER tools, the Non-Commercial User accepts these
unmodified terms and conditions. 

The term "Software" in these terms and conditions refers to the Vy-PER software
components, including upgrades, additions, copies and documentation as supplied
by the Institute of Clinical Molecular Biology. 
The term "Non-Commercial User" in these terms and conditions refers to the person
who uses the Software for non-commercial purposes. Commercial use is prohibited unless
expressly permitted in writing by the authors of BLAT (Jim Kent, UCSC, USA) and Phobos
(Christoph Mayer, Ruhr-Uni Bochum). 
The term "person" includes the organisation or group if the Software is used within
the activities of an organisation or group. 
The Institute of Clinical Molecular Biology grants to the Non-Commercial User a
non-exclusive license to use the Software free of charge. The Institute of Clinical
Molecular Biology grants to the Non-Commercial User a non-exclusive license to use,
store, publish and exploit free of charge and without any time limit any results which
the Non-Commercial User has generated with the help of the Software, provided that the
Non-Commercial User cites the Software and version, the web site www.ikmb.uni-kiel.de/vy-per
and the paper on Vy-PER. 

The Non-Commercial User acknowledges that the Software is a tool to be used by the
Non-Commercial User in his or her work and that the Non-Commercial User shall remain
solely responsible for the content and quality of any results produced by the Non-Commercial
User from the use or with the aid of the Software. The Non-Commercial User undertakes to
indemnify and hold the Institute of Clinical Molecular Biology harmless in the event of
any claim against the Institute of Clinical Molecular Biology arising from any incorrect
or misleading information that the Non-Commercial User has obtained with the aid or by the
use of the Software. 

These terms and conditions shall be governed by and construed and interpreted in accordance
with German Law and submitted to the sole jurisdiction of the court of Kiel (Germany). 


5) Citing:
===========

When publishing results obtained with the help of Vy-PER, please cite
a) the software: Vy-PER v0.3.5, 
b) the web site: www.ikmb.uni-kiel.de/vy-per
c) the paper:
Forster M*, Szymczak S*, Ellinghaus D, Hemmrich G, Ruehlemann M, Kraemer L, Mucha S, Wienbrandt L,
Stanulla M for the UFO Sequencing Consortium within the I-BFM Study Group and
Franke A. 
Vy-PER: eliminating false positive detection of virus integration events in next generation sequencing data. 
(Original submission date to Scientific Reports: 14 Oct 2014, revised manuscript submitted 10 Feb 2015)
