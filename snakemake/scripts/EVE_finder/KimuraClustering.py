# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 12:41:09 2016

@author: zwhitfield

This script merges TEs with their respective Kimura scores (both as found in the .fa.align file)
"""

import sys
import pandas as pd
import numpy as np
import scipy.stats as sp
#scipy.__version__
#np.__version__
import matplotlib.pyplot as plt

inputdir = "/path/to/output/from/ExtractKimuraScores.sh/
outputdir = "/path/for/outputFileToBeVisualized/"


#--------------------------------------------------------------
#Load TE and Kimura data, then merge
#--------------------------------------------------------------

#Read in file of headers from .fa.align that were extracted using ExtractKimuraScores.sh
#These can have some duplicates (at least by how I am using the files to merge), such as:
#236 22.98 4.62 6.25 000000F 73596 73725 (7827976) C rnd-3_family-132#LTR/Copia (768) 512 385 m_b2s001i28 118
#These stem from the original .fa.align file
pathToTEs = inputdir + "/onlyHeader_filtered_From_Aag2_Contigs.fa.align.txt"
allTEs =  pd.read_csv(pathToTEs,
                      header = None,
                      names = ["SWscore","percDiv","percDel","percIns","querySeq","qPosBeg","qPosEnd","qLeft", "strand","matchingRepeat","rLeft","rPosBeg","rPosEnd","dunno","id"],
                      sep = '\s+'
                      )

#Explanation of columns according to http://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html
#1320 15.6  6.2  0.0  HSU08988  6563 6781 (22462) C  MER7A   DNA/MER2_type    (0)  337  104  20
#
#  1320     = Smith-Waterman score of the match, usually complexity adjusted
#        The SW scores are not always directly comparable. Sometimes
#        the complexity adjustment has been turned off, and a variety of
#        scoring-matrices are used dependent on repeat age and GC level.

#  15.6     = % divergence = mismatches/(matches+mismatches) **
#  6.2      = % of bases opposite a gap in the query sequence (deleted bp)
#  0.0      = % of bases opposite a gap in the repeat consensus (inserted bp)
#  HSU08988 = name of query sequence
#  6563     = starting position of match in query sequence
#  6781     = ending position of match in query sequence
#  (22462)  = no. of bases in query sequence past the ending position of match
#  C        = match is with the Complement of the repeat consensus sequence
#  MER7A    = name of the matching interspersed repeat
#  DNA/MER2_type = the class of the repeat, in this case a DNA transposon 
#            fossil of the MER2 group (see below for list and references)
#  (0)      = no. of bases in (complement of) the repeat consensus sequence 
#             prior to beginning of the match (0 means that the match extended 
#             all the way to the end of the repeat consensus sequence)
#  337      = starting position of match in repeat consensus sequence
#  104      = ending position of match in repeat consensus sequence
#  20       = unique identifier for individual insertions 


#--------------------------------------------------------------
#Format the file of headers corresponding to all TEs
#--------------------------------------------------------------

#The original files leave strand either 'C' or nothing (ie one less column).
#Need to make sure name of matched TE(matchingRepeat column) contains the correct info.
allTEs.loc[allTEs.strand != 'C', 'matchingRepeat'] = allTEs['strand']
allTEs.loc[allTEs.strand != 'C', 'id'] = allTEs['dunno']

allTEs = allTEs.drop(["percDiv","percDel","percIns","qLeft", "strand","rLeft","rPosBeg","rPosEnd","dunno"],1)

#Only keep the specific name of the matched TE ([0]) or TE class/family ([1])
allTEs["matchingRepeat"]=allTEs["matchingRepeat"].str.split('#').str[1]

#add column of Kimura Scores that were extracted using ExtractKimuraScores.sh
pathToKimuraScores = inputdir + "/onlyKimuraScores_filtered_From_Aag2_Contigs.fa.align.txt"
allTEs_KimuraScores =  pd.read_csv(pathToKimuraScores,
                      names = ["text","KimuraScore"],
                      sep = '='
                      )
#add column of Kimura Scores that were extracted using ExtractKimuraScores.sh
#Previous scripts have ensured that number of rows is equal in each.
allTEs['KimuraDivergence'] = allTEs_KimuraScores['KimuraScore']

allTEs["matchingRepeat"]=allTEs["matchingRepeat"].str.split('/').str[0]

allTEs.to_csv(outputdir + "Aag2_dotFAdotAlign_WithKimura.txt",
                   sep = "\t",
                   index = False,
                   quoting = False)

#In previous version, I then merged the .align file (with Kimura scores) with the dot out file. I originally did this because the .out file's TE coordinates better matched
#the coordinates of the .bed file. However, in cases where the .align files coordinates do not perfectly match the .out/.bed's coordinates, the divergence scores are also
#subsequently off because the TE read is likely shorter. May not be far off, but not a fair comparison.
