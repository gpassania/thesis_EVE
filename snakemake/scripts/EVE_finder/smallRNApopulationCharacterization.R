# For a given small RNA population, graphs 2 outputs:
# 1. SeqLogo for + and - strands
# 2. Distribution of sizes
# NOTE: the input files are .map files from bowtie mapping (collapsed by fastx_collapser); but in condensed version. So need to manually 'expand' each piRNA by its counts
library(seqLogo)
library(Biostrings)
library(ggplot2)
library(dplyr)
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##---------------------------FUNCTIONS----------------------------------------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------

#*****Note: genomeSeq in following function means the genome sequence the piRNA mapped to. THIS INCLUDES MISMATCHES, so this is the rev. complement of the piRNA in the antisense case. It is the exact sequence in the sense case.
#piRNAstart is the 5' position IN THE REFERENCE the piRNA maps to. In the case of antisense, this is the position mapping to the 3' end of the piRNA
#piRNAfastaID is a unique ID asigned within the FASTA file (ie before mapping). This allows for comparison of piRNAs found in two different mapping/bowtie outputs

getData = function(pathToFile, filterBySize)
  {
  piRNAmappings = read.table(pathToFile,
                             sep = "\t",
                             header = FALSE,
                             col.names = c("piRNAfastaID","piRNAstrand","Contig","piRNAstart","genomeSeq"),
                             stringsAsFactors = FALSE)
  if (filterBySize == 'yes')
    {
    piRNAmappings = piRNAmappings[nchar(piRNAmappings$genomeSeq) >= 24 & nchar(piRNAmappings$genomeSeq) <= 30,] #Filter by size
    }
  # piRNAmappings$PhasedPair = "NONE" #Will eventually label sense-antisense piRNAs whose 5' ends are offset by 10 nucleotides
  piRNAmappings$Counts = as.numeric(unlist(lapply(strsplit(piRNAmappings$piRNAfastaID,split = "-"),'[[',2))) #Only take second part of piRNAfastaID collumn (everything afte the '-')
  piRNAmappings$piRNAfastaIDonly = as.numeric(unlist(lapply(strsplit(piRNAmappings$piRNAfastaID,split = "-"),'[[',1))) #Only take first part of piRNAfastaID collumn (everything afte the '-')
  
  
  return(piRNAmappings)
  }

#since mapping was done with fastx collapsed small RNA fasta files, generating a position weight matrix seems best way to get info for seqLogo
generatePWM = function(data)
  {
  countsMatrix= matrix(data = 0,
                       nrow = 24, 
                       ncol = 4,
                       dimnames = list(c(1:24), 
                                       c("A","C","G","T")))
  
  #iterate through each row (ie each unique small RNA sequence)
  for (currentRow in 1:nrow(data))
    {
    currentCounts = data$Counts[currentRow]
    strand = data$piRNAstrand[currentRow]
    
    if (strand == "+")
      {
      currentSequence = data$genomeSeq[currentRow]
      }
    
    #if the small RNA maps anti-sense, then need to take reverse complement using Biostrings package
    if (strand == "-")
      {
      currentSequence = data$genomeSeq[currentRow]
      currentSequence = DNAString(currentSequence)
      currentSequence = as.character(reverseComplement(currentSequence))
      }
    
    #Split small RNA sequence so can iterate through one base at a time
    splitSeq = unlist(strsplit(currentSequence, split = ""))
     for (currentPos in 1:24)
      {
      currentPos_sample = splitSeq[currentPos] #current base
      countsMatrix[currentPos, currentPos_sample] = countsMatrix[currentPos, currentPos_sample] + currentCounts #add counts for current position (1-24), and corresponding base
      }
    
    }
  freqMatrix = countsMatrix/sum(data$Counts) #turn into matrix of frequencies instead of counts
  # sum(countsMatrix[1,])
  # sum(freqMatrix[1,])
  # 
  # t(freqMatrix)
  pwm = makePWM(t(freqMatrix)) #seqLogo requires the rows (not columns) to be the bases
  
  return (pwm)
  # return (t(freqMatrix))
}

plotSizeDist = function(data)
  {
  data$piRNAlength = nchar(data$genomeSeq)
  forHist = data %>% 
    group_by (piRNAlength, piRNAstrand) %>% 
    mutate(forPlot = sum(Counts))
  forHist = forHist[!duplicated(forHist[,c("piRNAstrand","piRNAlength","forPlot")]),] #get one entry per group (RNA size, strand)
  forHist[forHist$piRNAstrand == '-', "forPlot"] = -forHist[forHist$piRNAstrand == '-', "forPlot"] #make antisense counts negative
  
  yMaxVal = max(abs(forHist$forPlot))
  #plot with ggPlot
  g = ggplot(forHist)
  g = g + geom_bar(aes(x = piRNAlength, y = forPlot), stat = "identity") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size=16),
          axis.title = element_text(size=20,face="bold")) +
    scale_x_continuous(breaks = seq(15, 30, by = 1)) +
    coord_cartesian(ylim = c(-yMaxVal, yMaxVal)) +
    ggtitle("Size distribution of small RNAs") +
    xlab("small RNA Length") + 
    ylab("Count")
  print (g)
  # return (forHist)
  }

plotSizeDist_3conditions = function(data)
  {
  data$piRNAlength = nchar(data$genomeSeq)

  #Group by small RNAs who share the same length, strand, and condition together, then get total counts for that group
  forHist = data %>% 
    group_by (piRNAlength, piRNAstrand, condition) %>% 
    mutate(forPlot = sum(Counts)/totalLibraryCounts)

  forHist = forHist[!duplicated(forHist[,c("piRNAstrand","piRNAlength","forPlot", "condition")]),] #get one entry per group (RNA size, strand, and conditon)
  forHist[forHist$piRNAstrand == '-', "forPlot"] = -forHist[forHist$piRNAstrand == '-', "forPlot"] #make antisense counts negative
  
  forHist$condition = factor(forHist$condition, levels = c("dsFluc-B_S2", "dsAgo3-B_S6", "dsPiwi4-B_S8")) #set order for graphing by ggplot using factors
  
  yMaxVal = max(abs(forHist$forPlot))
  
  g = ggplot(forHist)
  g = g + geom_bar(aes(x = piRNAlength, y = forPlot, fill = condition), color = "black",position = "dodge", stat = "identity") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size=16),
          axis.title = element_text(size=20,face="bold")) +
    scale_x_continuous(breaks = seq(15, 30, by = 1)) +
    coord_cartesian(ylim = c(-yMaxVal, yMaxVal)) +
    ggtitle("Size distribution of small RNAs") +
    xlab("small RNA Length") + 
    ylab("Fraction of total library")
  print (g)
  # return (forHist)
}

##----------------------------------------------------------------------------------------
##---------------------------Create relevant individual datasets--------------------------
##----------------------------------------------------------------------------------------

outputdir = "/path/to/directory/containing/folders/of/piRNAmapping/samples/" #This should also have the fastx collapsed small RNA fasta files
currentCondition = "dsFluc-B_S2"
currentSample = "PCLV_NC"

#small RNAs mapping to EVEs-------------
piRNAsMappingToEVEs_sizeFiltered = getData(paste(outputdir, "Aag_EVEs","/",currentCondition,"_c_", "Aag_EVEs", "_v1.map",sep = ""), "yes")
piRNAsMappingToEVEs_NOTsizeFiltered = getData(paste(outputdir, "Aag_EVEs","/",currentCondition,"_c_", "Aag_EVEs", "_v1.map",sep = ""), "no")

piRNAsMappingToEVEs_sizeFiltered_posStrand = piRNAsMappingToEVEs_sizeFiltered[piRNAsMappingToEVEs_sizeFiltered$piRNAstrand == "+",]
piRNAsMappingToEVEs_sizeFiltered_negStrand = piRNAsMappingToEVEs_sizeFiltered[piRNAsMappingToEVEs_sizeFiltered$piRNAstrand == "-",]

piRNAsMappingToEVEs_NOTsizeFiltered_posStrand = piRNAsMappingToEVEs_NOTsizeFiltered[piRNAsMappingToEVEs_NOTsizeFiltered$piRNAstrand == "+",]
piRNAsMappingToEVEs_NOTsizeFiltered_negStrand = piRNAsMappingToEVEs_NOTsizeFiltered[piRNAsMappingToEVEs_NOTsizeFiltered$piRNAstrand == "-",]

#small RNAs mapping to Virus-------------
piRNAmappingToVirus1Genome_sizeFiltered = getData(paste(outputdir, currentSample,"/",currentCondition,"_c_", currentSample, "_v3.map",sep = ""), "yes")
piRNAmappingToVirus1Genome_NOTsizeFiltered = getData(paste(outputdir, currentSample,"/",currentCondition,"_c_", currentSample, "_v3.map",sep = ""), "no")

piRNAmappingToVirus1Genome_sizeFiltered_posStrand = piRNAmappingToVirus1Genome_sizeFiltered[piRNAmappingToVirus1Genome_sizeFiltered$piRNAstrand == "+",]
piRNAmappingToVirus1Genome_sizeFiltered_negStrand = piRNAmappingToVirus1Genome_sizeFiltered[piRNAmappingToVirus1Genome_sizeFiltered$piRNAstrand == "-",]

piRNAmappingToVirus1Genome_NOTsizeFiltered_posStrand = piRNAmappingToVirus1Genome_NOTsizeFiltered[piRNAmappingToVirus1Genome_NOTsizeFiltered$piRNAstrand == "+",]
piRNAmappingToVirus1Genome_NOTsizeFiltered_negStrand = piRNAmappingToVirus1Genome_NOTsizeFiltered[piRNAmappingToVirus1Genome_NOTsizeFiltered$piRNAstrand == "-",]

#small RNAs mapping to Genome------------
piRNAsMappingToAag2Genome_uniqueMappers_NOTsizeFiltered = getData(paste(outputdir, "Aag2_PacBio_uniqueMappers","/",currentCondition,"_c_", "Aag2", "_m1v1.map",sep = ""), "no")

##----------------------------------------------------------------------------------------
##---------------------------Get total counts for each library condition------------------
##----------------------------------------------------------------------------------------
currentCondition = "dsAgo3-B_S6"

fullLibrary = read.table(paste(outputdir, currentCondition, "_c.fa",sep = ""))
fullLibrary2 = data.frame(fullLibrary[c(TRUE,FALSE),])
fullLibrary2$Counts = as.numeric(unlist(lapply(strsplit(as.character(fullLibrary2$fullLibrary.c.TRUE..FALSE....),split = "-",fixed = TRUE),'[[',2)))
sum(fullLibrary2$Counts)
##----------------------------------------------------------------------------------------------
##---------------------------Create combined datasets from multiple dsRNA treatmentconditions---
##----------------------------------------------------------------------------------------------

#For plotting small RNA size distributions for multiple conditions
conditions = c("dsFluc-B_S2", "dsAgo3-B_S6","dsPiwi4-B_S8")
correspondingLibraryTotalCounts = c(14008617, 6153485, 10108849) #calculated total library counts above
readAll=vector(mode="list",length=length(conditions)) #each entry in this dataset will be an entire small RNA mapping dataset corresponding to one condition in 'conditions'

for (currentCondition in 1:length(conditions))
  {
  #Load in one dataset----------------------------------
  #virus
  # data = getData(paste(outputdir, currentSample,"/", conditions[currentCondition], "_c_", currentSample, "_v3.map",sep = ""), "no")
  #OR-------
  #EVEs
  data = getData(paste(outputdir, "Aag_EVEs","/",conditions[currentCondition],"_c_", "Aag_EVEs", "_v1.map",sep = ""), "no")
  #OR-------
  #Aag2 Genome (unique mappers)
  # data = getData(paste(outputdir, "Aag2_PacBio_uniqueMappers","/",conditions[currentCondition],"_c_", "Aag2", "_m1v1.map",sep = ""), "no")
  #-----------------------------------------------------
  
  #Specify column giving current condition, for grouping later
  data$condition = conditions[currentCondition]
  
  #Since this is comparing distributions of 3 different libraries, will need to normalize counts to fraction of total library
  data$totalLibraryCounts = correspondingLibraryTotalCounts[currentCondition]
  
  readAll[[currentCondition]] = data
  }
combined_3conditions = do.call(rbind, readAll) #bind all datasets together from all conditions into one large data set
plotSizeDist_3conditions(combined_3conditions)
##----------------------------------------------------------------------------------------
##---------------------------Create Seq Logos (first 24 bp)-------------------------------
##----------------------------------------------------------------------------------------

#First, make a Position Weight Matrix for sample of interest
forLogo = generatePWM(piRNAsMappingToEVEs_sizeFiltered_posStrand)
forLogo = generatePWM(piRNAsMappingToEVEs_sizeFiltered_negStrand)

forLogo = generatePWM(piRNAmappingToVirus1Genome_sizeFiltered_posStrand)
forLogo = generatePWM(piRNAmappingToVirus1Genome_sizeFiltered_negStrand)

#A limitation of this package is that it doens't allow RNA logos to be made, but it does accept
#position weight matricies as input, which is best given the small RNA mapping files I'm using are collapsed (ie include counts and do not have the same sequences repeated over and over).
seqLogo(forLogo)
##----------------------------------------------------------------------------------------
##---------------------------piRNA size distributions-------------------------------------
##----------------------------------------------------------------------------------------

plotSizeDist(piRNAsMappingToAag2Genome_uniqueMappers_NOTsizeFiltered)

plotSizeDist(piRNAsMappingToEVEs_NOTsizeFiltered)

plotSizeDist(piRNAmappingToVirus1Genome_NOTsizeFiltered)
