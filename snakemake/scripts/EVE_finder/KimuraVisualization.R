
library(ggplot2)

#-----------------------------------------------------------------------------------
#Load .fa.align file (which I believe is what original RepeatMasker Kimura plot is based on)
#-----------------------------------------------------------------------------------
allTEs = read.table("/path/to/save/Aag2_dotFAdotAlign_WithKimura.txt",
                    header = TRUE,
                    na.strings = '',
                    sep = "\t",
                    stringsAsFactors = FALSE)
allTEs['percentGenome'] = ((allTEs$qPosEnd-allTEs$qPosBeg)/1723952533)*100
allTEs = allTEs[allTEs$matchingRepeat!= "Simple_repeat",]
head(allTEs,100)


#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-------------------------------------Scatterplot,Boxplots, and Violin Plots of Kimura Divergence-----------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#Histogram of kimura scores for all TEs. Very similar to default Repeat Masker output
#I need to specify the breaks, rather than just binwidth.
#If I did binwidth=1, it put the first bin from -0.5-0.5, causing many counts to be lost. Not sure why it defaults to that.

# allTEs = allTEs[allTEs$TEclass == "LTR",]
g = ggplot(allTEs, aes(x = as.numeric(as.character(KimuraDivergence))))
g = g + geom_histogram(aes(fill = factor(matchingRepeat, c("DNA","LINE","LTR","Unknown","MITEs","UD","Helitrons","SINE","Penelope","RC","Satellite","Simple_repeat")),
                           weight = percentGenome), color = "black",  breaks=seq(0,50, by=1), position = 'stack') +
# g = g + geom_histogram(aes(fill = matchingRepeat,
#                              weight = percentGenome), color = "black",  breaks=seq(0,50, by=1), position = 'stack') +
  scale_fill_brewer(palette = "Paired") +
  xlab("Kimura Divergence Score") + 
  ylab("Percent Genome") + 
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"))
g
ggsave(filename = "Aag2KimuraHistogram_LTRonly_Aag2_allTEs_KimuraWithPercDivWithDotOutNoDup.pdf",
       path = "/path/to/save/output/",
       dpi = 600)