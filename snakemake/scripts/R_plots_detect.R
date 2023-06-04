#!/usr/bin/env Rscript
library(ggplot2)

#read file as argument
args = commandArgs(trailingOnly=TRUE)

###individual plots
#set name
name = args[1]
name = substr(name,1,nchar(name)-4)
plotname = paste(name, ".png", sep="")

#import file as table
df = read.table(args[1], header=TRUE, sep="\t", row.names=NULL)
colnames(df) = colnames(df)[2:ncol(df)]
df = df[ , - ncol(df)]

#check if any results remain
filesize = file.info(args[1])$size
if (filesize < 100){
	print("no results remaining after filtering")
	file.create(plotname)
	quit(save="no")
}

filesize = file.info(args[1])$size
#plot table
plot = ggplot(df, aes(VirusName, Length, colour = Percentage_ID)) 
plot = plot + geom_jitter(width= 0.1, size = 3.5)
plot = plot + scale_x_discrete(guide=guide_axis(n.dodge=3))

#write to png
png(filename= plotname, width=760,  height=480, units="px", bg = "black")
print(plot)
dev.off()



###combined plot
#set name
name = args[2]
name = substr(name,1,nchar(name)-4)
plotname = paste(name, ".png", sep="")

#import file as table
df = read.table(args[2], header=TRUE, sep="\t", row.names=NULL)
#plot table
plot = ggplot(df, aes(VirusName, Length, colour = SampleID)) +
geom_jitter(width= 0.3, size = 3.5, aes(shape=AssemblyLevel)) +
scale_x_discrete(guide=guide_axis(n.dodge=4))

#write to png
png(filename= plotname, width=1040,  height=760, units="px", bg = "black")
print(plot)
dev.off()
