 #!/usr/bin/env Rscript

"""
author: Gian Passania

script for visualization of identification pipeline outputs

usage:
Rscript R_plot_ident.R [overlapping_regions.txt] [unique_regions.txt]
"""

#import statements
library(ggplot2)

#read file as argument
args = commandArgs(trailingOnly=TRUE)


###overlapping regions plot
#set name
name = args[1]
name = substr(name,1,nchar(name)-4)
plotname = paste(name, ".png", sep="")

#import file as table
df = read.table(args[1], header=TRUE, sep="\t", row.names=NULL)

#check if any results remain
filesize = file.info(args[1])$size
if (filesize < 100){
	print("no results remaining after filtering")
	file.create(plotname)
	quit(save="no")
}

#plot table
plot = ggplot(df, aes(Chr, Length, colour = fullVirusName)) 
plot = plot + geom_jitter(width= 0.1, size = 3.5)
plot = plot + scale_x_discrete(guide=guide_axis(n.dodge=4))

#write to png
png(filename= plotname, width=1040,  height=760, units="px", bg = "black")
print(plot)
dev.off()


###Unique regions plot
#set name
name = args[2]
name = substr(name,1,nchar(name)-4)
plotname = paste(name, ".png", sep="")
plotname
#import file as table
df = read.table(args[2], header=TRUE, sep="\t", row.names=NULL, fill = TRUE)
df
#plot table
plot = ggplot(df, aes(Chr, colour = fullVirusName)) 
plot = plot + geom_bar(stat="identity")

#write to png
png(filename= plotname, width=1040,  height=760, units="px", bg = "black")
print(plot)
dev.off()
