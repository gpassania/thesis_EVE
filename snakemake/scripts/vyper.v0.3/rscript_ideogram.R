#! /usr/bin/env Rscript

## R script to generate ideogram for Vy-Per results
## version 1.4
##
## Silke Szymczak, IKMB, Kiel, Germany
## February 2015
##
## Additional modification:
## Annotation file contains symbol in column pch
## 
## Modification:
## Add possibility to have plot on a black background
## Default color: RColorBrewer palette Dark2 or rainbow
## Larger plotting symbols
## Different plotting symbols

## Usage:
## rscript_ideogram.R <results file> <plot file> <background> <legend> <main> <bin size> <chimera threshold> [<col file>]
##
## results file: name of Vy-Per cluster results file
## plot file: name of file for pdf plot 
## background: should background be black (TRUE or FALSE)
## legend: should legend be shown on figure? (TRUE or FALSE)
## main: title of plot
## bin size: size of bins in bp
## chimera threshold: threshold for number of candidates
## anno file (optional): name of tab separated annotation file specifying colors, plotting symbols and names to use for different viruses (columns virus: virus name used in Vy-Per result file, color: R color, pch: plotting symbol s number, name: virus name or abbreviation that should be used in legend)

## check arguments
args = commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: rscript_ideogram.R <results file> <plot file> <background> <legend> <main> <bin size> <chimera threshold> [<anno file>] ...")
}


## --------------------------------------------------
## parameters
## --------------------------------------------------
res.file = args[1]
plot.file = args[2]
bg.black = as.logical(args[3])
legend = as.logical(args[4])
main = args[5]
bin.size = as.numeric(args[6])
threshold = as.numeric(args[7])
if (length(args) == 8) {
  anno.file = args[8]
} else {
  anno.file = NULL
}

## --------------------------------------------------
## check arguments
## --------------------------------------------------
if (!file.exists(res.file)) {
  stop(paste("Vy-Per results file", res.file, "does not exist!"))
}

if (is.na(bg.black)) {
  stop("background parameter needs to be TRUE or FALSE")
}

if (is.na(legend)) {
  stop("legend parameter needs to be TRUE or FALSE")
}

if (!is.null(anno.file) && ! file.exists(anno.file)) {
  stop(paste("Annotation file", anno.file, "does not exist!"))
}

## --------------------------------------------------
## functions
## --------------------------------------------------
library(quantsmooth)
library(RColorBrewer)

plot.ideogram <- function(loc.hits, col.virus, pch.virus, bleach = 0, main = "") {
  
  plot.points <- function(chrompos, ind = NULL, offset, col, pch = 19) {
    if (is.null(ind)) ind = 1:nrow(chrompos)
    points(chrompos[ind, 2], chrompos[ind, 1] + offset,
           pch = pch, cex = 1.25, col = col)
  }
  
  ## plot ideograms
  chrompos = data.frame(CHR = loc.hits$chr,
                        MapInfo = loc.hits$pos, stringsAsFactors = FALSE)
  
  chrompos = prepareGenomePlot(chrompos = chrompos,
                               paintCytobands = TRUE,
                               bleach = bleach,
                               organism = "hsa",
                               units = "bases",
                               sexChromosomes = TRUE,
                               topspace = 0,
                               main = main,
                               cex.main = 2)
  
  ## subheading with parameter settings
  if (bin.size >= 10^6) {
    bin.size = bin.size/(10^6)
    unit = "Mb"
  } else if (bin.size >= 10^3) {
    bin.size = bin.size/(10^3)
    unit = "kb"
  } else {
    unit = "bp"
  }
  text = paste("chimera threshold:", threshold, "paired-end(s)", "    ",
               "bin size:", bin.size, unit)
  #   mtext(text, side = 3, line = -3.5, outer = TRUE)
  mtext(text, side = 1, line = -0.5)
  
  ## plot hits virus 1
  plot.points(chrompos = chrompos, offset = 0.2, col = col.virus[loc.hits$virus1],
              pch = pch.virus[loc.hits$virus1])
  
  ## plot additional virus hits
  for (i in 2:3) {
    name = paste("virus", i, sep = "")
    ind = which(!is.na(loc.hits[, name]))
    if (length(ind) > 0) {
      plot.points(chrompos = chrompos, ind = ind, offset = 0.2*i,
                  col = col.virus[loc.hits[ind, name]],
                  pch = pch.virus[loc.hits[ind, name]])
    }
  }
}


## --------------------------------------------------
## prepare data of cluster hits
## --------------------------------------------------
loc.hits = read.delim(file = res.file, header = TRUE, as.is = TRUE)
## remove chr part of chromosomes
loc.hits$chr = substr(loc.hits$Chr, 4, 99)
## calculate middle point of cluster
loc.hits$pos = loc.hits$cluster_start + (loc.hits$cluster_end - loc.hits$cluster_start)/2
## select hits on autosomes and X, Y
loc.hits = subset(loc.hits, chr %in% c(1:22, "X", "Y"))

if (nrow(loc.hits) == 0) 
  stop("No cluster hits found!")
print(paste("Loaded", nrow(loc.hits), "cluster hits on chromosomes 1-22, X, Y"))


## --------------------------------------------------
## colors
## --------------------------------------------------
## extract unique set of viruses
virus = sort(unique(na.omit(unlist(lapply(grep("virus", colnames(loc.hits)), function(x) {
  loc.hits[, x]})))))

if (length(virus) == 0)
  stop("No virus found!")
print(paste(length(virus), "viruses found"))

if (!is.null(anno.file)) {
  info = read.table(file = anno.file, header = TRUE, as.is = TRUE, sep = "\t")
  if (!all(virus %in% info$virus)) {
    stop(paste("Some viruses not found in annotation file:",
               paste(setdiff(virus, info$virus), collapse = ", ")))
  }
  col.virus = info$color
  names(col.virus) = info$virus
  col.virus = col.virus[virus]
  
  pch.virus = info$pch
  names(pch.virus) = info$virus
  pch.virus = pch.virus[virus]
  
  names.virus = info$name
  names(names.virus) = info$virus
  names.virus = names.virus[virus]
  
} else {
  print("using internal colors und symbols")
  
  #    col.virus = c("black", "cyan", "magenta", "gray", "orange", "green", "purple", "blue",
  #        "darkgreen", "gold2")
  #    col.virus = c(col.virus, rainbow(length(virus) - length(col.virus)))[1:length(virus)]
  if (length(virus) < 3) {
    col.virus = c("red", "blue")[1:length(virus)]
  } else if (length(virus) <= 8 && length(virus) >= 3) {
    col.virus = brewer.pal(length(virus), "Dark2")
  } else {
    col.virus = rainbow(length(virus))
  }
  names(col.virus) = virus
  names.virus = virus
  
  pch.virus = rep(15:25, 5)[1:length(virus)]
  names(pch.virus) = names(col.virus)
}


## --------------------------------------------------
## plot
## --------------------------------------------------
if (legend) {
  no.lines = ceiling(length(virus) / 2)
  layout.m = matrix(1:2, ncol = 1)
  prop = 0.05 * no.lines
  layout.h = c(1 - prop, prop)
} else {
  layout.m = matrix(1, nrow = 1)
  layout.h = c(1)
  
}
pdf(file = plot.file, width = 9, height = 9, useDingbats = FALSE)
if (bg.black) {
  op = par(fg = "white", bg = "black",
           col = "white", col.axis = "white",
           col.lab = "white", col.main = "white",
           col.sub = "white")
  if (col.virus[1] == "black") col.virus[1] = "white"
  bleach = 0.5
} else {
  op = par()
  bleach = 0
}
layout(layout.m, height = layout.h)
plot.ideogram(loc.hits = loc.hits, col.virus = col.virus, pch.virus = pch.virus, 
              bleach = bleach, main = main)

## legend
if (legend) {
  old.par = par(mar=c(0.1, 0.5, 0, 0))
  plot.new()
  plot.window(c(-1, 1), c(0.5, no.lines + 0.5))
  
  for (i in 1:no.lines) {
    end = 2 * i
    if (i == no.lines) end = min(end, length(virus))
    ind = (2*i-1):end
    v = names.virus[ind]
    temp = legend(0, no.lines - i + 1, v, horiz = TRUE,
                  col = col.virus[ind], pch = pch.virus[ind], bty = "n", xjust = 0.5, yjust = 0.5,
                  text.width = 1)
  }
  par(old.par)
}
suppressWarnings(par(op))
dev.off()

print("plotting finished")
