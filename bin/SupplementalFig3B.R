#Recruitment Plot for FR-HIT outputs for Exiguobacterium 
#RAWIII 
#First created Nov 20, 2013
#Last update Aug 1st, 2018
#Supplemental Figure 3B

#Load Libraries
library("ggplot2")
library("reshape")
library("gridExtra")

#Set working Directory
setwd("~/exi/")

#Load Data
all.data  <- read.delim("PLMB20_ExiPav.txt")
forward.1  <- all.data[all.data[,"Strand"]=="Negative", ]
reverse.2  <- all.data[all.data[,"Strand"]=="Positive", ]

all.data1  <- read.delim("PLMB20_ExispS17.txt")
all.data2  <- read.delim("PLMB20_ExispAT1b.txt")

#Plot in ggplot2
#Option1 - Color (Show both Strands)
svg('~/exi/plot.svg', width = 12, height = 7, pointsize = 16)
p1.leg <- ggplot(all.data, aes(x=Begin.1, y=Identity, colour=Strand)) + geom_point() 
p1  <- ggplot(all.data, aes(x=Begin.1, y=Identity, colour=Strand)) + geom_point() + theme(legend.position="none")
p2  <- ggplot(all.data1, aes(x=Begin.1, y=Identity, colour=Strand)) + geom_point() + theme(legend.position="none")
p3  <- ggplot(all.data2, aes(x=Begin.1, y=Identity, colour=Strand)) + geom_point() + theme(legend.position="none")

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

leg<-g_legend(p1.leg)

grid.arrange(arrangeGrob(p1,p2,p3, leg, ncol=4, widths=c(5/7, 5/7, 5/7, 3/16)))
dev.off()



