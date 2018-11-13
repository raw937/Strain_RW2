#Exiguobacterium Comperative MetaCyc Venn Diagram (Low Genome similarity)
#Written by RAWIII 
#First created March 24, 2014
#Last update Aug 1st, 2018
#Supplemental Figure 2B

# Load Libraries
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")

# Set working directory
setwd('/home/rawiii/Desktop/TestDrive/')

# Source for R program operation
source("venn_diagram4wSvg.r")

# Load MetaCyc v17 database
meta17 <- read.table("meta_17.txt", sep="\t", header=T, row.names=1)

# Load MetaCyc v17 Hierarchy and set levels
meta_17_hier <- read.table("meta_17_hierarchy.txt", sep="\t", row.names=1)
meta_17_hier$V3 <- factor(meta_17_hier$V3, levels=union(meta_17_hier$V3,meta_17_hier$V3))
meta_17_hier$V4 <- factor(meta_17_hier$V4, levels=union(meta_17_hier$V4,meta_17_hier$V4))

# Load Datasets (Delim)
# read in the datasets (tab delimited)
RW2 <- read.delim("e_pavilionensis_pgdb.txt", row.names=1)
B7 <- read.delim("e_antarcticum_B7_pgdb.txt", row.names=1)
MH3 <- read.delim("e_sp_mh3_pgdb.txt", row.names=1)
S255 <- read.delim("e_sibiricum_pgdb.txt", row.names=1)

# check dimensions of loaded tables match with number of lines - 1
dim(RW2)
dim(B7)
dim(MH3)
dim(S255)

# a simpler header
header <- c("long_name", "num_rxn", "num_covered", "num_orfs")

colnames(RW2) <- header
colnames(B7) <- header
colnames(MH3) <- header
colnames(S255) <- header

# trim from the top pathways that do not have 1 ORFs
RW2 <- droplevels(RW2[RW2["num_orfs"] >1,])
B7 <- droplevels(B7[B7["num_orfs"] >1,])
MH3 <- droplevels(MH3[MH3["num_orfs"] >1,])
S255 <- droplevels(S255[S255["num_orfs"] >1,])

normalize_length = TRUE # normalize by pathway length
normalize_orfs = TRUE # normalize to number of ORFs in sample

# total ORF predicted for each sample
num_orfs_RW2 = 1381
num_orfs_B7 = 1334
num_orfs_MH3 = 1417
num_orfs_S255 = 1396

# normalize all pathway counts to length
if (normalize_length == TRUE) {
  RW2[,"num_orfs"] <- RW2[,"num_orfs"] / RW2[,"num_rxn"]
  B7[,"num_orfs"] <- B7[,"num_orfs"] / B7[,"num_rxn"]
  MH3[,"num_orfs"] <- MH3[,"num_orfs"] / MH3[,"num_rxn"]
  S255[,"num_orfs"] <- S255[,"num_orfs"] / S255[,"num_rxn"]
}

if (normalize_orfs == TRUE) {
  RW2[,"num_orfs"] <- (RW2[,"num_orfs"] / num_orfs_RW2) * 100
  B7[,"num_orfs"] <- (B7[,"num_orfs"] / num_orfs_B7) * 100
  MH3[,"num_orfs"] <- (MH3[,"num_orfs"] / num_orfs_MH3) * 100
  S255[,"num_orfs"] <- (S255[,"num_orfs"] / num_orfs_S255) * 100
}

#Venn diagram plot
Exi_venn <- venn_diagram4(rownames(RW2), 
                               rownames(B7), 
                               rownames(MH3), 
                               rownames(S255),
                               "RW2",
                               "B7",
                               "MH3",
                               "255-15",
                               name_output="default"
)
