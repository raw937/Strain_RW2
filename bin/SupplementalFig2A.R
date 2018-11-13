#Exiguobacterium Comperative MetaCyc Venn Diagram (High Genome similarity)
# Written by RAWIII March 24, 2014
#First created Nov 20, 2013
#Last update Aug 1st, 2018
#Supplemental Figure 2A

# Load Libraries
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")

# Set working directory
setwd('~/exi/')

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
S17 <- read.delim("e_sp_S17_pgdb.txt", row.names=1)
S8111 <- read.delim("e_sp_8111_pgdb.txt", row.names=1)
AT1B <- read.delim("e_sp_AT1b_pgdb.txt", row.names=1)

# check dimensions of loaded tables match with number of lines - 1
dim(RW2)
dim(S17)
dim(S8111)
dim(AT1B)

# a simpler header
header <- c("long_name", "num_rxn", "num_covered", "num_orfs")

colnames(RW2) <- header
colnames(S17) <- header
colnames(S8111) <- header
colnames(AT1B) <- header

# trim from the top pathways that do not have 1 ORFs
RW2 <- droplevels(RW2[RW2["num_orfs"] >1,])
S17 <- droplevels(S17[S17["num_orfs"] >1,])
S8111 <- droplevels(S8111[S8111["num_orfs"] >1,])
AT1B <- droplevels(AT1B[AT1B["num_orfs"] >1,])

normalize_length = TRUE # normalize by pathway length
normalize_orfs = TRUE # normalize to number of ORFs in sample

# total ORF predicted for each sample
num_orfs_RW2 = 1381
num_orfs_S17 = 1270
num_orfs_S8111 = 1267
num_orfs_AT1B = 1299

# normalize all pathway counts to length
if (normalize_length == TRUE) {
  RW2[,"num_orfs"] <- RW2[,"num_orfs"] / RW2[,"num_rxn"]
  S17[,"num_orfs"] <- S17[,"num_orfs"] / S17[,"num_rxn"]
  S8111[,"num_orfs"] <- S8111[,"num_orfs"] / S8111[,"num_rxn"]
  AT1B[,"num_orfs"] <- AT1B[,"num_orfs"] / AT1B[,"num_rxn"]
}

if (normalize_orfs == TRUE) {
  RW2[,"num_orfs"] <- (RW2[,"num_orfs"] / num_orfs_RW2) * 100
  S17[,"num_orfs"] <- (S17[,"num_orfs"] / num_orfs_S17) * 100
  S8111[,"num_orfs"] <- (S8111[,"num_orfs"] / num_orfs_S8111) * 100
  AT1B[,"num_orfs"] <- (AT1B[,"num_orfs"] / num_orfs_AT1B) * 100
}

#Venn diagram plot
Exi_venn <- venn_diagram4(rownames(RW2), 
                               rownames(S17), 
                               rownames(S8111), 
                               rownames(AT1B),
                               "RW2",
                               "S17",
                               "8111",
                               "AT1B",
                               name_output="default"
)
