#!/bin/bash

#Written RAWIII
#First created Aug 1, 2018

##Align with mafft
~/mafft-7.394-without-extensions/core/./mafft --localpair --maxiterate 1000 fig2_exi_16S_sequences.fasta >fig2_exi_16S_alignment.clust

##Build phylogeny with iqtree
~/iqtree-1.6.2-Linux/bin/iqtree -s fig2_exi_16S_alignment.clust -st DNA -m TEST -bb 1000 -alrt 1000

## Yields contree file unannotated. Rest annotated by hand
