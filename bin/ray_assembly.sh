#!/bin/bash

#Ray assembly of strain RW2
#Written by RAWIII 
#First created Feb 1st, 2013
#Last update Aug 1st, 2018

#Install ray assembly dependences on ubuntu 12.04
sudo apt-get remove libcr-dev mpich2 mpich2-doccd 
sudo apt-get install libopenmpi-dev openmpi-bin openmpi-doc

#Install ray
#change the max k-mer size in the config file from 34 to 150.
sudo cp ray-build/* /usr/local/bin/
make PREFIX=/usr/local/bin/ray-build

#Run ray assembler
mpiexec -n 16 Ray -k 55 -i exi_RR_SS.fastq -o exi55k/
