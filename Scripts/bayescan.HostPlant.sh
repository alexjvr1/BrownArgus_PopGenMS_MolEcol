#!/bin/bash
#PBS -N bayescan1
#PBS -l nodes=1:ppn=1 #1nodes, 1 processor per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##20 hours wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
##PBS -t 1-10

#run job in working directory
cd $PBS_O_WORKDIR  #run job in working directory
#cd 1a_Aricia_agestis_PopGenomics/Bayescan

# set a local temporary folder.  Useful if you create a lot of temp files and do not want to receive an email saying that you are filling up all the memory of the common temporary folder. It happened to me!
export TMPDIR=$HOME/.local

echo "START ----------------------------"

#load your program if it is installed globally or the modules you used to install your program locally (compilers, etc) 
#Check what is available with module avail
#module load languages/R-3.0.2

# Run program

~/software/BayeScan2.1/binaries/BayeScan2.1_linux64bits AA251.HostPlant.bayescan -output AA251.HostPlant.out -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 20


move the output on server to 
/newhome/aj18951/1a_Aricia_agestis_PopGenomics/Bayescan/Out_HostPlant
