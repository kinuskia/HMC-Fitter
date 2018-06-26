#!/usr/bin/tcsh
#PBS -l walltime=05:00:00
#PBS -l mem=64mb
cd /lpt/xrxjktfa/HMC-Fitter
./fitter $j
mv /lpt/xrxjktfa/HMC-Fitter/data$j.txt /lpt/xrxjktfa/HMC-Fitter/LPT-Cluster/data$j.txt
