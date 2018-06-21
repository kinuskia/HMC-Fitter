#!/usr/bin/tcsh
#PBS -l walltime=15:00:00
#PBS -l mem=128mb
cd /lpt/xrxjktfa/HMC-Fitter
./fitter $j
mv /lpt/xrxjktfa/HMC-Fitter/data$j.txt /lpt/xrxjktfa/HMC-Fitter/LPT-Cluster/data$j.txt
