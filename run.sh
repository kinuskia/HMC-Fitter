#!/usr/bin/tcsh

foreach i (`seq 1 1 250`)
	qsub -v j=$i -q run64bit HMC_fitting.sh
end