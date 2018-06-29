#!/usr/bin/tcsh

foreach i (`seq 1 1 500`)
	qsub -v j=$i -q run64bit HMC_fitting.sh
	sleep 1
end