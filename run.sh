#!/usr/bin/tcsh

foreach i (`seq 1 1 400`)
	qsub -v j=$i -q run64bit HMC_fitting.sh
	sleep 1
end