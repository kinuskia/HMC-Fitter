#!/usr/bin/tcsh

foreach i (`seq 1 1 52`)
	qsub -v j=$i -q run64bit script.sh
end