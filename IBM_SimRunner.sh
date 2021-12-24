#/mnt/md0/GitHubProjects/IBM-2020/IBM_SimRunner.sh

#!/bin/bash
## Set heritability values for size and growth
## for sH in 0.1 0.2 0.3 0.4 0.5
for sH in 0.2
	do
	s=${sH: -1}
##	for gH in 0.1 0.2 0.3 0.4 0.5
  for gH in 0.2
		do
		g=${gH: -1}
		for nD in {0..5}
			do
			for simNum in {1..100}
				do
				FILE="g"$g"s"$s"Sim"$simNum"Results"$nD"Dams.csv"
				if [ -f "$FILE" ]; then
					echo "$FILE exists."
				else
					Rscript IBM_header.R $simNum $sH $gH $nD
					d=`date`
					echo "----- Simulation $simNum with $nD dams complete at $d -----"
				fi
				done
			done
		done
	done
exit 0
