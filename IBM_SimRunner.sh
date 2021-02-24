#/mnt/md0/GitHubProjects/IBM-2020/IBM_SimRunner.sh

#!/bin/bash
## Set heritability values for size and growth
sH=0.2
s=${sH: -1}
gH=0.2
g=${gH: -1}
## Set number of dams
nD=5
for simNum in {1..100}
	do
	FILE="g"$g"s"$s"Sim"$simNum"Results"$nD"Dams.csv"
	if [ -f "$FILE" ]; then
		echo "$FILE exists."
	else
		Rscript IBM_header.R $simNum $sH $gH $nD
		d=`date`
		echo "----- Simulation $simNum complete at $d -----"
	fi
	done
exit 0
