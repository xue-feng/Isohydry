#!/bin/bash
# Job name:
#SBATCH --job-name=sa_sampling
#
# Account:
#SBATCH --account=fc_drought
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
## Command(s) to run:
git checkout master
source ~/VirtualEnvirons/sobolenv/bin/activate
x=$(python testing_SAVIO.py)
DATE='date +%Y-%m-%d:%H:%M:%S'
TITLESTR='SOBOL_$DATE'
sendmail xue.feng@berkeley.edu << EOF
subject:$TITLESTR
$x
EOF
