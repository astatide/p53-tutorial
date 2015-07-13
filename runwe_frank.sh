#!/bin/bash
#PBS -N P53.TUTORIAL
#PBS -S /bin/bash
#PBS -j oe
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=64
#PBS -q dist_amd
#PBS -m ae
#PBS -A westpa2015

set -x
cd $PBS_O_WORKDIR
source env.sh || exit 1

cd $WEST_SIM_ROOT

$WEST_ROOT/bin/w_run --work-manager=processes --n-workers=$PBS_NUM_PPN &> west-$PBS_JOBID.log &

wait
