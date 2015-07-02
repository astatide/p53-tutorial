#!/bin/bash
#PBS -N P53.TUTORIAL
#PBS -S /bin/bash
#PBS -j oe
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=25
#PBS -q dist_amd
#PBS -m ae
#PBS -A westpa2015

set -x
cd $PBS_O_WORKDIR
source env.sh || exit 1

cd $WEST_SIM_ROOT

# start server in processes mode.
$WEST_ROOT/bin/w_run --work-manager=processes --n-workers=$PBS_NUM_PPN &> west-$PBS_JOBID.log &

# wait on host info file up to one minute
#for ((n=0; n<60; n++)); do
#    if [ -e $SERVER_INFO ] ; then
#        echo "== server info file $SERVER_INFO =="
#        cat $SERVER_INFO
#        break
#    fi
#    sleep 1
#done

# exit if host info file doesn't appear in one minute
#if ! [ -e $SERVER_INFO ] ; then
#    echo 'server failed to start'
#    exit 1
#fi

# start clients, with the proper number of cores on each
#pbsdsh -v -u $PWD/node.sh --work-manager=zmq --zmq-mode=client --n-workers=32 --zmq-info=$SERVER_INFO &

wait
