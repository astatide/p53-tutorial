#!/bin/bash -l

set -x

cd $PBS_O_WORKDIR
(source env.sh
cd $WEST_SIM_ROOT

$WEST_ROOT/bin/w_run "$@" ) &> west-$PBS_NODENUM-node.log
