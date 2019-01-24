#!/bin/bash
export mcfost=${HOME}'/../cpinte/mcfost_cigri/mcfost'
export MCFOST_UTILS=${HOME}'/mcfost_utils/'
export OMP_NUM_THREADS=8
    
mkdir "$OAR_JOBID"

for param in *.par;do
if [ -f $param ]; then
      new_dir="`ls "$param" | sed s/.par//`"
      mv "$new_dir".par "$OAR_JOBID"
      cd `echo $OAR_JOBID`
      mkdir "$new_dir"
      $mcfost "$new_dir".par -rt > sed.log
      mv sed.log data_th
      mv data_* "$new_dir"
      cd ..
    fi
done

