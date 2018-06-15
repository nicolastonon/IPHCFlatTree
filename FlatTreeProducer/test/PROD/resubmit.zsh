#!/bin/env zsh

tlist=($(ls . | grep "crab_"))
for i in ${tlist}
do
  isfail=$(crab status -d ${i} | grep failed)
  if [[ ${isfail} != "" ]]; then
    echo "${i}: RESUBMIT"
    res=$(crab resubmit -d ${i})
  else
    echo "${i}: DONE"
  fi
done
