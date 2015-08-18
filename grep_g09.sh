#!/bin/bash

Usage(){
  echo "Usage:
    $(basename $0) <file> step
    $(basename $0) <file> conv

Prints Number of optimisation steps or convergence criteria 
in Gaussian outfile"
  exit 0
}

if [ -z "$1" ] || [ -z "$2" ]; then
  Usage
fi

file=$1

if [ "$2" == "step" ]; then
  cat $file | grep "Step number" | tail -n 1
elif [ "$2" == "conv" ]; then
  cat $file | grep "Converged" -A 4
else
  Usage
fi
