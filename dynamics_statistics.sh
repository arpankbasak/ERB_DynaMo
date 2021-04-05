#!/bin/bash

# conda activate erb_dynamo;
WKDIR=$0
PATH_TO_WKDIR=$1
PATH_TO_metadata=$2

cd $WKDIR/ERB_dynamics/;

echo "Calculation and statistical comparison of displacement parameter ..." ;

./script/4_moment_calculation.R $PATH_TO_WKDIR $PATH_TO_metadata;

echo "JOB COMPLETED" ;