#!/bin/bash

# conda activate erb_dynamo;

echo "Segregation initiated ..." ;

WKDIR=$0
PATH_TO_IMAGE=$1
PATH_TO_PARAM=$2
cd $WKDIR/ERB_morphology/;


./script/0_segregate.R $PATH_TO_IMAGE $PATH_TO_PARAM;

echo "DONE!" ;

echo "Segmentation of ERBs ..." ;

./script/1_segmentation.R $PATH_TO_PARAM;

echo "DONE!" ;

echo "Feature matrix build ..." ;

./script/2_featurematrix.R;

echo "DONE!" ;

echo "JOB COMPLETED" ;