#!/bin/bash

# conda activate erb_dynamo;

echo "Segregation initiated ..." ;

WKDIR=$0
PATH_TO_IMAGE=$1
PATH_TO_PARAM=$2

cd $WKDIR/ERB_dynamics/;

# ./script/0_segregate.R $PATH_TO_IMAGE $PATH_TO_PARAM;

./script/0_segregate.R ./images/raw_images/ ./data/global_parameters.txt;
echo "DONE!" ;

echo "Segmentation of ERBs ..." ;

# ./script/1_segmentation_Cells_ERBs.R $PATH_TO_PARAM;
./script/1_segmentation_Cells_ERBs.R ./data/global_parameters.txt;

echo "DONE!" ;

echo "Feature matrix build ..." ;

./script/2_featurematrix.R;

echo "DONE!" ;

echo "Building video for time series projection ..." ;

./script/3_MomentProjection.R;

echo "DONE!" ;

echo "JOB COMPLETED" ;