#!/bin/bash

# conda activate erb_dynamo;

echo "Segregation initiated ..." ;

WKDIR=$0
PATH_TO_IMAT=$1
PATH_TO_metadata=$2

cd $WKDIR/ERB_morphology/;

echo "Image-wise analysis of morphological features ..." ;

./script/4_1_Image.R $PATH_TO_IMAT $PATH_TO_metadata;

echo "Analysing morhology of the sub-cellular compartments among the segmented cells ..." ;

./script/4_2_SegmentedCells.R $PATH_TO_IMAT $PATH_TO_metadata;

echo "Pooling the segmented cells that are distinct from the background ..." ;

./script/5_Pool_Features.R $PATH_TO_IMAT;

echo "Analysising the sub-cellular features of the clustered cells that are distinct from background ..." ;

./script/6_Clustering_features.R $PATH_TO_IMAT $PATH_TO_metadata;

echo "DONE!" ;

echo "JOB COMPLETED" ;