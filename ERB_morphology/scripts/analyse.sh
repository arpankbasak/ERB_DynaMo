#!/bin/bash

# conda activate erb_dynamo;

echo "Segregation initiated ..." ;

./script/0_segregate.R ./images/raw_images/ ./data/global_parameters.txt;

echo "DONE!" ;

echo "Segmentation of ERBs ..." ;

./script/1_segmentation.R ./data/global_parameters.txt;

echo "DONE!" ;

echo "Feature matrix build ..." ;

./script/2_featurematrix.R;

echo "DONE!" ;

echo "Computing statistics for morphological parameters :: Image-wise analysis ..." ;

./script/4_1_Image.R;

echo "Computing statistics for morphological parameters :: Segmented Cell-wise analysis ..." ;

./script/4_2_SegmentedCells.R;

echo "Pooling clustered features ..." ;

./script/5_Pool_Features.R;

echo "Building video for time series projection ..." ;

./script/6_Clustering_features.R;

echo "DONE!" ;

echo "JOB COMPLETED" ;