# Read Me

## ERB-DYNAMO (DYNamics And MOrphology)
Analysis of confocal micrographs using unbiased approach.

---
## Summary

### Title
Texture feature extraction from microscope images enables robust estimation of ER body morphology in Arabidopsis

### Background
Cellular components are controlled by genetic or physiological factors that define their shape and size. However, quantitively capturing the morphological characteristics of cellular organelles from the micrograph images is challenging, because the analysis deals with complexities of images that frequently lead inaccuracy in the estimation of the features. Here we show a unique quantitative method to overcome biases and inaccuracy of biological samples from confocal micrographs. 

### Results
We generated 2D images of cell walls and spindle-shaped cellular organelles, namely ER bodies, with a maximum contrast projection of 3D confocal fluorescent microscope images. The projected images are further processed and segmented by adaptive thresholding of fluorescent levels in the cell walls. Micrographs are composed of pixels, which have the information of position and intensity. From the pixel information, we calculated three types of features (spatial, intensity and Haralick) in each ER bodies corresponding to segmented cells. The spatial features include basic information of morphology, e.g., surface area and perimeters. The intensity features include information of mean, standard deviation and quantile of intensities within an ER body. Haralick features describe the texture features, which can be calculated mathematically from the interrelationship between the pixel information. Together these parameters are subjected to multivariate analysis to estimate the morphological diversity. We captured similar morphological diversity within ER body phenotype on several microscopy experiments performed in different settings and scanned under different objectives. Consequently, we describe morphological differences of ER bodies between A. thaliana wild type and the mutants deficient in ER body-related genes based on the morphological parameters and determine the morphological variations within the cellular features. 

### Conclusion
The findings unexpectedly revealed multiple genetic factors that are involved in the shape and size of ER bodies in A. thaliana. This is the first report showing morphological characteristics of cellular components are quantitatively measurable to distinguish plant phenotypes even in plants that shows similar cellular components. The estimation of morphological diversity is independent of cell staining method and objective lends used in microscopy. Hence, our study enables the robust estimation of plant phenotypes by recognizing small differences of complex cell organelle shapes, which is further beneficial in a comprehensive analysis of the molecular mechanism for cell organelle formation independent of technical variations. 

### Keywords
_ER-Body, confocal microscopy, morphology, Haralick feature, quantitative analysis_

---

__Applications__
* Morphological comparison between control and treatment or wildtype and mutants
* Counting objects lying within the cellular space
* Compute statistics from micrograph
* Clustering of cells having similar features
* Multivariate analysis
* Dynamic behaviour of the features in red/green channel along time

---

> *"One confocal image, ~100 cells, ~10 to ~15 sub-cellular features within the cells. ERB_dynamo computes the dynamics and morphology of the sub-cellular features in the hierarchy data-structure"*

### Lets see what `ERB_dynamo` can do

* This is an assembly of tools used to quantify apparent difference between sub-cellular features dedicated in a given colour channel
* One may use this pipeline to quantify morphological features and movement of the sub-cellular feature across time point experiments.

### How to execute

1. Set your parameters for background and object in distinct RGB channels

2. Make your metadata that should have the details of the cells with the filenames as unique ID

3. Images will be segregated using `maxContrastProjection` and cell segmentation will be done by `EBImage` module

4. A feature matrix will be generated - Here onwards you can either follow our instruction or anlyse by your own.

5. (Recommended) Data analysis - Image wise analysis

6. (Recommended) Data analysis - Dynamics of the orgenelle

7. (Recommended) Data analysis - Segmented Cell analysis

8. (Recommended) Data analysis - Feature based analysis

9. (Recommended) Data analysis - Clustering of features that have similar morphology


### Installation

We have a conda encapsulation of this pipeline in www.anaconda.org/arpankbasak/erb_dynamo

`conda create -n erb_dynamo r-base -y`

`conda env update --name erb_dynamo --file ./erb_dynamo.yml `

### Step by step Instructions

1. After installation, activate the conda environment

`conda activate erb_dynamo` or `source activate erb_dynamo` 

2. To quantify the dynamics 

`./dynamics.sh /path_of_the_raw_images/ /path_of_the_global_parameters/`

3. To quantify the morphology

`./morphology.sh /path_of_the_raw_images/ /path_of_the_global_parameters/`

> *"Its as simple as that !!"*

4. (Recommended yet optional) To compute statistics,

- of the dynamics of the features

`./dynamics_statistics.sh /path_to_feature_matrix_data/ /metadata/`

- and, for morphology

`./morphology_statistics.sh /path_to_feature_matrix_data/ /metadata/`

> *"Its as simple as that !!"*



### User Input

- Specify the path of image file and their location. 
- Make a metadata for statistical analysis using the `metadata.txt` file in the output

### Tidy filename

* Avoid spaces and special characters in the filename, use indices to track the metadata.
`<date>_<genotype>_<tissue>_<days_after_germination>_<batch>_<replicate>_<depth>_<time>`

---
## FAQs

1. _How many samples samples should I consider for morphological diversity analysis?_

Atleast 3 technical replicates and 2 biological replicates of the corresonding genotypes.

2. _What if there are more factors?_

One can modify the script in the corresponding directories `./scripts/` for including more factors and customise the statistical analysis. We have considered, the follwoing factors:
	- Date of experiment
	- Genotype
	- Tissue type
	- Days after germination or incubation
	- Batch as biological replicate
	- Tehnical replicate
	- Number of z-stacks
	- Time points taken (only for time-series eperiment)

3. _What if the global parameters does not fit the my object of interest?_

You can actuall make different list of parameters following the `default` parameter in the global parameter table and run the pipeline for each parameters to find the __optimum__ parameter for the object that you want to detect. In that case, you can run the scripts in parallel.

---
