# MicroarrayPipeline

## Requirements
R with these packages installed: affy, limma, sva, org.Hs.eg.db, ggplot2, [hgu133plus2hsrefseqcdf](http://mbni.org/customcdf/22.0.0/refseq.download/hgu133plus2hsrefseqcdf_22.0.0.tar.gz) (Or download from [here](http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/22.0.0/refseq.asp))

## Usage
This is the command line that to excute this script:
```
Rscript affy_limma.R data_path design_matrix_path result_path 1 0.05
```
*Only affy_limma.R is available now.*  

The first argument is the data files path, the second is the design matrix path.  

The design matrix should include two columns and delaminated by tab. The first column is GSM IDs of RNAseq samples. The second column is the design of experiment. Usually control group shoud be upper and treatment group should be at the bottom. For example:
```
GSM10001	Control
GSM10002	Control
GSM10003	Treatment
GSM10004	Treatment
```
The third argument is the results path or output path, and is supposed to be a folder.  

The fourth is fold change and the fifth is adjust P value.
