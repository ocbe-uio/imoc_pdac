# Biomarker discovery and patient stratification in pancreatic cancer using incomplete multi-omics data
This repository contains the data and scripts used in the paper "Biomarker discovery and patient stratification in pancreatic cancer using incomplete multi-omics data". 

## Data download
To download the data used in this study, run the R script `/data/download_data.R`. To process the data, run the script `/data/data_preprocessing.py` twice, once with `complete_sample_set = False` and `complete_sample_set = True` to generate the entire data used in the study.

## Benchmark steps
To obtain the results found in the `/results/cluster_analysis/benchmarking_files` folder, run the file `/src/scripts/generating_indxs.py` first, followed by `/src/scripts/incomplete_algorithms_evaluation.py`. The settings must be changed for every benchmark in the `settings.py` file. However, this is not adviced as it can take several hours to run every script. 


### First benchmark
```
best_combination = False
run_amputation = False
select_datasets = ['patients_with_all_views']
n_clusters = [2, 3, 4, 5]   # can be done for each number of clusters individually, e.g.: 2 clusters, 3 clusters, etc.
runs_per_alg = np.arange(10)
sampling = True
```

### Second benchmark
```
best_combination = ["CNA", "Methylation"]
run_amputation = True
select_datasets = ['patients_with_all_views']
n_clusters = [2]
runs_per_alg = np.arange(10)
sampling = False
```

### Final clusters
```
best_combination = ["CNA", "Methylation"]
run_amputation = False
select_datasets = ['all_patients']
n_clusters = [2]
runs_per_alg = np.arange(1024)
sampling = True
```

## Cluster analysis
The analysis of the final clusters can be reproduced by running the jupyter notebooks in the `/notebooks` folder. This will generate the data used for the figures and supplementary data. The files 

## Omics Analysis

To reproduce the omics analysis, open the `.Rproj` file in R. Then, run the R scripts in the `/omics_analysis` folder. This will generate the data used for the figures and supplementary data.

To reproduce the manuscript's figures and supplementary data, run the script `/omics_analysis/05_Results_Figure.R`