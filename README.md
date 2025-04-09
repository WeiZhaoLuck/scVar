scVar:
===========


Install
===========



Dependency data
===========
The dependency data can be downloaded from these links and should be placed in the <reference_path> folder.
- **Vep**:
- **Annovar**:


Usage
===========
## Preparatory work
Please prepare the input files in advance. Follow the instructions below based on your data type.

### 1. Raw Data (FASTQ Files)
If you are providing raw sequencing data, ensure the FASTQ files are named according to the following convention and move the files  to the `<data_path>` directory, like : `./Demo/Raw/*.fastq.gz`

[SampleName]S1_L00[LaneNumber][ReadType]_001.fastq.gz  

### 2. Processed Data (BAM, h5ad, or rds Files)

If you are skipping upstream analysis, prepare the following processed files:
- **BAM files**: Ensure the BAM files are named according to the following convention and move the files to the `<data_path>` directory, like : `./Demo/Processed/*.bam`
- **h5ad files**: Ensure the h5ad files are named according to the following convention and move the files to the `<data_path>` directory, like : `./Demo/Processed/*.h5ad`
- **Seurat files**: Ensure the Seurat files are named according to the following convention and move the files to the `<data_path>` directory, like : `./Demo/Processed/*.rds`

## Run scVar
### 1. Setup Configuration
Create a configuration file `./scVar/config.yaml` in the `<results_path>` directory. This file should contain the following information:

```yaml
# Configuration file for scVar analysis
# Please fill in the required information below
project:
    <ProjectName> # Project name
samples:
    <SampleName>: SamplePath # Sample name and path to the sample data
result_path:
    <results_path> # Path to save the results
SNV_filter_vaf:
    100
threads:
    7
genotype_filter:
    0
genotype_mapq:
    0
genotype_baseq:
    0
```
### 2. Copy Snakefile
Copy a snakefile file in the `<results_path>` directory.
- **For Raw Data**: Copy `./scVar/Snakefile_Raw` and rename it to Snakefile.
- **For Processed Data**: Copy `./scVar/Snakefile` and rename it to Snakefile.

### 3. Run scVar
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && cd <results_path> && snakemake --cores <cores_number>'
```

- `<reference_path>`: Path to the reference genome directory.
- `<results_path>`: Path to the results directory.
- `<data_path>`: Path to the data directory.
- `<cores_number>`: Number of cores to use for the analysis.

## Upstream Analysis
### 1. Mutation signature, TMB, Entropy
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Analysis_SigTMBEntroy.sh <project_path> <sample_name>'
```
- `<project_path>`: Path to the project directory.
- `<sample_name>`: Name of the sample to analyze.

### 2. Celltype Specific Mutations
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Calculate_Specific_Mutations.sh <project_path> <sample_name>'
```
- `<project_path>`: Path to the project directory.
- `<sample_name>`: Name of the sample to analyze.

### 3. Mutation Clustering
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/MutationCluster.sh --path <project_path> --sample <sample_name> --method <method> --flag <flag> --number <number> --clustermethod <clustermethod>'
```
- `<project_path>`: Path to the project directory.
- `<sample_name>`: Name of the sample to analyze.
- `--method`: Method for clustering (default: TF_IDF).
- `--flag`: Flag for clustering (default: 0).
- `--number`: Number of clusters (default: 100).
- `--clustermethod`: Clustering method (default: complete).

### 4. GO and Oncogenic Pathway
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/GOonco.sh --path <project_path> --sample <sample_name> --pCutoff <pCutoff> --qCutoff <qCutoff>'
```
- `<project_path>`: Path to the project directory.
- `<sample_name>`: Name of the sample to analyze.
- `--pCutoff`: P-value cutoff for GO analysis (default: 0.05).
- `--qCutoff`: Q-value cutoff for GO analysis (default: 0.2).

### 5. Psedutime Analysis
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Pseudotime.sh --path <project_path> --sample <sample_name> --mutation <mutations>'
```
- `<project_path>`: Path to the project directory.
- `<sample_name>`: Name of the sample to analyze.
- `--mutation`: Mutation IDs connected by commas, e.g.,12_25245350_C_G,17_47592542_A_G.

### 6. Statistic
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Statistic.sh --path <project_path> --sample <sample_name>'
```
- `<project_path>`: Path to the project directory.
- `<sample_name>`: Name of the sample to analyze.

For detailed information on inputs, outputs, scripts, and parameters of each module, refer to [documentation/section].

## Example Demo
### Preparatory work
- **Copy Input Files**: Copy files from either `./Demo/Raw/` or `./Demo/Processed/` to the <data_path> directory.
- **Copy Configuration File**: If you copy data from `./Demo/Raw/`, copy `./Demo/Raw_config/config.yaml` and `./Demo/Raw_config/Snakefile` to the `<results_path>` directory. If you copy data from `./Demo/Processed/`, copy `./Demo/Processed_config/config.yaml` and `./Demo/Processed_config/Snakefile` to the `<results_path>` directory.


### Run scVar
    
```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && cd /results && snakemake --cores 1'
```

### Upstream Analysis
#### 1. Mutation signature,TMB,Entropy

```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Analysis_SigTMBEntroy.sh /results/Example Demo'
```

#### 2. Celltype Specific Mutations

```shell
docker run -it -v <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Calculate_Specific_Mutations.sh /results/Example Demo'
```


#### 3. Mutation Clustering

```shell
docker run -it <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/MutationCluster.sh --path /results/Example --sample Demo --method TF_IDF --flag 0 --number 100 --clustermethod complete'
```

#### 4. GO and Oncogenic Pathway

```shell
docker run -it <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/GOonco.sh --path /results/Example --sample Demo  --pCutoff 0.05 --qCutoff 0.2'
```

#### 5. Psedutime Analysis

```shell
docker run -it <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Pseudotime.sh --path /results/Example --sample Demo --mutation 7_6004027_A_G,12_25245350_C_G,17_47592542_A_G,22_20708085_C_T,2_28942361_C_A,3_49684128_A_C,3_49684173_G_A,3_49684565_A_G,7_56116040_A_C,9_137717161_G_T,9_137717162_C_A,X_100662268_C_G'
```


#### 6. Statistic

```shell
docker run -it <reference_path>:/reference -v <results_path>:/results  -v <data_path>:/data scvar /bin/bash -c 'source /opt/miniconda/bin/activate scVar && bash /codes/Statistic.sh /results/Example Demo'
```


## Citation:
