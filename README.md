# BABACHI pipeline

Nextflow pipeline for calling background allelic copies maps from SNP calls and 

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- babachi (v2.0.14)

## Pipeline overview

BAD maps reconstruction from SNPs with BABACHI and calling of allele-specific events.

## Usage
Before using, check params paths in ```babachi_params.config```.

The pipeline consists of several parts:
- Extract and filter (very time consuming, do only once and specify in filteredVCFs variable, see ```babachi_params.config```)
- BAD calling 
- P-value calculation
- CAVs filtering
- BAD calling on data with filtered out CAVs
- P-value calculation with new BADmaps

To extract from ```all.vcf``` by INDIV ID (according to the metadata file) and filter with BABACHI:
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius
```

To run all the remaining stages of the pipeline use:
```
nextflow run main.nf -c nextflow.config -profile Altius
```

