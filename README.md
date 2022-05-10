# BABACHI pipeline

Nextflow pipeline for calling background allelic copies maps from SNP calls and 

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- babachi 

## Pipeline overview

BAD maps reconstruction from SNPs with BABACHI and calling of allele-specific events.

## Usage
Before using, check params paths in ```babachi_params.config```.

The pipeline consists of three parts:
- Extract and filter
- BAD calling
- P-value calculation

To run all stages of the pipeline use
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius
```

<br>
Below you can find detailed description of each pipeline stage:<br><br> 

### Extract and filter
Extract vcfs by INDIV (according to metadata file) and filter with babachi
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius
```

Extract all samples from vcf file (for advanced users)
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius -entry extractAllSamples --outdir <out_dir>
```
Filter all samples using babachi filter (for advanced users)
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius -entry extractAllSamples --outdir <out_dir>
```
### BAD calling
Estimate BAD with BABACHI and intersect resulting BADmaps with SNP calls
```
nextflow run bad_estimation.nf -c nextflow.config -profile Altius
```
### P-value estimation

