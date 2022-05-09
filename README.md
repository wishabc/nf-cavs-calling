# BABACHI pipeline

Nextflow pipeline for calling background allelic copies maps from SNP calls and 

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- babachi 

## Pipeline overview

BAD maps reconstruction from SNP calls with BABACHI and P-value estimation of allele-specific events

## Usage
### Extract and filter
Extract vcfs by INDIV (according to metadata file) and filter with babachi
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius
```

Extract all samples from vcf file (for advanced users)
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius -entry extractAllSamples
```

