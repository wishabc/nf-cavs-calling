# BABACHI pipeline

Nextflow pipeline for calling background allelic copies maps from SNP calls

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- babachi 

## Pipeline overview

Samples are extracted by corresponding individual and then used for BAD maps reconstruction. 

## Usage
```
nextflow run main.nf -config nextflow.config -profile Altius
```