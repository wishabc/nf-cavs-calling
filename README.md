# BABACHI pipeline

Nextflow pipeline for calling background allelic copies maps from SNP calls

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- babachi 

## Pipeline overview

Samples are extracted by corresponding individual and then used for a ``bcftools``-based genotyping pipeline.

## Usage
```
nextflow run main.nf -config nextflow.config -profile Altius
```