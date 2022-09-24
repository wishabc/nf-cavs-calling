# BABACHI pipeline

Nextflow pipeline for BAD maps reconstruction from SNVs with BABACHI and calling of chromatin altering variants

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- bedtools (https://bedtools.readthedocs.io/)
- babachi (v2.0.17) (https://github.com/autosome-ru/BABACHI/tree/2.0-dev)



## Usage
Before using, fill in params paths in ```babachi_params.config```. Please find detailed explanation of the parameters in the [Config section](#config).

The pipeline consists of several parts:
- Extract individual VCF files from ```all.vcf.gz``` and pre-filter SNVs for BABACHI
- BAD calling with BABACHI
- P-value calculation, aggregation and FDR correction
- CAVs filtering
- BAD calling on data with filtered out CAVs
- P-value calculation with new BADmaps, aggregation and FDR correction

To extract samples from ```all.vcf.gz``` file according to the ```indiv_id``` field of the metadata file (see ```samplesFile``` in the [Config section](#config)) and pre-filter calls for BABACHI:
```
nextflow run extract_and_filter.nf -c nextflow.config -profile Altius
```

To run all the remaining stages of the pipeline use:
```
nextflow run main.nf -c nextflow.config -profile Altius
```

## Config
There are two config files in the repository.
- ```nextflow.config``` - contains enviornment configuration. Detailed explanation can be found at https://www.nextflow.io/docs/latest/config.html. 
- ```babachi_params.config``` - specifies thresholds and paths to input files.

Following parameters should be present in ```babachi_params.config```. Each option can be specified either in ```babachi_params.config``` file or in command line using double dash (```--```) before the param name (e.g. ```--vcf_file```).
- ```vcf_file``` - path to input ```all.vcf.gz``` file with SNVs called in all samples

- ```outdir``` - directory to save results into.

- ```filteredVcfs``` - directory for results of ```extract_and_filter.nf``` script. If ```filteredVcfs``` is null or empty string, saves results to ```outdir```

- ```samplesFile``` - tab-delimited file with metadata for samples. The file must contain a header and the following columns (other columns are permitted and ignored)
    - ```indiv_id``` - unique individual ID; many samples can refer to one individual
    - ```ag_number``` - unique identifier of the sample in ```vcf_file```.<br><br>

- ```allele_tr``` - allelic reads threshold, SNVs with less than ```allele_tr``` reads on one of the alleles are filtered out

- ```fdrCovTr``` - coverage threshold. If <b>all</b> SNVs on the same position have <b>lower coverage</b> than specified, they are excluded from the analysis.

- ```excludeFdrTr``` - FDR threshold below which SNVs are called CAVs and excluded from the set of SNVs used for badmaps reestimation

Advanced params, change only if you know what you are doing
- ```states, prior, geometric_prior```,  - allowed states, prior type and coefficient for geometric prior, see https://github.com/autosome-ru/BABACHI/tree/2.0-dev for more details,

- ```esMethod``` - method of effect size calculation. Can be either ```exp``` (calculated as expected value), ```odds``` (calculated as odds ration) and ```cons``` (calculated with conservative model)
- ```recalcW``` - if true recalculates weights of the modes in the distributions mixture

