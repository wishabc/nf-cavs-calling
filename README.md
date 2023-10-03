# BABACHI pipeline

Nextflow pipeline for BAD maps reconstruction from SNVs with BABACHI and calling of chromatin altering variants

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- bedtools (https://bedtools.readthedocs.io/)
- babachi (https://github.com/autosome-ru/BABACHI)



## Usage
 1) Create conda environment from `environment.yml` file with ```conda env create -n cavs-calling -f environment.yml```
 2) Modify `nextflow.config` to computing enviroment specifications
 3) Fill in params paths in ```params.config```. You can also specify parameters in command line. Please find detailed explanation of the parameters in the [Config section](#config).
 4) Run the pipeline with `nextflow run main.nf -profile Altius`
The pipeline consists of several parts:
- BAD calling with BABACHI
- P-value calculation, aggregation by individual and FDR correction
- CAVs filtering
- BAD calling on data with filtered out CAVs
- P-value calculation with new BADmaps
- Aggregation by custom aggregation key

To run all stages of the pipeline use:
```
nextflow run main.nf -profile Altius
```

To run just the last, aggregation, step (expected to run previous command first):
```
nextflow run aggregation.nf -profile Altius --main_run_outdir <path to output folder of main.nf script>
```
The `--main_run_outdir` param can be omitted if you are running the pipeline in the same folder as `nextflow run main.nf -profile Altius`

## Config
There are two config files in the repository.
- ```nextflow.config``` - contains enviornment configuration. Detailed explanation can be found at https://www.nextflow.io/docs/latest/config.html. 
- ```params.config``` - specifies thresholds and paths to input files.

Following parameters should be present in ```params.config```. Each option can be specified either in ```params.config``` file or with a command line.

- ```outdir``` - directory to save results into.
- ```conda``` - path to installed conda (from environment.yml)
- ```samples_file``` - tab-delimited file with metadata for samples. The file must contain a header and the following columns (other columns are permitted and ignored)
    - ```indiv_id``` - unique individual ID; many samples can refer to one individual. Samples with the same individual ID are merged for BADmaps calculation.
    - ```ag_id``` - unique identifier of the sample.<br><br>
    - `snps_file` - (not required for aggregation workflow) path to the bed-formatted file with SNVs and their readcounts. See [babachi](https://github.com/autosome-ru/BABACHI) for more details.
    - `footprints_path` - (optional, not required for aggregation workflow) Path to bed-formatted footprint calls
    - `hotspots_path` - (optional, not required for aggregation workflow) Path to bed-formatted peak calls

- `aggregation_key` - column name in `samples_file`. Samples are grouped according to values in that columns and aggregated. Use `"all"` to aggregate all the data;

- ```max_coverage_tr``` - coverage threshold. SNPs with less than `max_coverage_tr` are not considered to be CAV candidates and excluded from aggregation.

- ```fdr_tr``` - FDR threshold below which SNVs are called CAVs and excluded from the set of SNVs used for badmaps re-estimation

BABACHI params, change only if you know what you are doing
- ```states, prior, geometric_prior```  - allowed states, prior type and coefficient for geometric prior, see https://github.com/autosome-ru/BABACHI/ for more details
- ```allele_tr``` - exclude SNPs with less than `alelle_tr` reads on each of the alleles for BAD maps estimation, tested on `allele_tr=5`
- ```babachi_maf_tr``` - exclude SNPs with less than `babachi_maf_tr` AF on each of the alleles for BAD maps estimation, tested on `babachi_maf_tr=0.05`
