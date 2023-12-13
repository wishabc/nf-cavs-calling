# CAV calling pipeline

Nextflow pipeline for calling of chromatin altering variants (CAVs).

## Requirements
- Nextflow (https://www.nextflow.io/)
- bcftools (http://www.htslib.org/)
- bedtools (https://bedtools.readthedocs.io/)
- babachi (https://github.com/autosome-ru/BABACHI)


## Usage
 0) (Optional) Create conda environment from `environment.yml` file with ```conda env create -n cavs-calling -f environment.yml```
 1) Modify `nextflow.config` to computing enviroment specifications
 2) Fill in params paths in ```params.config```. You can also specify parameters in command line. Please find detailed explanation of the parameters in the [Config section](#config).
 3) Run the pipeline with `nextflow run main.nf -profile Altius`
The pipeline consists of several parts:
- First round BAD calling with BABACHI
- P-value calculation
- BAD calling on data with filtered out first round CAVs
- P-value calculation with second round of BAD maps
- Aggregation by custom aggregation key
- Sorting, tabix indexing and collecting of statistics for output files

To run all stages of the pipeline use:
```
nextflow run main.nf -profile Altius -resume
```

To run just the last, aggregation, step (expected to run previous command first):
```
nextflow run aggregation.nf -profile Altius --main_run_outdir <path to output folder of main.nf script>
```
The `--main_run_outdir` param can be omitted if you are running the pipeline in the same folder as `nextflow run main.nf -profile Altius`. The output files are named according to aggregation key to avoid name collisions.

## Config
There are two config files in the repository.
- ```nextflow.config``` - contains enviornment configuration. Detailed explanation can be found at https://www.nextflow.io/docs/latest/config.html. 
- ```params.config``` - specifies thresholds and paths to input files.

Following parameters should be present in ```params.config```. Each option can be specified either in ```params.config``` file or with a command line.

- ```outdir``` - directory to save results into. Defaults to `output` folder in the launch directory
- ```conda``` - (optional) path to installed conda (from environment.yml). If not present, nextflow creates environment from environment.yml (was not tested).
- ```samples_file``` - tab-delimited file with metadata for samples. The file must contain a header and the following columns (other columns are permitted and ignored)
    - ```indiv_id``` - unique individual ID; many samples can refer to one individual. Samples with the same individual ID are merged for BADmaps calculation.
    - ```ag_id``` - unique identifier of the sample.<br><br>
    - `snps_file` - (not required for aggregation workflow) path to the bed-formatted file with SNVs and their readcounts. See [babachi](https://github.com/autosome-ru/BABACHI) for more details.
    - `footprints_path` - (optional) Path to bed-formatted footprint calls
    - `hotspot_peaks_point1per` - (optional) Path to bed-formatted peak calls

- `aggregation_key` - column name in `samples_file`. Samples are grouped according to values in that columns and aggregated. Use `"all"` to aggregate all the data. Samples with NA values in that column are ignored.

- ```max_coverage_tr``` - coverage threshold. Group of SNVs with coverage less than `max_coverage_tr` at *all* of them is excluded from aggregation.

- ```fdr_tr``` - FDR threshold below which SNVs are called CAVs. Same threshold is used to excluded first level CAVs for second round of badmaps estimation.


BABACHI params, change only if you know what you are doing:
- ```states, prior, geometric_prior```  - allowed states, prior type and coefficient for geometric prior, see https://github.com/autosome-ru/BABACHI/ for more details
- ```babachi_allele_tr``` - exclude SNPs with less than `alelle_tr` reads on each of the alleles for BAD maps estimation, tested on `allele_tr=5`
- ```babachi_maf_tr``` - exclude SNPs with less than `babachi_maf_tr` AF on each of the alleles for BAD maps estimation, tested on `babachi_maf_tr=0.05`
