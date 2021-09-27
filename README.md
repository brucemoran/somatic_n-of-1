# Somatic-Germline n-of-1 DNAseq Analysis
## Call SNV and CNA from somatic and matched germline sample data
### How to Setup
#### Dependencies:
[NextFlow](https://www.nextflow.io/index.html#GetStarted), [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
#### Reference Generation:
Run the included `download-references.nf`
```
nextflow run brucemoran/somatic_n-of-1/download-references.nf \
  -profile        [str]     singularity,refs
  --assembly      [str]     GRCh37 or GRCh38
  --exomeTag      [str]     naming for exome/panel outputs when supplied; tag is then used in somatic_n-of-1 and batch_somatic main.nf pipelines to select relevant exome reference data

  one of:
  --exomeBedURL   [str]     URL to exome bed file for intervals
  --exomeBedFile  [str]     locally downloaded exome bed file for intervals

General Optional Arguments:

  --exomeAssembly [str]     what assembly is the exome BED based on? (default: GRCH37)
  --cosmicUser    [str]     COSMIC login credentials, user email
  --cosmicPass    [str]     COSMIC login credentials, password
```
### Somatic-Germline n-of-1 Pipeline
#### About the pipeline:
This pipeline was developed to analyse and report on clinical cancer research. To this end we have tried to make a useful and clinician-readable output. This has been achieved largely on the back of [PCGR](https://github.com/sigven/pcgr)/[CPSR](https://github.com/sigven/cpsr) which provide really excellent HTML reports and annotation from multiple clinically relevant sources.

Variant calling uses an 'ensemble' approach with 3 callers (MuTect2, Manta/Strelka2 and Lancet (for exome only)) currently. More may be added in the future. The outputs are combined using our [somenone::variantconsensus](https://github.com/brucemoran/somenone/blob/master/R/somenone_variantconsensus.R) method. This takes parsed VCF output and combines calls, requiring support from at least 2 callers. It also parses 'raw' (unfiltered) calls, testing those variants that do not have support from 2 callers to see if they are contained in any single callers filtered VCF. In this way we attempt to retain as much true-positive calls as possible.
#### To run the pipeline:
```
nextflow run brucemoran/somatic_n-of-1 -profile standard,singularity
```
N.B. that currently Manta/Strelka2 and Lancet are not available through Conda. Because of this Singularity is required to run the pipeline. This may change and we will update, but Singularity is good so we recommend it anyway.
#### Arguments:
```
nextflow run brucemoran/somatic_n-of-1 --help
```
#### Formats of --sampleCsv, --sampleCat
To input specific fastq files per samples:
```
type,sampleID,meta,read1,read2
germline,germ1,whole_blood,/full/path/to/germ1.R1.fastq.gz,/full/path/to/germ1.R2.fastq.gz
somatic,soma1,primary_tumour,/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
somatic,soma2,metastasis,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
germsoma,"precancerous lesion",metastasis,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
```
To input directories containing multiple fastqs named {sampleID}*{ext}
```
type,sampleID,meta,dir,ext
germline,germ1,whole_blood,/full/path/to/dir/with/germ1/fqs,R1.fastq.gz;R2.fastq.gz
somatic,soma1,primary_tumour,/full/path/to/dir/with/soma1/fqs,R1.fastq.gz;R2.fastq.gz
somatic,soma2,metastasis,/full/path/to/dir/with/soma2/fqs,R1.fastq.gz;R2.fastq.gz
```
The `meta` column is used for reporting in PCGR, CPSR where `sampleID` may include clinical/personal data, for example.

Headers of `sample.csv` file must match above exactly, and you should have only one germline/normal sample per run.

Column `type` must be `germline` for one sample only. In case of 2+ germline samples, either run with each as `germline` in turn specifying others as `somatic`, or run one 'true' germline (e.g. blood) and other samples as 'germsoma'. This option means results for both germline CPSR results via HaplotypeCaller, and somatic PCGR results via ensemble calling, are returned for the sample. 

####Singularity Hub and containers
Unfortunately shub where the containers were housed is gone, so you will need to build the conatiners. We are working on a quay.io repo. Please pull our github Singularity repo (brucemoran/Singularity) and use our recipe_builder.sh script:
```
git clone https://github.com/brucemoran/Singularity
cd Singularity
sudo sh ./utilities/recipe_builder.sh <recipe_name> <$NXF_SINGULARITY_CACHE> shub
```
This builds a correctly named container into the $NXF_SINGULARITY_CACHE so the pipeline will run.
