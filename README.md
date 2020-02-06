# Somatic-Germline n-of-1 DNAseq Analysis
## Call SNV and CNA from somatic and matched germline sample data
### How to Setup
#### Dependencies:
[NextFlow](https://www.nextflow.io/index.html#GetStarted), [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
#### Reference Generation:
See [DNAseq_references repo](https://github.com/brucemoran/DNAseq_references/tree/master/GRCh_37-38)
### Somatic-Germline n-of-1 Pipeline
#### To run the pipeline:
```
nextflow run brucemoran/somatic_n-of-1
```
#### Mandatory Arguments:
```
--sampleCsv      STRING      CSV format, headers: type(either "germline" or "somatic"),sampleID,meta,/path/to/read1.fastq.gz,/path/to/read2.fastq.gz
--refDir      STRING      dir in which reference data and required indices held; recommended to run associated reference creation NextFlow, DNAseq_references
```
#### Optional Arguments:
```
-profile      STRING      standard, singularity (recommended), conda (not tested/supported)
--germline      STRING      run HaplotypeCaller on germline sample and annotate with Cancer Predisposition Sequencing Reporter ([CPSR](https://github.com/sigven/cpsr))
--includeOrder      STRING      in final plots, use this ordering of samples (if multiple somatic samples); comma-separated, no spaces
--multiqcConfig      STRING      config file for multiqc
--seqLevel      STRING      "WGS" or "exome" (default)
```
#### Formats and Other Info:
For `--sampleCsv`, the format is:
```
type,sampleID,meta,read1,read2
germline,germ1,whole_blood,/full/path/to/germ1.R1.fastq.gz,/full/path/to/germ1.R2.fastq.gz
somatic,soma1,primary_tumour,/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
somatic,soma2,metastasis,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
```
The `meta` column is used for reporting in PCGR, CPSR where `sampleID` may include clinical/personal data, for example.

Headers of `sample.csv` file must match above exactly, and you should have only one germline/normal sample per run.

Column `type` must be `germline` for one sample only. In case of 2+ germline samples, run with each as `germline` in turn, specifying others as `somatic`.
