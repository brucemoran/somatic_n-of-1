#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                          SOMATIC_N-OF-1 PIPELINE
  -----------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/somatic_n-of-1

  Mandatory arguments:

    -profile        [str]       Configuration profile
                                (required: standard,singularity)

    --sampleCsv     [file]      CSV format, headers: type (either "germline" or
                                "somatic"), sampleID, meta, read1 (e.g.
                                /path/to/read1.fastq.gz), read2
                                (e.g. /path/to/read2.fastq.gz); use meta for
                                naming in PCGR, CPSR reports

    --runID         [str]       Name for run, used to tag outputs

    --refDir        [file]      Path of dir in which reference data are held;
                                this should be created by download-references.nf
                                and contain dir <assembly>

    --assembly      [str]       Either GRCh37 or GRCh38 (default), as per
                                download-references.nf

    --email         [str]       Email address to send reports

  General Optional Arguments:

    --germline      [bool]      Run HaplotypeCaller on germline sample and
                                annotate with CPSR (default: true)

    --germCNV      [bool]       Run CNVkit germine CNV caller (default: false)

    --scatGath      [int]       Number of pieces to divide intervalList into for
                                scattering to variant calling processes
                                (default: 20 for exome, 100 for WGS)

    --incOrder      [str]       In final plots, use this ordering of samples
                                (if multiple somatic samples); comma-separated,
                                no spaces (default: alphanumeric sort)

    --sampleCat     [str]       File used when fastq data is in multiple files
                                which are cat'ed; replaces --sampleCsv; headers:
                                type (either germline or somatic), sampleID,
                                meta, dir (contains fastq to be cat'ed), ext
                                (extension scheme for parsing read1, read2 e.g.
                                _1.fq.gz;_2.fq.gz).

    --bamCsv        [file]      CSV format, headers as sampleCsv except read1,
                                read2 are swapped for bam which is sent to
                                duplicate marking

    --multiqcConfig [str]       Config file for multiqc
                                (default: bin/somatic_n-of-1.multiQC_config.yaml)

    --seqLevel      [str]       WGS or exome (default: WGS)

    --exomeTag      [str]       Tag used for exome kit when running download-references.nf

    --cosmic        [bool]      set this to specify output of COSMIC CGC genes
                                only (somatic only; based on download and supply
                                of CGC file in download_references.nf)

    --phylogeny     [bool]      conduct subclonal phylogeny reconstruction with
                                pairtree
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

//Test Mandatory Arguments
if(params.sampleCsv && params.sampleCat){
  exit 1, "Please include only one of --sampleCsv or --sampleCat, see --help for format"
}

if(params.sampleCsv == null && params.sampleCat == null && params.bamCsv == null){
  exit 1, "Please include one of --sampleCsv, --bamCsv or --sampleCat, see --help for format"
}

if(!Channel.from(params.runID, checkIfExists: true)){
    exit 1, "Please include --runID <your_runID>"
}

if(!Channel.from(params.refDir, checkIfExists: true)){
  exit 1, "Please include --refDir <path> see github.com/brucemoran/somatic_n-of-1/ for how to run download-references.nf"
}

if(!Channel.from(params.assembly, checkIfExists: true)){
    exit 1, "Please include --assembly <GRCh3x>"
}

if(!params.email){
    exit 1, "Please include --email your@email.com"
}

if(params.seqLevel == "exome" && params.exomeTag == null){
    exit 1, "Please define --exomeTag when using --seqLevel exome"
}
//Global Variables based on input
params.outDir = "${params.seqLevel}_output"
params.seqlevel = "${params.seqLevel}".toLowerCase()

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

//Reference data as value channels and reusable therefore
reference = [
    grchvers: false,
    fa: false,
    fai: false,
    dict: false,
    bwa: false,
    hc_dbs: false,
    dbsnp: false,
    gridss: false,
    pcgrbase: false,
    intlist: false,
    seqlevel: false,
    bbres: false
]

reference.grchvers  = Channel.fromPath("${params.refDir}/${params.assembly}/pcgr/data/*", type: 'dir').getVal()
reference.fa = Channel.value(file(params.genomes[params.assembly].fa))
reference.fai = Channel.value(file(params.genomes[params.assembly].fai))
reference.dict = Channel.value(file(params.genomes[params.assembly].dict))
reference.bwa = Channel.value(file(params.genomes[params.assembly].bwa))
reference.hc_dbs = Channel.value(file(params.genomes[params.assembly].hc_dbs))
reference.dbsnp = Channel.value(file(params.genomes[params.assembly].dbsnp))
reference.gridss = Channel.value(file(params.genomes[params.assembly].gridss))
reference.pcgrbase = Channel.value(file(params.genomes[params.assembly].pcgr))
reference.refflat = Channel.value(file(params.genomes[params.assembly].refflat))

if(params.microbiome){
  reference.pathseq = Channel.value(file(params.genomes[params.assembly].pathseq))
}

//if seqlevel is exome, there is a dir per exome already parsed according to exomeTag
reference.seqlevel = params.seqlevel == "wgs" ? Channel.value(file(params.genomes[params.assembly].wgs)) : Channel.value(file(params.genomes[params.assembly].exome))

//set cosmic
reference.cosmic = params.cosmic == true ? Channel.value(file(params.genomes[params.assembly].cosmic)) : null

//setting of intlist, bed based on seqlevel and exomeTag
reference.intlist = params.seqlevel == "wgs" ? Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/wgs.bed.interval_list").getVal() : Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/${params.exomeTag}/${params.exomeTag}.bed.interval_list").getVal()
reference.bed = params.seqlevel == "wgs" ? Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/wgs.bed").getVal() : Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/${params.exomeTag}/${params.exomeTag}.bed").getVal()

/*
================================================================================
                          -0. PREPROCESS INPUT SAMPLE FILE
================================================================================
*/
/* 0.00: Input using sample.csv, bam.csv, sample_cat
*/
if(params.sampleCsv){
  Channel.fromPath("${params.sampleCsv}")
         .splitCsv( header: true )
         .map { row -> [row.type, row.sampleID, row.meta, file(row.read1), file(row.read2)] }
         .into { bbduking; ubaming }
}

if(params.bamCsv){
  Channel.fromPath("${params.bamCsv}")
         .splitCsv( header: true )
         .map { row -> [row.type, row.sampleID, row.meta, file(row.bam)] }
         .set { which_bam }

  process bam_input {

    label 'low_mem'

    input:
    tuple val(type), val(sampleID), val(meta), file(bam) from which_bam

    output:
    tuple val(type), val(sampleID), val(meta), file(bam), file("*.bai") into dup_marking

    script:
    """
    #! bash
    samtools index ${bam}
    """
  }
}

if(params.sampleCat){
  Channel.fromPath("${params.sampleCat}")
         .splitCsv( header: true )
         .map { row -> [row.type, row.sampleID, row.meta, row.dir, row.ext] }
         .set { samplecating }

  process samplecat {

    label 'low_mem'
    publishDir "${params.outDir}/samples/${sampleID}/cat", mode: "copy"

    input:
    tuple val(type), val(sampleID), val(meta), val(dir), val(ext) from samplecating

    output:
    tuple val(type), val(sampleID), val(meta), file(read1), file(read2) into ( bbduking, ubaming )

    script:
    rd1ext = "${ext}".split(';')[0]
    rd2ext = "${ext}".split(';')[1]
    read1 = "${sampleID}.R1.fastq.gz"
    read2 = "${sampleID}.R2.fastq.gz"
    """
    #! bash
    cat \$(find ${dir} | grep ${rd1ext} | sort) > ${read1}
    cat \$(find ${dir} | grep ${rd2ext} | sort) > ${read2}
    """
  }
}

// 0.01: Create a uBAM for Pathseq
if(!params.bamCsv){
  process ubam {

    label 'med_mem'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from ubaming

    output:
    tuple val(type), val(sampleID), val(meta), file('*.bam'), file('*.bai') into pathseqing

    when:
    params.microbiome == true

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    """
    DATE=\$(date +"%Y-%m-%dT%T")
    mkdir tmp
    {
    picard ${taskmem} FastqToSam \
      FASTQ=${read1} \
      FASTQ2=${read2} \
      OUTPUT=${sampleID}.unaligned.bam \
      READ_GROUP_NAME=${sampleID} \
      SAMPLE_NAME=${sampleID} \
      LIBRARY_NAME=LANE_X \
      PLATFORM_UNIT=IL_X \
      PLATFORM=ILLUMINA \
      SEQUENCING_CENTER=UCD \
      RUN_DATE=\$DATE \
      TMP_DIR="tmp"
    } 2>&1 | tee > ${sampleID}.FastqToSam.log.txt

    samtools index ${sampleID}.unaligned.bam

    rm -r tmp
    """
  }

  /*
  ================================================================================
                            0. PREPROCESS INPUT FASTQ
  ================================================================================
  */
  // 0.1: Input trimming
  process bbduk {

    label 'med_mem'
    publishDir path: "${params.outDir}/samples/${sampleID}/bbduk", mode: "copy", pattern: "*.txt"

    input:
    tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from bbduking

    output:
    file('*.txt') into log_bbduk
    tuple val(type), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming
    tuple val(type), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz'), file(read1), file(read2) into fastping
    tuple val(type), val(sampleID), val(meta), file(read1), file(read2) into fastqcing

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    """
    {
    sh bbduk.sh ${taskmem} \
      in1=${read1} \
      in2=${read2} \
      out1=${sampleID}".bbduk.R1.fastq.gz" \
      out2=${sampleID}".bbduk.R2.fastq.gz" \
      k=31 \
      mink=5 \
      hdist=1 \
      ktrim=r \
      trimq=20 \
      qtrim=rl \
      maq=20 \
      ref=/opt/miniconda/envs/somatic_n-of-1/opt/bbmap-adapters.fa \
      tpe \
      tbo \
      stats=${sampleID}".bbduk.adapterstats.txt" \
      overwrite=T
    } 2>&1 | tee > ${sampleID}.bbduk.runstats.txt
    """
  }

  // 0.2: fastp QC of pre-, post-bbduk
  process fastp {

    label 'low_mem'
    publishDir "${params.outDir}/samples/${sampleID}/fastp", mode: "copy", pattern: "*.html"

    input:
    tuple val(type), val(sampleID), val(meta), file(preread1), file(preread2), file(postread1), file(postread2) from fastping

    output:
    file('*.html') into fastp_html
    file('*.json') into fastp_multiqc

    script:
    """
    fastp -w ${task.cpus} -h ${sampleID}"_pre.fastp.html" -j ${sampleID}"_pre.fastp.json" --in1 ${preread1} --in2 ${preread2}

    fastp -w ${task.cpus} -h ${sampleID}"_post.fastp.html" -j ${sampleID}"_post.fastp.json" --in1 ${postread1} --in2 ${postread2}
    """
  }

  // 0.3: fastQC of per, post-bbduk
  process fastqc {

    label 'low_mem'
    publishDir "${params.outDir}/samples/${sampleID}/fastqc", mode: "copy", pattern: "*.html"

    input:
    tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from fastqcing

    output:
    file('*.html') into fastqc_multiqc

    script:
    """
    #!/bin/bash
    fastqc -t ${task.cpus} --noextract -o ./ ${read1} ${read2}
    """
  }

  /*
  ================================================================================
                          1. ALIGNMENT AND BAM PROCESSING
  ================================================================================
  */
  // 1.0: Input alignment
  process bwamem {

    label 'high_mem'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from bwa_memming
    file(bwa) from reference.bwa

    output:
    tuple val(type), val(sampleID), val(meta), file('*.bam'), file('*.bai') into (cramming, dup_marking)

    script:
    def fa = "${bwa}/*fasta"
    """
    DATE=\$(date +"%Y-%m-%dT%T")
    RGLINE="@RG\\tID:${sampleID}\\tPL:ILLUMINA\\tSM:${sampleID}\\tDS:${type}\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

    bwa mem \
      -t${task.cpus} \
      -M \
      -R \$RGLINE \
      ${fa} \
      ${read1} ${read2} | \
      samtools sort -T "tmp."${sampleID} -o ${sampleID}".sort.bam"
    samtools index ${sampleID}".sort.bam"
    """
  }

  // 1.1: CRAM alignment and output
  // TODO: output upload schema for ENA/EGA
  process cram {

    label 'low_mem'
    publishDir path: "${params.outDir}/samples/${sampleID}/bwa", mode: "copy", pattern: "*.cra*"

    input:
    tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from cramming
    file(bwa) from reference.bwa

    output:
    tuple file('*.cram'), file('*.crai') into completedcram

    script:
    """
    samtools view -hC -T ${bwa}/*fasta ${sampleID}".sort.bam" > ${sampleID}".sort.cram"
    samtools index ${sampleID}".sort.cram"
    """
  }
}

// 1.2: MarkDuplicates
process mrkdup {

  label 'high_mem'
  publishDir path: "${params.outDir}/samples/${sampleID}/picard", mode: "copy", pattern: "*.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from dup_marking

  output:
  file('*.txt') into mrkdup_multiqc
  tuple val(type), val(sampleID), val(meta), file('*.md.bam'), file('*.md.bam.bai') into ( gatk4recaling, gridssing )

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  OUTBAM=\$(echo ${bam} | sed 's/bam/md.bam/')
  OUTMET=\$(echo ${bam} | sed 's/bam/md.metrics.txt/')
  {
  picard ${taskmem} \
    MarkDuplicates \
    TMP_DIR=./ \
    INPUT=${bam} \
    OUTPUT=/dev/stdout \
    COMPRESSION_LEVEL=0 \
    QUIET=TRUE \
    METRICS_FILE=\$OUTMET \
    REMOVE_DUPLICATES=FALSE \
    ASSUME_SORTED=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    VERBOSITY=ERROR | samtools view -Shb - > \$OUTBAM

  samtools index \$OUTBAM
  } 2>&1 | tee > ${sampleID}.picard_markDuplicates.log.txt
  """
}

// 1.3: GATK4 BestPractices
process gtkrcl {

  label 'high_mem'
  publishDir path: "${params.outDir}/samples/${sampleID}/gatk4/bestpractice", mode: "copy", pattern: "*.GATK4_BQSR.log.txt "

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from gatk4recaling
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp
  file(intlist) from reference.intlist

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( germfiltering, gmultimetricing, mosdepthing)
  tuple val(type), val(sampleID), val(meta), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( hc_germ, cnvgerm )
  tuple val(sampleID), val(meta) into metas_pcgr
  file("${sampleID}.GATK4_BQSR.log.txt") into bqsr_log

  script:
  def dbsnp = "${dbsnp_files}/*gz"
  """
  {
  ##ensure seq dict from BAM has same regions as bed interval list
  head -n1 ${intlist} > use.interval_list
  samtools view -H ${bam} | grep "@SQ" | cut -f 2 > bam_sq.txt
    cat bam_sq.txt | while read SEQ; do
      export seq=\${SEQ};
      perl -ane 'if(\$F[1] eq \$ENV{'seq'}){print \$_;}' ${intlist};
    done >> use.interval_list
  sed 's/SN://g' bam_sq.txt | while read CHR; do
    export chr=\${CHR};
    perl -ane 'if(\$F[0] eq \$ENV{'chr'}){print \$_;}' ${intlist};
  done >> use.interval_list

  gatk BaseRecalibrator \
    -R ${fasta} \
    -I ${bam} \
    --known-sites \$(echo ${dbsnp}) \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    --disable-sequence-dictionary-validation true \
    -L use.interval_list

  #ApplyBQSR
  OUTBAM=\$(echo ${bam} | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R ${fasta} \
    -I ${bam} \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L use.interval_list

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee >  ${sampleID}.GATK4_BQSR.log.txt
  """
}

// 1.4: scatter-gather implementation for mutect2, lancet
process scat_gath {

  label 'low_mem'

  input:
  file(intlist) from reference.intlist

  output:
  file('lancet.scatgath.*.bed') into lancet_bedding
  file('mutect2.scatgath.*.bed.interval_list') into mutect2_bedding
  file('hc.scatgath.*.bed.interval_list') into hc_bedding

  script:
  def sgcount = params.scatGath
  if (params.scatGath == null){
    if (params.seqlevel == "exome"){
      sgcount = 20
    }
    else {
      sgcount = 100
    }
  }
  """
  ##strip out all but chromosomes in the interval_list (no decoys etc)
  CHRS=\$(grep -v "@" ${intlist} | cut -f 1 | uniq)
  for CHR in \$CHRS; do
    grep "SN:\$CHR\\s" ${intlist} >> used.interval_list
  done
  grep -v "@" ${intlist} >> used.interval_list

  ##generate scatters
  picard IntervalListTools \
    I=used.interval_list \
    SCATTER_COUNT=${sgcount} \
    O=./

  ##rename scatters and parse into appropriate format for tools
  ls temp*/* | while read FILE; do
    COUNTN=\$(dirname \$FILE | perl -ane '@s=split(/\\_/); print \$s[1];');
    mv \$FILE mutect2.scatgath.\${COUNTN}.bed.interval_list;
    cp mutect2.scatgath.\${COUNTN}.bed.interval_list hc.scatgath.\${COUNTN}.bed.interval_list
    grep -v @ mutect2.scatgath.\${COUNTN}.bed.interval_list | \
      cut -f 1,2,3,5 > lancet.scatgath.\${COUNTN}.bed
  done
  """
}

process Mosdepth {

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from mosdepthing
  file(bed) from reference.bed

  output:
  file('*') into mosdepth_multiqc

  script:
  """
  mosdepth \
    --no-per-base \
    --by ${bed} \
    ${sampleID} \
    ${bam}
  """
}
/*
================================================================================
                            2.  MUTATION CALLING
================================================================================
*/
// 2.0: GATK4 Germline Haplotypecaller
// Groovy to combine scatter-gather BEDs with bam file for germline
hcbedding = hc_bedding.flatten()
hc_germ
  .map { it -> [it[0],it[1],it[2],it[3],it[4]] }
  .combine(hcbedding)
  .set { hcgermbedding }

process haplotypecaller {

  label 'med_mem'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai), file(intlist) from hcgermbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp
  file(hc_dbs_files) from reference.hc_dbs

  output:
  tuple val(sampleID), val(meta), file('*sort.hc.vcf') into hc_gt
  tuple val(type), val(sampleID) into ver_germID

  when:
  type == "germline" & params.germline != false \
   | type == "germsoma" & params.germline != false

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  def dbsnp = "${dbsnp_files}/*gz"
  def omni = "${hc_dbs_files}/KG_omni*.gz"
  def kgp1 = "${hc_dbs_files}/KG_phase1*.gz"
  def hpmp = "${hc_dbs_files}/hapmap*.gz"
  """
  SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print \$s[2];')
  gatk ${taskmem} HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp \$(echo ${dbsnp}) \
    --native-pair-hmm-threads ${task.cpus} \
    -O ${sampleID}".\${SCATGATHN}.hc.vcf" \
    --disable-sequence-dictionary-validation true \
    -L ${intlist}

  picard SortVcf \
    I=${sampleID}".\${SCATGATHN}.hc.vcf" \
    O=${sampleID}".\${SCATGATHN}.sort.hc.vcf" \
    SD=${dict}
  """
}

//in case of germsoma, verify the true germID
//only one with germline as type
ver_germID
  .filter( { it[0] == "germline" } )
  .map( { it -> it[1] } )
  .into { gridssgermID; vcfGRaID }

//group those outputs
hc_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1].unique(), it[2..-1].flatten()) }
  .set { hc_fm }

/* 2.0.5: GATK germline CNV
*/

process germCnvkit {

  publishDir path: "${params.outDir}/samples/${sampleID}/cnvkit", mode: "copy"
  label 'med_mem'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from cnvgerm
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(refflat) from reference.refflat
  file(bed) from reference.bed

  output:
  file('*') into germCnvkit_comp
  tuple file("${sampleID}.call.cns"), file("${sampleID}.scatter.pdf"), file("${sampleID}.diagram.pdf") into sendmail_cnvkit

  when:
  type == "germline" & params.germline != false & params.germCNV != false \
   | type == "germsoma" & params.germline != false & params.germCNV != false

  script:
  seqlev = params.seqlevel == "wgs" ? "wgs" : "hybrid"
  """
  cnvkit.py access ${fasta} -o access.bed
  cnvkit.py autobin -f ${fasta} \
                    -m ${seqlev} \
                    -t access.bed \
                    --annotate ${refflat} \
                    --short-names \
                    --target-output-bed target.bed \
                    --antitarget-output-bed antitarget.bed \
                    ${bam}

  # For each sample...
  cnvkit.py coverage -o ${sampleID}.targetcoverage.cnn \
                     ${bam} target.bed
  cnvkit.py coverage -o ${sampleID}.antitargetcoverage.cnn \
                     ${bam} antitarget.bed

  # reference
  cnvkit.py reference -f ${fasta} \
                      -o reference.cnn \
                      ${sampleID}.antitargetcoverage.cnn \
                      ${sampleID}.targetcoverage.cnn

  # cnr
  cnvkit.py fix -o ${sampleID}.cnr ${sampleID}.targetcoverage.cnn ${sampleID}.antitargetcoverage.cnn reference.cnn

  # cns
  cnvkit.py segment -o ${sampleID}.cns ${sampleID}.cnr

  #call
  cnvkit.py call -o ${sampleID}.call.cns ${sampleID}.cns

  # Optionally, with --scatter and --diagram
  cnvkit.py scatter -s ${sampleID}.cns -o ${sampleID}.scatter.pdf ${sampleID}.cnr
  cnvkit.py diagram -s ${sampleID}.cns -o ${sampleID}.diagram.pdf ${sampleID}.cnr
  """
}

/* 2.1.1: GRIDSS SV calling in WGS
* because we do not know order or number of samples, create tuple with
* dummy as first and all others as list of elements
* then count them inside process
*/
gridssing
  .collect()
  .map { it -> tuple(it.flatten()) }
  .set { gridssin }

process gridss {

  label 'max_mem'
  publishDir path: "${params.outDir}/combined/gridss", mode: "copy", pattern: "*.[!bam, vcf.gz]"

  input:
  file(listbams) from gridssin
  val(germlineID) from gridssgermID.collect().flatten().unique()
  file(bwa) from reference.bwa
  file(gridss_files) from reference.gridss

  output:
  file('*') into completegridss
  tuple val(germlineID), file("tumords.txt"), file("${params.runID}.output.vcf.gz") into gridssfilter

  when:
  params.seqlevel == "wgs"

  script:
  def jvmheap_taskmem = task.memory == null ? "" : "--jvmheap " + javaTaskmem("${task.memory}")
  def fasta = "${bwa}/*fasta"
  def gridss_blacklist = "${gridss_files}/gridss_blacklist.noChr.bed"
  def gridss_props = "${gridss_files}/dbs/gridss/gridss.properties"
  """
  GERMLINEBAM=\$(ls | grep ${germlineID} | grep bam\$ | grep -v bai)
  BAMFILES=\$(echo -n \$GERMLINEBAM" "\$(ls *.bam | grep -v \$GERMLINEBAM))
  LABELS=\$(echo -n ${germlineID}" "\$(ls *bam | grep -v ${germlineID} | grep -v assembly | cut -d "." -f1) | sed 's/\\s */,/g')
  TUMORDS=\$(echo \$LABELS | perl -ane '@s=split(/\\,/);for(\$i=2;\$i<=@s;\$i++){push(@o,\$i);} print join(",",@o[0..\$#o]) . "\\n";')
  TASKCPUS=\$(( ${task.cpus} / 4 )) ##"preprocessing will use up to 200-300% CPU per thread"
  echo \$TUMORDS > tumords.txt

  gridss.sh \
    --reference \$(echo ${fasta}) \
    --output ${params.runID}".output.vcf.gz" \
    --assembly ${params.runID}".assembly.bam" \
    --threads \$TASKCPUS \
    --jar /opt/gridss/gridss-2.9.4-gridss-jar-with-dependencies.jar \
    --workingdir ./ ${jvmheap_taskmem} \
    --blacklist ${gridss_blacklist} \
    --steps All \
    --configuration ${gridss_props} \
    --maxcoverage 50000 \
    --labels \$LABELS \
    \$BAMFILES
  """
}

/* 2.1.2: GRIDSS SV filtering
*/
process gridss_filter {

  label 'max_mem'
  publishDir path: "${params.outDir}/combined/gridss", mode: "copy"

  input:
  tuple val(germlineID), file(tumords), file("${params.runID}.output.vcf.gz") from gridssfilter

  output:
  file('*') into gridssfilterd
  tuple val(germlineID), file("${params.runID}.somatic_filter.vcf.bgz") into gridsspp

  when:
  params.seqlevel == "wgs"

  script:
  """
  Rscript --vanilla /opt/gridss/gridss_somatic_filter.R \
    --input ${params.runID}".output.vcf.gz" \
    --output ${params.runID}".somatic_filter.vcf" \
    --plotdir ./ \
    --scriptdir /opt/gridss \
    --normalordinal 1 \
    --tumourordinal \$(cat $tumords)
  """
}

//2.1.3 GRIDSS parse and plot
process gridss_vcf_pp {

  label 'low_mem'
  errorStrategy 'retry'
  maxRetries 3
  publishDir path: "${params.outDir}/combined/gridss", mode: "copy", pattern: "*.[pdf, tsv, png, vcf.gz]"

  input:
  tuple val(germlineID), file(vcf) from gridsspp
  file(bwa) from reference.bwa

  output:
  file('*') into completegridsspp
  file("*.pdf") into sendmail_gridss_pdf
  file("*.tsv") into sendmail_gridss_tsv

  when:
  params.seqlevel == "wgs"

  script:
  def dict = "${bwa}/*dict"
  def which_genome = params.assembly == "GRCh37" ? "hg19" : "hg38"
  """
  tabix ${vcf}
  Rscript -e "somenone::gridss_parse_plot(vcf = \\"${params.runID}.somatic_filter.vcf.bgz\\", germline_id = \\"${germlineID}\\", dict_file = \$(echo \\"${dict}\\"), which_genome = \\"${which_genome}\\", output_path = NULL)"
  """
}

// 2.2: HaplotypeCaller merge
process hc_merge {

  label 'high_mem'
  publishDir path: "${params.outDir}/samples/${sampleID}/haplotypecaller", mode: "copy", pattern: '*.vcf.*'

  input:
  tuple val(sampleID), val(meta), file(rawvcfs) from hc_fm

  output:
  tuple val(sampleID), val(meta), file("${sampleID}.hc.merge.vcf.gz"), file("${sampleID}.hc.merge.vcf.gz.tbi") into ( cpsr_vcf, vep_hc_vcf )

  script:
  """
  ls *.sort.hc.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".hc.merge.vcf"
  bgzip ${sampleID}".hc.merge.vcf"
  tabix ${sampleID}".hc.merge.vcf.gz"
  """
}

// 2.3: CPSR annotation of GATK4 Germline
process cpsrreport {

  label 'med_mem'

  publishDir "${params.outDir}/reports/cpsr", mode: "copy", pattern: "*.html"
  publishDir "${params.outDir}/samples/${sampleID}/cpsr", mode: "copy", pattern: "*[!.html]"

  input:
  tuple val(sampleID), val(meta), file(vcf), file(tbi) from cpsr_vcf
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  file('*') into cpsr_vcfs
  file("${metaid}.cpsr.${grchv}.html") into sendmail_cpsr

  script:
  grchv = "${grchver}".split("\\/")[-1]
  metaid = "${meta}".replaceAll("\\s *", "_").replaceAll("[\\[\\(\\)\\]]","").replaceAll("\"","")
  """
  {
  ##CPSR v0.6.1
  cpsr.py \
    --no-docker \
    --no_vcf_validate \
    --panel_id 0 \
    --query_vcf ${vcf} \
    --pcgr_dir ${pcgrbase} \
    --output_dir ./ \
    --genome_assembly ${grchv} \
    --conf ${pcgrbase}/data/${grchv}/cpsr_configuration_default.toml \
    --sample_id ${metaid}
  } 2>&1 | tee > ${sampleID}.cpsr.log.txt
  """
}

// 2.4: PicardTools metrics suite for MultiQC HTML report
process mltmet {

  label 'med_mem'

  publishDir "${params.outDir}/samples/${sampleID}/metrics", mode: "copy"

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from gmultimetricing
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(intlist) from reference.intlist

  output:
  file('*.txt') into multimetrics_multiqc

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  if [[ ${params.seqlevel} == "exome" ]]; then
  picard ${taskmem} CollectHsMetrics \
    I=${bam} \
    O=${sampleID}".hs_metrics.txt" \
    TMP_DIR=./ \
    R=${fasta} \
    VALIDATION_STRINGENCY=LENIENT \
    BAIT_INTERVALS=${intlist}  \
    TARGET_INTERVALS=${intlist}
  fi
  picard ${taskmem} CollectAlignmentSummaryMetrics \
    I=${bam} \
    O=${sampleID}".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    VALIDATION_STRINGENCY=LENIENT \
    R=${fasta}

  picard ${taskmem} CollectMultipleMetrics \
    I=${bam} \
    O=${sampleID}".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    VALIDATION_STRINGENCY=LENIENT \
    R=${fasta}

  picard ${taskmem} CollectSequencingArtifactMetrics \
    I=${bam} \
    O=${sampleID}".artifact_metrics.txt" \
    TMP_DIR=./ \
    VALIDATION_STRINGENCY=LENIENT \
    R=${fasta}

  picard ${taskmem} CollectInsertSizeMetrics \
    I=${bam} \
    O=${sampleID}".insert_size_metrics.txt" \
    H=${bam}".histogram.pdf" \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=./

  } 2>&1 | tee > ${sampleID}.picard.metrics.log
  """
}

// 2.5: filter germline channel, tap into somatic channels for all processes subsequent
if(!params.germOnly){

  def germfilter = branchCriteria {
                    germfiltered: it[0] == "germline"
                    return it

                    somafiltered: true
                    return it
                   }
  germfiltering
      .branch(germfilter)
      .set { somagerm }

  somagerm.somafiltered
      .map { [it[1], it[2..-1]] }
      .tap { somatap }

  somagerm.germfiltered
      .map { [it[1], it[2..-1]] }
      .tap { germtap }
  somatap.combine(germtap).tap{ somagermtap }

  somagermtap
    .map { it -> tuple(it[0],
                       it[1][0..1],
                       it[2],
                       it[3][0..1]).flatten() }
    .into { mutect2somaticing; mutect2_contam; facetsomaing; mantastrelka2ing; lanceting; gpling }

  //2.6.1: SCNA with facets CSV snp-pileup
  process fctcsv {

    label 'med_mem'

    publishDir "${params.outDir}/samples/${sampleID}/facets", mode: "copy"

    input:
    tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from facetsomaing
    file(dbsnp_files) from reference.dbsnp

    output:
    tuple file("${sampleID}.fit_cncf_jointsegs.tsv"), file("${sampleID}.fit_ploidy_purity.tsv") into ( facets_consensusing, facets_pyclone )
    tuple val(sampleID), file("${sampleID}.cncf_jointsegs.pcgr.tsv"), file("${sampleID}.fit_ploidy_purity.pcgr.tsv") into facets_pcgr
    file("${sampleID}.facets.log.txt") into facets_log

    script:
    def dbsnp = "${dbsnp_files}/*gz"
    """
    { snp-pileup \
        \$(echo ${dbsnp}) \
        -r 10 \
        -p \
        ${sampleID}.facets.r10.csv \
        ${germlinebam} \
        ${tumourbam}

      Rscript -e "somenone::facets_cna_call(\\"${sampleID}.facets.r10.csv\\")"

      tail -n+2 ${sampleID}.fit_ploidy_purity.tsv > ${sampleID}.fit_ploidy_purity.pcgr.tsv
    } 2>&1 | tee > ${sampleID}.facets.log.txt
    """
  }

  // 2.6.2: SCNA consensus from facets
  process fctcon {
    label 'med_mem'

    publishDir "${params.outDir}/combined/facets", mode: "copy"

    input:
    file(filesn) from facets_consensusing.collect()
    file(cosmicbed) from reference.cosmic
    file(dict) from reference.dict

    output:
    file('*') into complete_facets
    file(filesn) into pairtee_facets
    file("${params.runID}.ENS.facets.CNA.master.tsv") into pairtree_facet
    file('*.pdf') into sendmail_facets

    script:
    if( !params.cosmic )
      """
      { Rscript -e "somenone::facets_cna_consensus(\\"fit_cncf_jointsegs.tsv\\", \\"${dict}\\", \\"${params.runID}\\")"
      } 2>&1 | tee > facets_cons.log.txt
      """
    else
      """
      { Rscript -e "somenone::facets_cna_consensus(\\"fit_cncf_jointsegs.tsv\\", \\"${dict}\\", \\"${params.runID}\\", \\"${cosmicbed}\\")"
      } 2>&1 | tee > facets_cons.log.txt
      """
  }

  mutect2bedding = mutect2_bedding.flatten()
  mutect2somaticing
    .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5]]}
    .combine(mutect2bedding)
    .set { mutect2somaticbedding }

  // 2.7.1: MuTect2
  // NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error

  process mutct2_sg {

    label 'med_mem'

    input:
    tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(intlist) from mutect2somaticbedding
    file(fasta) from reference.fa
    file(fai) from reference.fai
    file(dict) from reference.dict

    output:
    tuple val(sampleID), file('*sort.mutect2.vcf') into mutect2_gt
    tuple val(sampleID), file('*.vcf.stats') into mutect2_st
    tuple val(sampleID), file('*mutect2.f1r2.tar.gz') into mutect2_f1r2

    script:
    def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
    """
    SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print\$s[2];')
    gatk ${taskmem} \
      Mutect2 \
      --native-pair-hmm-threads ${task.cpus} \
      --reference ${fasta} \
      --input ${germlinebam} \
      --input ${tumourbam} \
      --normal-sample ${germlineID} \
      --tumor-sample ${sampleID} \
      --output ${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
      --disable-sequence-dictionary-validation true \
      --f1r2-tar-gz \${SCATGATHN}".mutect2.f1r2.tar.gz" \
      -L ${intlist}

    picard SortVcf \
      I=${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
      O=${sampleID}"."\${SCATGATHN}".sort.mutect2.vcf" \
      SD=${dict}
    """
  }

  // 2.7.2: MuTect2_merge
  mutect2_gt
    .groupTuple()
    .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
    .set { mutect2_fm }

  process mutct2_concat {

    label 'med_mem'

    input:
    tuple val(sampleID), file(rawvcfs) from mutect2_fm

    output:
    tuple val(sampleID), file('*mutect2.merge.vcf') into mutect2_merge

    script:
    """
    ls *.sort.mutect2.vcf > vcf.list
    picard MergeVcfs I=vcf.list O=${sampleID}".mutect2.merge.vcf"
    """
  }

  mutect2_st
    .groupTuple()
    .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
    .set { mutect2_sm }

  // 2.7.3: MuTect2 Concatenate VCFs
  process mutct2_concstat {

    label 'med_mem'

    input:
    tuple val(sampleID), file(stats) from mutect2_sm

    output:
    tuple val(sampleID), file('*mutect2.merge.vcf.stats') into mutect2_stats

    script:
    """
    STATS=\$(ls *stats | perl -ane 'foreach \$k (@F){print "--stats \$k ";}')
    gatk MergeMutectStats --output ${sampleID}".mutect2.merge.vcf.stats" \$STATS
    """
  }

  // 2.7.4: MuTect2 Concatenate VCFs
  mutect2_f1r2.groupTuple()
              .map { it -> [it[0], it[1..-1].flatten()] }
              .set { mutect2_f1r2_set }

  process mutct2_f1r2_comb {

    label 'med_mem'

    input:
    tuple val(sampleID), file(mutect2_ro) from mutect2_f1r2_set

    output:
    tuple val(sampleID), file("${sampleID}.mutect2.f1r2.tar.gz") into mutect2_f1r2_comb

    script:
    """
    ALL_F1R2_INPUT=\$(for x in *.mutect2.f1r2.tar.gz; do echo -n "-I \$x "; done)
    gatk LearnReadOrientationModel \$ALL_F1R2_INPUT -O ${sampleID}.mutect2.f1r2.tar.gz
    """
  }

  // 2.7.5: MuTect2 Contamination
  mutect2_contam
    .join(mutect2_merge)
    .join(mutect2_stats)
    .join(mutect2_f1r2_comb)
    .groupTuple()
    .map { it -> [it[0], it[1..5].flatten(), it[6], it[7], it[8]].flatten() }
    .set { mutect2_contam_merge }

  process mutct2_contam_filter {

    label 'med_mem'

    publishDir path: "${params.outDir}/samples/${sampleID}/mutect2", mode: "copy", overwrite: true

    input:
    tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(mergevcf), file(statsvcf), file(readorient) from mutect2_contam_merge
    file(fasta) from reference.fa
    file(fai) from reference.fai
    file(dict) from reference.dict
    file(gps_files) from reference.seqlevel
    file(intlist) from reference.intlist

    output:
    file('*snv_indel.pass.vcf') into mutect2_veping mode flatten
    file('*.raw.vcf') into mutect2_rawVcf
    file('*') into completedmutect2call

    script:
    def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
    hg = params.assembly == "GRCh37" ? "hg19" : "hg38"
    gpsgz = params.seqlevel == "exome" ? "${gps_files}/${params.exomeTag}/af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf.gz" : "${gps_files}/af-only-gnomad.wgs.${hg}.noChr.vcf.gz"
    """
    gatk ${taskmem} \
      GetPileupSummaries \
      -I ${tumourbam} \
      -V ${gpsgz} \
      -O ${sampleID}".getpileupsummaries.table" \
      -L ${intlist}

    gatk CalculateContamination \
      -I ${sampleID}".getpileupsummaries.table" \
      -O ${sampleID}".calculatecontamination.table"

    CONTAM=\$(tail -n+2 ${sampleID}.calculatecontamination.table | cut -f 2 | cut -d "." -f 1)
    if [[ \$CONTAM != 0 ]]; then
      touch ${sampleID}".CONTAMINATION.WARNING.txt"
    fi

    gatk IndexFeatureFile \
      --input ${mergevcf}

    gatk ${taskmem} \
      FilterMutectCalls \
      --reference ${fasta} \
      --contamination-table ${sampleID}".calculatecontamination.table" \
      --interval-padding 5 \
      --output ${sampleID}".mutect2.FilterMutectCalls.vcf" \
      --unique-alt-read-count 3 \
      --variant ${mergevcf} \
      --stats ${statsvcf} \
      --disable-sequence-dictionary-validation true \
      --ob-priors ${readorient} \
      -L ${intlist}

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=${sampleID} \
      DP=14 \
      MD=2 \
      VCF=${sampleID}".mutect2.FilterMutectCalls.vcf"
    """
  }

  // 2.9: Manta output is a pre-req for Strelka2, so call both here
  process mntstr {

    label 'high_mem'

    publishDir path: "${params.outDir}/samples/${sampleID}/manta-strelka2", mode: "copy"

    input:
    tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
    file(fasta) from reference.fa
    file(fai) from reference.fai
    file(dict) from reference.dict
    file(bed_files) from reference.seqlevel

    output:
    file("${sampleID}.strelka2.snv_indel.pass.vcf") into strelka2_veping
    file("${sampleID}.strelka2.raw.vcf") into strelka2_rawVcf
    file('*.txt') into log_mantastrelka

    script:
    def bedgz = params.seqlevel == "wgs" ? "${bed_files}/wgs.bed.gz" : "${bed_files}/${params.exomeTag}/${params.exomeTag}.bed.gz"
    def callRegions = params.seqlevel == "exome" ? "--exome --callRegions ${bedgz}" : "--callRegions ${bedgz}"
    """
    {
      configManta.py ${callRegions} --referenceFasta=${fasta} --normalBam=${germlinebam} --tumourBam=${tumourbam} --runDir=manta

      manta/runWorkflow.py -m local -j ${task.cpus}

      configureStrelkaSomaticWorkflow.py ${callRegions} --referenceFasta=${fasta} --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz --normalBam=${germlinebam} --tumorBam=${tumourbam} --runDir=strelka2

      strelka2/runWorkflow.py -m local -j ${task.cpus}

      ##merge into raw snv_indel
      gatk MergeVcfs -I strelka2/results/variants/somatic.snvs.vcf.gz -I strelka2/results/variants/somatic.indels.vcf.gz -O tmp.strelka2.snv_indel.vcf

      ${workflow.projectDir}/bin/manta_strelka2_rename_filter.sh  tmp.strelka2.snv_indel.vcf tmp2.strelka2.snv_indel.vcf ${sampleID} ${germlineID}

      perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
          ID=${sampleID} \
          DP=14 \
          MD=2 \
          VCF=tmp2.strelka2.snv_indel.vcf

    } 2>&1 | tee > ${sampleID}.manta-strelka2.log.txt
    """
  }

  // 2.10.1: Lancet
  lancetbedding = lancet_bedding.flatten()
  lanceting
    .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5]]}
    .combine(lancetbedding)
    .set { lancetsbedding }

  process lancet_sg {

    label 'med_mem'

    input:
    tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(bed) from lancetsbedding
    file(fasta) from reference.fa
    file(fai) from reference.fai
    file(dict) from reference.dict

    output:
    tuple val(sampleID), file('*.sort.lancet.vcf') into lancet_gt

    when:
    params.seqlevel == "exome"

    script:
    scatgathn = "${bed}".split("\\.")[2]
    """
    lancet \
      --num-threads ${task.cpus} \
      --ref ${fasta} \
      --bed ${bed} \
      --tumor ${tumourbam} \
      --normal ${germlinebam} | \
      perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
        \$_=~s/TUMOR/${sampleID}/;
        \$_=~s/NORMAL/${germlineID}/;
        print \$_;}
      else{print \$_;}' > ${sampleID}"."${scatgathn}".lancet.vcf"

    picard SortVcf \
      I=${sampleID}"."${scatgathn}".lancet.vcf" \
      O=${sampleID}"."${scatgathn}".sort.lancet.vcf" \
      SD=${dict}
    """
  }

  // 2.10.2: Lancet Merge
  lancet_gt
    .groupTuple()
    .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
    .set { lancet_fm }

  process lancet_concat {

    label 'med_mem'

    input:
    tuple val(sampleID), file(rawvcf) from lancet_fm

    output:
    tuple val(sampleID), file('*lancet.merge.vcf') into lancet_merge

    script:
    """
    ls *.sort.lancet.vcf > vcf.list
    picard MergeVcfs I=vcf.list O=${sampleID}".lancet.merge.vcf"
    """
  }

  /* 2.10.2: Lancet Filter
  */
  process lancet_filter {

    label 'med_mem'

    publishDir path: "${params.outDir}/samples/${sampleID}/lancet", mode: "copy"

    input:
    tuple val(sampleID), file(mergevcf) from lancet_merge

    output:
    file('*.pass.vcf') into lancet_veping mode flatten
    file('*.raw.vcf') into lancet_rawVcf
    file('*') into completedlancetcall

    script:
    """
    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=${sampleID} \
      DP=14 \
      MD=2 \
      VCF=${mergevcf}
    """
  }

  /*
  ================================================================================
                            3.  ANNOTATION AND REPORTING
  ================================================================================
  */

  // 3.01: HC_merge VEP
  process vepHC {

    label 'low_mem'

    publishDir path: "${params.outDir}/samples/${sampleID}/haplotypecaller", mode: "copy"

    input:
    tuple val(sampleID), val(meta), file(vcf), file(tbi) from vep_hc_vcf
    file(fasta) from reference.fa
    file(fai) from reference.fai
    file(dict) from reference.dict
    file(grchver) from reference.grchvers
    file(pcgrbase) from reference.pcgrbase

    output:
    file("${sampleID}.hc.merge.vep.vcf") into hc_vepd

    script:
    grch_vers = "${grchver}".split("\\/")[-1]
    """
    vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
      --offline \
      --assembly ${params.assembly} \
      --vcf_info_field ANN \
      --symbol \
      --species homo_sapiens \
      --check_existing \
      --cache \
      --fork ${task.cpus} \
      --af_1kg \
      --af_gnomad \
      --vcf \
      --input_file ${vcf} \
      --output_file ${sampleID}.hc.merge.vep.vcf \
      --format "vcf" \
      --fasta ${fasta} \
      --hgvs \
      --canonical \
      --ccds \
      --force_overwrite \
      --verbose
    """
  }

  // 3.02: Annotate Vcfs
  ALLVCFS = lancet_veping
            .mix( mutect2_veping )
            .mix( strelka2_veping )

  process vepann {

    label 'med_mem'

    publishDir path: "${params.outDir}/samples/${sampleID}/${caller}", mode: "copy", pattern: "${sampleID}.${caller}.snv_indel.pass.vep.vcf"

    input:
    each file(vcf) from ALLVCFS
    file(fasta) from reference.fa
    file(fai) from reference.fai
    file(dict) from reference.dict
    file(grchver) from reference.grchvers
    file(pcgrbase) from reference.pcgrbase

    output:
    file('*.vcf') into runGRanges

    script:
    def grch_vers = "${grchver}".split("\\/")[-1]
    def vcf_anno = "${vcf}".replaceAll(".vcf", ".vep.vcf")
    sampleID = "${vcf}".split("\\.")[0]
    caller = "${vcf}".split("\\.")[1]
    """
    vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
      --offline \
      --assembly ${params.assembly} \
      --vcf_info_field ANN \
      --symbol \
      --species homo_sapiens \
      --check_existing \
      --cache \
      --fork ${task.cpus} \
      --af_1kg \
      --af_gnomad \
      --vcf \
      --input_file ${vcf} \
      --output_file ${vcf_anno} \
      --format "vcf" \
      --fasta ${fasta} \
      --hgvs \
      --canonical \
      --ccds \
      --force_overwrite \
      --verbose
    """
  }

  /* 3.1 RData GRanges from processed VCFs
  * take publishDir and check for number of files therein
  * each sample has 6 associated (raw, pass per caller)
  * NB increment if adding callers!
  */
  runGRanges
   .mix(lancet_rawVcf)
   .mix(mutect2_rawVcf)
   .mix(strelka2_rawVcf)
   .collect()
   .set { allvcfs }

  //separate out impacts processing as WGS cango above walltime
  impacts = ["HIGH,MODERATE,MODIFIER,LOW"]

  process vcfGRa {

    label 'max_mem'

    publishDir "${params.outDir}/combined/variant_consensus", mode: "copy"

    input:
    file(grangesvcfs) from allvcfs
    val(germlineID) from vcfGRaID.unique()
    each impact from impacts

    output:
    file('*impacts.pcgr.tsv.vcf') into vcfs_pcgr
    file('*') into completedvcfGRangesConsensus
    file('*.pdf') into sendmail_vcfGRa

    script:
    def inc_ord = params.incOrder ? params.incOrder : "noord"
    def which_genome = params.assembly == "GRCh37" ? "hg19" : "hg38"
    """
    Rscript -e "somenone::variant_consensus(germline_id = \\"${germlineID}\\", vep_vcf_pattern = \\"snv_indel.pass.vep.vcf\\", raw_vcf_pattern = \\"raw.vcf\\", tag = \\"${params.runID}\\", which_genome = \\"${which_genome}\\", included_order = \\"${inc_ord}\\", impacts = \\"${impact}\\")"
    """
  }

  // 3.2 Create VCF for PCGR from consensus (reheader based on vcf 4.2 standards)
  vcfs_pcgrd = vcfs_pcgr
                .collect()
                .flatten()

  process pcgr_vcf {

    label 'low_mem'

    input:
    file(vcf) from vcfs_pcgrd

    output:
    tuple val(sampleID), file(ovcf) into snvpass_pcgr

    when:
    vcf =~ "HMML_impacts.pcgr.tsv.vcf"

    script:
    sampleID = "${vcf}".split("\\.")[0]
    ovcf = "${vcf}".replace("pcgr.tsv.vcf", "snv_indel.pass.pcgr.vcf")
    """
    cat ${workflow.projectDir}/assets/vcf42.head.txt > $ovcf
    head -n1 $vcf >> $ovcf
    tail -n+2 $vcf | sort -V >> $ovcf
    """
  }

  snvpass_pcgr
    .join(metas_pcgr)
    .join(facets_pcgr)
    .map { it -> tuple(it[0], it[1..-1]).flatten() }
    .set { pcgr_inputs }

  // 3.3 PCGR report
  // take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR
  process pcgrreport {

    label 'low_mem'
    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.outDir}/reports/pcgr", mode: "copy", pattern: "*html"
    publishDir "${params.outDir}/samples/${sampleID}/pcgr", mode: "copy", pattern: "*[!.html]"

    input:
    tuple val(sampleID), file(vcf), val(meta), file(jointsegs), file(ploidpur) from pcgr_inputs
    file(grchver) from reference.grchvers
    file(pcgrbase) from reference.pcgrbase
    file(exomebase) from reference.seqlevel

    output:
    file('*') into completed_pcgr
    file("*.pcgr_acmg.${grch_vers}.html") into sendmail_pcgr

    script:
    grch_vers = "${grchver}".split("\\/")[-1]
    config = params.seqlevel == "exome" ? "${exomebase}/${params.exomeTag}/pcgr_configuration_${params.exomeTag}.toml" : "${pcgrbase}/data/${grch_vers}/pcgr_configuration_default.toml"
    assay = params.seqlevel == "exome" ? "WES" : "WGS"
    """
    {
    ##want META to allow spaces, remove non-alphanum
    META=\$(echo ${meta} | perl -ane '\$_=~s/[^a-zA-Z0-9_ \\n]//g; print \$_;' | sed 's/\\s */_/g')

    PLOIDY=""; PURITY="";
    if [[ \$(cut -f 1 ${ploidpur}) != "NA" ]]; then
      PLOIDY="--tumor_ploidy \$(cut -f 1 ${ploidpur})"
    fi
    if [[ \$(cut -f 2 ${ploidpur}) != "NA" ]]; then
      PURITY="--tumor_purity \$(cut -f 2 ${ploidpur})"
    fi

    ##PCGR 0.9.1
    pcgr.py \
      --pcgr_dir ${pcgrbase} \
      --output_dir ./ \
      --genome_assembly ${grch_vers} \
      --conf ${config} \
      --sample_id \$META \
      --input_vcf ${vcf} \
      --input_cna ${jointsegs} \$PLOIDY \$PURITY \
      --no-docker \
      --force_overwrite \
      --no_vcf_validate \
      --estimate_tmb \
      --estimate_msi_status \
      --estimate_signatures \
      --include_trials \
      --assay ${assay}

    } 2>&1 | tee > ${sampleID}.pcgr.log.txt
    """
  }

  if(params.phylogeny){
    //3.41 PycloneVI / Pairtree
    process pairtree_check {

      input:
      file(rdata) from completedvcfGRangesConsensus

      output:
      file("${params.runID}.HMML_impacts.master_consensus_all.RData") into pairtree_rdata

      when:
      params.phylogeny == true

      script:
      """
      """
    }

    process pairtree_setup {

      label 'low_mem'

      publishDir "${params.outDir}/combined/pylogeny", mode: "copy"

      input:
      file(rdata) from pairtree_rdata
      file(filesn) from pairtee_facets
      file(cna_master) from pairtree_facet

      output:
      tuple file("${params.runID}.pairtree.psm"), file("${params.runID}.in_params_fn.json") into pairtree_in

      when:
      params.phylogeny == true

      script:
      def which_genome = params.assembly == "GRCh37" ? "hg19" : "hg38"
      """
      Rscript -e "somenone::make_pairtree_input(rdata_input = \\"${rdata}\\", cn_master = \\"${cna_master}\\", cn_pattern = \\"fit_cncf_jointsegs.tsv\\", pp_pattern = \\"fit_ploidy_purity.tsv\\", which_genome = \\"${which_genome}\\", tag = \\"${params.runID}\\")"
      """
    }

    //3.42
    // process pyclonevi {
    //
    //   label 'low_mem'
    //
    //   publishDir "${params.outDir}/combined/pylogeny", mode: "copy"
    //
    //   input:
    //   tuple file(pyclone_input), file(pairtree_psm) from pyclone_in
    //
    //   output:
    //   tuple file("${params.runID}.pyclone.results.tsv"), file("${pairtree.psm}") into pyclone_res
    //
    //   script:
    //   """
    //   pyclone-vi fit -i ${pyclone_input} \
    //                  -o ${params.runID}.pyclonevi.output.tsv
    //   pyclone-vi write-results-file -i ${params.runID}.pyclonevi.output.tsv \
    //                                 -o ${params.runID}.pyclone.results.tsv
    //   """
    // }

    //3.42: pairtree run
    process pairtree_run {

      label 'low_mem'

      publishDir "${params.outDir}/combined/phylogeny/pairtree/${model}_${concn}", mode: "copy"

      input:
      tuple file(pairtree_psm), file(pairtree_json) from pairtree_in
      each concn from Channel.from("-2","-1","0.5","1.5")
      each model from Channel.from("pairwise","linfreq")

      output:
      file('*') into pairtree_res
      file('*.html') into sendmail_pairtree

      when:
      params.phylogeny == true

      script:
      def concnt = "${concn}" < 0 ? "${concn}".replaceAll("-", "minus") : "${concn}"
      """
      cut -f 1,2,3,4,5 ${pairtree_psm} > ${params.runID}.pairtree_${model}_${concnt}.ssm

      clustervars --model ${model} \
                  --concentration ${concn} \
                  ${params.runID}.pairtree_${model}_${concnt}.ssm \
                  ${pairtree_json} \
                  ${params.runID}.out_params_${model}_${concnt}.json

      touch ${params.runID}.rmvaf_params_${model}_${concnt}.json
      python /opt/miniconda/envs/pairtree/share/pairtree/util/remove_high_vaf.py \
             ${params.runID}.pairtree_${model}_${concnt}.ssm \
             ${params.runID}.out_params_${model}_${concnt}.json \
             ${params.runID}.rmvaf_params_${model}_${concnt}.json

      WCLTEST=\$(wc -l ${params.runID}.rmvaf_params_${model}_${concnt}.json | perl -ane 'print \$F[0];')
      if [[ \$WCLTEST < 1 ]]; then
        OUTPARAMS=${params.runID}.out_params_${model}_${concnt}.json
      else
        OUTPARAMS=${params.runID}.rmvaf_params_${model}_${concnt}.json
      fi

      pairtree --params \$OUTPARAMS \
               ${params.runID}.pairtree_${model}_${concnt}.ssm \
               ${params.runID}.res_${model}_${concnt}.npz

      plottree ${params.runID}.pairtree_${model}_${concnt}.ssm \
               \$OUTPARAMS \
               ${params.runID}.res_${model}_${concnt}.npz \
               ${params.runID}.pairtree_${model}_${concnt}.results.html
      """
    }

    sendmail_soma
      .mix(sendmail_pairtree)
      .set { sendmail_soma }
  }

  //3.5: Pathseq
  if(!params.bamCsv){
    if(params.microbiome){
      process Pathseq {

        label 'high_mem'

        publishDir "${params.outDir}/samples/${sampleID}/pathseq", mode: "copy"

        input:
        tuple val(type), val(sampleID), val(meta), file(ubam), file(ubai) from pathseqing
        file(pathseq_refs) from reference.pathseq

        output:
        file('*') into  pathseq_res
        file('*.txt') into sendmail_pathseq
        when:
        params.microbiome == true

        script:
        def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""

        """
         gatk ${taskmem} PathSeqPipelineSpark  \
           --input ${ubam} \
           --kmer-file ${pathseq_refs}/pathseq_host.bfi \
           --filter-bwa-image ${pathseq_refs}/pathseq_host.fa.img \
           --microbe-bwa-image ${pathseq_refs}/pathseq_microbe.fa.img \
           --microbe-dict ${pathseq_refs}/pathseq_microbe.dict \
           --taxonomy-file ${pathseq_refs}/pathseq_microbe_taxonomy.db \
           --min-clipped-read-length 60 \
           --min-score-identity 0.90 \
           --identity-margin 0.02 \
           --scores-output ${sampleID}.pathseq.scores.txt \
           --output ${sampleID}.pathseq.output_reads.bam \
           --filter-metrics ${sampleID}.pathseq.filter_metrics.txt \
           --score-metrics ${sampleID}.pathseq.score_metrics.txt
        """
      }

    sendmail_soma
      .mix(sendmail_pathseq)
      .set { sendmail_soma }
    }
  }

  fastp_multiqc.collect()
    .mix(fastqc_multiqc.collect())
    .mix(mrkdup_multiqc.collect())
    .set{ mrkdup_multiqc_col }

  sendmail_pcgr
    .mix(sendmail_vcfGRa)
    .mix(sendmail_facets)
    .set { sendmail_soma }
} else {
  //germOnly
  mrkdup_multiqc
    .collect()
    .set { mrkdup_multiqc_col }
}

/*
================================================================================
                          4.  MULTIQC AND CLOSEOUT
================================================================================
*/
// 4.0 Run multiQC to finalise report
process MultiQC {

  label 'low_mem'
  publishDir path: "${params.outDir}/reports/multiQC", mode: "copy"

  input:
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()
  file(mrkdups) from mrkdup_multiqc_col
  file(mosdepth) from mosdepth_multiqc.collect()

  output:
  file('*') into completedmultiqc
  file("*.html") into sendmail_multiqc

  script:
  """
  multiqc . -i ${params.runID}".somatic_n-of-1" --tag DNA -f -c ${params.multiqcConfig}
  """
}

// 4.1.1: somatic_n-of-1 container software versions
process somenone_software_vers {

  label 'low_mem'
  publishDir "pipeline_info", mode: 'copy'
  publishDir path: "${params.outDir}/reports/software_vers", mode: "copy"

  output:
  file 'somenone_software_versions.yaml' into ch_somenone_software_vers

  script:
  """
  conda env export > somenone_software_versions.yaml
  """
}

//4.1.2: gridss software version numbers
//out of use as using Docker from GRIDSS now
process gridss_software_vers {

  label 'low_mem'
  publishDir "${params.outDir}/pipeline_info", mode: 'copy'

  output:
  file 'gridss_software_versions.txt' into ch_gridss_software_vers

  when:
  params.seqlevel == "wgs"

  script:
  """
  ls -l /opt/gridss > gridss_software_versions.txt
  """
}

//4.1.3: pcgr software version numbers
//out of use as using Docker from pcgr now
process pcgr_software_vers {

  label 'low_mem'
  publishDir "${params.outDir}/pipeline_info", mode: 'copy'

  output:
  file 'pcgr_software_versions.txt' into ch_pcgr_software_vers


  script:
  """
  pcgr.py --version > pcgr_software_versions.txt
  cpsr.py --version >> pcgr_software_versions.txt
  """
}

// 4.19: ZIP for sending on sendmail
sendmail_cpsr
  .mix(sendmail_multiqc)
  .mix(sendmail_cnvkit)
  .set { sendmail_germ }

if(!params.germOnly) {
  sendmail_germ
    .mix(sendmail_soma)
    .mix(sendmail_gridss_pdf)
    .mix(sendmail_gridss_tsv)
    .set { sendmail_all }
}
if(params.germOnly) {
  sendmail_germ
  .set { sendmail_all }
}

process zipup {

  label 'low_mem'
  publishDir "${params.outDir}/reports", mode: 'copy'

  input:
  file(send_all) from sendmail_all.collect()

  output:
  file("${params.runID}.somatic_n-of-1.zip") into send_zip

  script:
  """
  mkdir html_reports && mv *html ./html_reports/
  mkdir other && mv *.* ./other/
  zip -r ${params.runID}.somatic_n-of-1.zip *
  """
}

// 4.2: Completion e-mail notification

workflow.onComplete {
  sleep(1000)
  def subject = """\
    [brucemoran/somatic_n-of-1] SUCCESS: $params.runID [$workflow.runName]
    """
    .stripIndent()
  if (!workflow.success) {
      subject = """\
        [brucemoran/somatic_n-of-1] FAILURE: $params.runID [$workflow.runName]
        """
        .stripIndent()
  }

  def msg = """\
    Pipeline execution summary
    ---------------------------
    RunID       : ${params.runID}
    RunName     : ${workflow.runName}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
    .stripIndent()

  def attachments = send_zip.toList().getVal()

  sendMail(to: "${params.email}",
           subject: subject,
           body: msg,
           attach: attachments)
}
