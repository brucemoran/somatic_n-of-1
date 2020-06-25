#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '-----------------------------------------------------------------------'
  log.info 'NEXTFLOW 19.10 FASTQ QC, TRIM, ALIGN, SOMATIC+GERMLINE SNV, CNA, REPORT'
  log.info '-----------------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run brucemoran/somatic_n-of-1'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    -profile    Configuration profile (required: standard,singularity)'
  log.info '    --sampleCsv      STRING      CSV format, headers: type (either "germline" or "somatic"),sampleID,meta,/path/to/read1.fastq.gz,/path/to/read2.fastq.gz; use meta for naming in PCGR, CPSR report'
  log.info '    --runID      STRING      name for run, used to tag outputs'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --refDir        STRING      dir in which reference data and required indices are held; if not specified, this is created by brucemoran/DNAseq_references/GRCh_37-38 (default: work/GRCh38)'
  log.info '    --germline      STRING      run HaplotypeCaller on germline sample and annotate with CPSR (true/false, default: true)'
  log.info '    --scatGath      NUM         number of pieces to divide intervalList into for scattering to variant calling processes (default: 20 for exome, 100 for WGS)'
  log.info '    --incOrder      STRING      in final plots, use this ordering of samples (if multiple somatic samples); comma-separated, no spaces (default: alphanumeric sort)'
  log.info '    --multiqcConfig      STRING      config file for multiqc (default: bin/somatic_n-of-1.multiQC_config.yaml)'
  log.info '    --seqLevel      STRING      WGS or exome (default: WGS)'
  log.info ''
  exit 1
}

/* -2 Test if refDir is defined, if not run DNAseq_references pipeline under defaults
*/
if(!params.refDir){
  exit 1, "Please run: nextflow run brucemoran/DNAseq_references --outDir work -profile standard,singularity, then specify: nextflow run brucemoran/somatic_n-of-1 --refDir work/GRCh38"
}

/* -1: Global Variables
*/
params.outDir = "analysis/${params.seqLevel}"

// Reference data as params, and reusable therefore
params.fasta = Channel.fromPath("$params.refDir/*fasta").getVal()
params.fai = Channel.fromPath("$params.refDir/*fasta.fai").getVal()
params.dict = Channel.fromPath("$params.refDir/*dict").getVal()

params.amb = Channel.fromPath("$params.refDir/*fasta.amb").getVal()
params.ann = Channel.fromPath("$params.refDir/*fasta.ann").getVal()
params.bwt = Channel.fromPath("$params.refDir/*fasta.bwt").getVal()
params.pac = Channel.fromPath("$params.refDir/*fasta.pac").getVal()
params.sa = Channel.fromPath("$params.refDir/*fasta.sa").getVal()

params.twobit = Channel.fromPath("$params.refDir/*fasta.2bit").getVal()

params.seqlevel = "$params.seqLevel".toLowerCase()

params.intlist = Channel.fromPath("$params.refDir/${params.seqlevel}/*.bed.interval_list").getVal()
params.bed = Channel.fromPath("$params.refDir/${params.seqlevel}/*.bed").getVal()
params.bedgz = Channel.fromPath("$params.refDir/${params.seqlevel}/*.bed.gz").getVal()
params.bedgztbi = Channel.fromPath("$params.refDir/${params.seqlevel}/*.bed.gz.tbi").getVal()

params.dbsnp = Channel.fromPath("$params.refDir/dbsnp*.gz").getVal()
params.dbsnptbi = Channel.fromPath("$params.refDir/dbsnp*.tbi").getVal()

params.omni = Channel.fromPath("$params.refDir/KG_omni*.gz").getVal()
params.otbi = Channel.fromPath("$params.refDir/KG_omni*.gz.tbi").getVal()
params.kgp1 = Channel.fromPath("$params.refDir/KG_phase1*.gz").getVal()
params.ktbi = Channel.fromPath("$params.refDir/KG_phase1*.gz.tbi").getVal()
params.hpmp = Channel.fromPath("$params.refDir/hapmap*.gz").getVal()
params.htbi = Channel.fromPath("$params.refDir/hapmap*.gz.tbi").getVal()

params.gps = Channel.fromPath("$params.refDir/${params.seqlevel}/af-only-gnomad.*.vcf.gz").getVal()
params.gpstbi = Channel.fromPath("$params.refDir/${params.seqlevel}/af-only-gnomad.*.vcf.gz.tbi").getVal()

//PCGR, CPSR version and base data dir
params.grchvers = Channel.fromPath("$params.refDir/pcgr/data").getVal()
params.pcgrdir = Channel.fromPath("$params.refDir/pcgr").getVal()

//GRIDSS params
params.blacklist = Channel.fromPath("$params.refDir/gridss_blacklist.noChr.bed").getVal()
params.gridssProps = Channel.fromPath("$params.refDir/gridss.properties").getVal()
params.gridssponbedpe = Channel.fromPath("$params.refDir/gridss_pon_breakpoint.bedpe").getVal()
params.gridssponsinbed = Channel.fromPath("$params.refDir/gridss_pon_single_breakend.bed").getVal()

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

//somaticVariantConsensus scripts
process cons_scripts {

  executor 'local'

  publishDir "${workflow.projectDir}/bin"

  output:
  tuple file('variants_GRanges_consensus_plot.call.R'), file('variants_GRanges_consensus_plot.func.R') into vcfGRa_Scripts

  script:
  """
  wget "https://raw.githubusercontent.com/brucemoran/somaticVariantConsensus/master/scripts/variants_GRanges_consensus_plot.func.R"
  wget "https://raw.githubusercontent.com/brucemoran/somaticVariantConsensus/master/scripts/variants_GRanges_consensus_plot.call.R"
  """
}

/* 0.00: Input using sample.csv
*/
Channel.fromPath("$params.sampleCsv")
       .splitCsv( header: true )
       .map { row -> [row.type, row.sampleID, row.meta, file(row.read1), file(row.read2)] }
       .set { bbduking }

/* 0.01: Input trimming
*/
process bbduk {

  label 'med_mem'

  publishDir path: "$params.outDir/samples/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from bbduking

  output:
  file('*.txt') into log_bbduk
  tuple val(type), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming
  tuple val(type), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz'), file(read1), file(read2) into fastping
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) into fastqcing

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  {
  sh bbduk.sh -Xmx$taskmem \
    in1=$read1 \
    in2=$read2 \
    out1=$sampleID".bbduk.R1.fastq.gz" \
    out2=$sampleID".bbduk.R2.fastq.gz" \
    k=31 \
    mink=5 \
    hdist=1 \
    ktrim=r \
    trimq=20 \
    qtrim=rl \
    maq=20 \
    ref=/opt/miniconda/envs/somatic_n-of-1/opt/bbmap-38.57-0/resources/adapters.fa \
    tpe \
    tbo \
    stats=$sampleID".bbduk.adapterstats.txt" \
    overwrite=T
  } 2>&1 | tee > ${sampleID}.bbduk.runstats.txt
  """
}

/* 0.2: fastp QC of pre-, post-bbduk
*/
process fastp {

  label 'low_mem'

  publishDir "$params.outDir/samples/$sampleID/fastp", mode: "copy", pattern: "*.html"

  input:
  tuple val(type), val(sampleID), val(meta), file(preread1), file(preread2), file(postread1), file(postread2) from fastping

  output:
  file('*.html') into fastp_html
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h $sampleID"_pre.fastp.html" -j $sampleID"_pre.fastp.json" --in1 $preread1 --in2 $preread2

  fastp -w ${task.cpus} -h $sampleID"_post.fastp.html" -j $sampleID"_post.fastp.json" --in1 $postread1 --in2 $postread2
  """
}

/* 0.3: fastQC of per, post-bbduk
*/
process fastqc {

  label 'low_mem'

  publishDir "$params.outDir/samples/$sampleID/fastqc", mode: "copy", pattern: "*.html"

  input:
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from fastqcing

  output:
  file('*.html') into fastqc_multiqc

  script:
  """
  #!/bin/bash
  fastqc -t ${task.cpus} --noextract -o ./ $read1 $read2
  """
}

/* 1.0: Input alignment
*/
process bwamem {

  label 'high_mem'

  input:
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from bwa_memming
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])

  output:
  tuple val(type), val(sampleID), val(meta), file('*.bam'), file('*.bai') into (cramming, dup_marking)

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:$type\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    $fasta \
    $read1 $read2 | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"
  samtools index $sampleID".sort.bam"
  """
}

/* 1.1: CRAM alignment and output
* TODO: output upload schema for ENA/EGA
*/
process cram {

  label 'low_mem'

  publishDir path: "$params.outDir/samples/$sampleID/bwa", mode: "copy", pattern: "*.cra*"

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from cramming
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple file('*.cram'), file('*.crai') into completedcram

  script:
  """
  samtools view -hC -T $fasta $sampleID".sort.bam" > $sampleID".sort.cram"
  samtools index $sampleID".sort.cram"
  """
}

/* 1.2: MarkDuplicates
*/
process mrkdup {

  label 'high_mem'

  publishDir path: "$params.outDir/samples/$sampleID/picard", mode: "copy", pattern: "*.txt]"

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from dup_marking

  output:
  file('*.txt') into mrkdup_multiqc
  tuple val(type), val(sampleID), val(meta), file('*.md.bam'), file('*.md.bam.bai') into ( gatk4recaling, gridssing )

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  OUTBAM=\$(echo $bam | sed 's/bam/md.bam/')
  OUTMET=\$(echo $bam | sed 's/bam/md.metrics.txt/')
  {
  picard -Xmx$taskmem \
    MarkDuplicates \
    TMP_DIR=./ \
    INPUT=$bam \
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

/* 1.3: GATK4 BestPractices
*/
process gtkrcl {

  label 'high_mem'

  publishDir path: "$params.outDir/samples/$sampleID/gatk4/bestpractice", mode: "copy", pattern: "*.GATK4_BQSR.log.txt "

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from gatk4recaling
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(intlist) from Channel.value(params.intlist)

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( germfiltering, gmultimetricing)
  tuple val(type), val(sampleID), val(meta), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into hc_germ
  tuple val(sampleID), val(meta) into metas_pcgr
  file("${sampleID}.GATK4_BQSR.log.txt") into bqsr_log

  script:
  """
  {
  gatk BaseRecalibrator \
    -R $fasta \
    -I $bam \
    --known-sites $dbsnp \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    --disable-sequence-dictionary-validation true \
    -L $intlist

  #ApplyBQSR
  OUTBAM=\$(echo $bam | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R $fasta \
    -I $bam \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L $intlist

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee >  ${sampleID}.GATK4_BQSR.log.txt
  """
}

/* 1.31: scatter-gather implementation for mutect2, lancet
*/
process scat_gath {

  label 'low_mem'

  input:
  file(intlist) from Channel.value(params.intlist)

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
  picard IntervalListTools \
    I=$intlist \
    SCATTER_COUNT=$sgcount \
    O=./
  ls temp*/* | while read FILE; do
    COUNTN=\$(dirname \$FILE | perl -ane '@s=split(/\\_/); print \$s[1];');
    mv \$FILE mutect2.scatgath.\${COUNTN}.bed.interval_list;
    cp mutect2.scatgath.\${COUNTN}.bed.interval_list hc.scatgath.\${COUNTN}.bed.interval_list
    grep -v @ mutect2.scatgath.\${COUNTN}.bed.interval_list | \
      cut -f 1,2,3,5 > lancet.scatgath.\${COUNTN}.bed
  done
  """
}

/* 1.4: GATK4 Germline Haplotypecaller
*/
hcbedding = hc_bedding.flatten()
hc_germ
  .map { it -> [it[0],it[1],it[2],it[3],it[4]] }
  .combine(hcbedding)
  .set { hcgermbedding }

process haplotypecaller {

  label 'med_mem'

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai), file(intlist) from hcgermbedding
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  tuple file(omni), file(otbi), file(kgp1), file(ktbi), file(hpmp), file(htbi) from Channel.value([params.omni, params.otbi, params.kgp1, params.ktbi, params.hpmp, params.htbi])

  output:
  tuple val(sampleID), file('*sort.hc.vcf') into hc_gt
  tuple val(sampleID), val(meta) into hc_mv
  val(sampleID) into ( gridssgermID, vcfGRaID )

  when:
  type == "germline" & params.germline != false

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  SCATGATHN=\$(echo $intlist | perl -ane '@s=split(/\\./);print \$s[2];')
  gatk --java-options -Xmx$taskmem HaplotypeCaller \
    -R $fasta \
    -I $bam \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp $dbsnp \
    --native-pair-hmm-threads ${task.cpus} \
    -O $sampleID".\${SCATGATHN}.hc.vcf" \
    --disable-sequence-dictionary-validation true \
    -L $intlist

  picard SortVcf \
    I=$sampleID".\${SCATGATHN}.hc.vcf" \
    O=$sampleID".\${SCATGATHN}.sort.hc.vcf" \
    SD=$dict
  """
}


/* 1.5: GRIDSS SV calling in WGS
* because we do not know order or number of samples, create tuple with
* dummy as first and all others as list of elements
* then count them inside process
*/
gridssing
  .collect()
  .map { it -> tuple(it.flatten()) }
  // .println { it }
  .set { gridssin }



process gridss {

  label 'high_mem'

  publishDir path: "$params.outDir/output/gridss", mode: "copy"

  input:
  file(listbams) from gridssin
  val(germlineID) from gridssgermID.collect().flatten().unique()
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])
  tuple file(gblist), file(gprops), file(gbedpe), file(gbedse) from Channel.value([params.blacklist, params.gridssProps, params.gridssponbedpe, params.gridssponsinbed])

  output:
  file('*') into completegridss

  when:
  params.seqlevel == "wgs" & params.grchvers == "grch37"

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  GERMLINEBAM=\$(ls | grep $germlineID | grep bam | grep -v bai)
  BAMFILES=\$(echo -n \$GERMLINEBAM" "\$(ls *.bam | grep -v \$GERMLINEBAM))
  LABELS=\$(echo -n $germlineID" "\$(ls *bam | grep -v $germlineID | cut -d "." -f1) | sed 's/\\s */,/g')
  TUMORDS=\$(echo \$LABELS | perl -ane '@s=split(/\\,/);for(\$i=2;\$i<=@s;\$i++){push(@o,\$i);} print join(",",@o[0..\$#o]) . "\\n";')

  gridss \
    --reference $fasta \
    --output ${params.runID}".output.vcf.gz" \
    --assembly ${params.runID}".assembly.bam" \
    --threads ${task.cpus} \
    --jar /opt/miniconda/envs/gridss/share/gridss-2.9.3-0/gridss.jar\
    --workingdir ./ \
    --jvmheap $taskmem \
    --blacklist $gblist \
    --steps All \
    --configuration $gprops \
    --maxcoverage 50000 \
    --labels \$LABELS \
    \$BAMFILES

  Rscript --vanilla /opt/miniconda/envs/gridss/share/gridss/scripts/gridss_somatic_filter.R \
    --input ${params.runID}".output.vcf.gz" \
    --output ${params.runID}".somatic_filter.vcf.bgz" \
    --plotdir ./ \
    --scriptdir /opt/miniconda/envs/gridss/share/gridss-2.9.3-0 \
    --normalordinal 1 \
    --tumourordinal \$TUMORDS
  """
}

// 2.42: HC_merge
hc_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1..-1].flatten()) }
  .set { hc_fm }

process hc_merge {

  label 'high_mem'

  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*.vcf.*'

  input:
  tuple val(sampleID), file(rawvcfs) from hc_fm
  tuple val(sampleID), val(meta) from hc_mv
  output:
  tuple val(sampleID), val(meta1), file("${sampleID}.hc.merge.vcf.gz"), file("${sampleID}.hc.merge.vcf.gz.tbi") into cpsr_vcf

  script:
  meta1 = "${meta}"
  """
  ls *.sort.hc.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=$sampleID".hc.merge.vcf"
  bgzip $sampleID".hc.merge.vcf"
  tabix $sampleID".hc.merge.vcf.gz"
  """
}

/* 1.25: CPSR annotation of GATK4 Germline
*/
process cpsrreport {

  label 'med_mem'

  publishDir "$params.outDir/reports", mode: "copy", pattern: "*.html"
  publishDir "$params.outDir/output/cpsr", mode: "copy", pattern: "*[!.html]"

  input:
  tuple val(sampleID), val(meta), file(vcf), file(tbi) from cpsr_vcf

  output:
  file('*') into cpsr_vcfs

  script:
  """
  {
  CONFIG=\$(readlink -e ${params.pcgrdir}/data/*/cpsr_configuration_default.toml)
  META=\$(echo $meta | sed 's/\\s */_/g' | sed 's/[()]//g')
  VERS=\$(ls ${params.pcgrdir}/data)
  ##CPSR v0.5.2.2
  cpsr.py \
    --no-docker \
    --no_vcf_validate \
    --panel_id 0 \
    $vcf \
    ${params.pcgrdir} \
    ./ \
    \$VERS \
    \$CONFIG \
    \$META
  } 2>&1 | tee > ${sampleID}.cpsr.log.txt
  """
}

/* 2.0: filter germline channel, tap into somatic channels for all processes subsequent
*/
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

/* 2.1: PicardTools metrics suite for MultiQC HTML report
*/
process mltmet {

  label 'med_mem'

  publishDir "$params.outDir/samples/$sampleID/metrics"

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from gmultimetricing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(intlist) from Channel.value(params.intlist)

  output:
  file('*.txt') into multimetrics_multiqc

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  {
  if [[ ${params.seqlevel} == "exome" ]]; then
  picard -Xmx$taskmem CollectHsMetrics \
    I=$bam \
    O=$sampleID".hs_metrics.txt" \
    TMP_DIR=./ \
    R=$fasta \
    BAIT_INTERVALS=$intlist  \
    TARGET_INTERVALS=$intlist
  fi
  picard -Xmx$taskmem CollectAlignmentSummaryMetrics \
    I=$bam \
    O=$sampleID".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    R=$fasta

  picard -Xmx$taskmem CollectMultipleMetrics \
    I=$bam \
    O=$sampleID".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    R=$fasta

  picard -Xmx$taskmem CollectSequencingArtifactMetrics \
    I=$bam \
    O=$sampleID".artifact_metrics.txt" \
    TMP_DIR=./ \
    R=$fasta

  picard -Xmx$taskmem CollectInsertSizeMetrics \
    I=$bam \
    O=$sampleID".insert_size_metrics.txt" \
    H=$bam".histogram.pdf" \
    TMP_DIR=./

  } 2>&1 | tee > ${sampleID}.picard.metrics.log
  """
}

/*2.21: SCNA with facets CSV snp-pileup
*/

process fctcsv {

  label 'med_mem'

  publishDir "$params.outDir/samples/$sampleID/facets"

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from facetsomaing
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])

  output:
  file("${sampleID}.cncf-jointsegs.pcgr.tsv") into facets_consensusing
  tuple val(sampleID), file("${sampleID}.cncf-jointsegs.pcgr.tsv"), file("${sampleID}.fit_ploidy-purity.pcgr.tsv") into facets_pcgr
  file("${sampleID}.facets.log.txt") into log_facets

  script:
  """
  {
  snp-pileup \
    $dbsnp \
    -r 10 \
    -p \
    ${sampleID}.facets.r10.csv \
    $germlinebam \
    $tumourbam

  Rscript --vanilla  ${workflow.projectDir}/bin/facets_cna.call.R ${sampleID}.facets.r10.csv

  echo -e "Chromosome\\tStart\\tEnd\\tSegment_Mean" > $sampleID".cncf-jointsegs.pcgr.tsv"
  tail -n+2 $sampleID".fit_cncf-jointsegs.tsv" | awk '{print \$1"\\t"\$10"\\t"\$11"\\t"\$5}' >> $sampleID".cncf-jointsegs.pcgr.tsv"

  tail -n+2 $sampleID".fit_ploidy-purity.tab" > $sampleID".fit_ploidy-purity.pcgr.tsv"
  } 2>&1 | tee > ${sampleID}.facets.log.txt
  """
}

/* 2.22: SCNA consensus from facets
*/
process fctcon {

  label 'med_mem'

  publishDir "$params.outDir/output/scna/facets"

  input:
  file(filesn) from facets_consensusing.collect()
  file(dict) from Channel.value(params.dict)

  output:
  file('*') into complete_facets

  script:
  """
  {
  Rscript --vanilla  ${workflow.projectDir}/bin/facets_cna_consensus.call.R \
    $dict \
    ./ \
    ${workflow.projectDir}/bin/facets_cna_consensus.func.R
  } 2>&1 | tee > facets_cons.log.txt
  """
}

mutect2bedding = mutect2_bedding.flatten()
mutect2somaticing
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5]]}
  .combine(mutect2bedding)
  .set { mutect2somaticbedding }

/* 2.41: MuTect2
* NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error
*/
process mutct2_sg {

  label 'med_mem'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(intlist) from mutect2somaticbedding
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple val(sampleID), file('*sort.mutect2.vcf') into mutect2_gt
  tuple val(sampleID), file('*.vcf.stats') into mutect2_st
  tuple val(sampleID), file('*mutect2.f1r2.tar.gz') into mutect2_f1r2

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  SCATGATHN=\$(echo $intlist | perl -ane '@s=split(/\\./);print\$s[2];')
  gatk --java-options -Xmx$taskmem \
    Mutect2 \
    --native-pair-hmm-threads ${task.cpus} \
    --reference $fasta \
    --input $germlinebam \
    --input $tumourbam \
    --normal-sample $germlineID \
    --tumor-sample $sampleID \
    --output $sampleID"."\${SCATGATHN}".mutect2.vcf" \
    --disable-sequence-dictionary-validation true \
    --f1r2-tar-gz \${SCATGATHN}".mutect2.f1r2.tar.gz" \
    -L $intlist

  picard SortVcf \
    I=$sampleID"."\${SCATGATHN}".mutect2.vcf" \
    O=$sampleID"."\${SCATGATHN}".sort.mutect2.vcf" \
    SD=$dict
  """
}

// 2.42: MuTect2_merge
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
  picard MergeVcfs I=vcf.list O=$sampleID".mutect2.merge.vcf"
  """
}

mutect2_st
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
  .set { mutect2_sm }

/* 2.43: MuTect2 Concatenate VCFs
*/
process mutct2_concstat {

  label 'med_mem'

  input:
  tuple val(sampleID), file(stats) from mutect2_sm

  output:
  tuple val(sampleID), file('*mutect2.merge.vcf.stats') into mutect2_stats

  script:
  """
  STATS=\$(ls *stats | perl -ane 'foreach \$k (@F){print "--stats \$k ";}')
  gatk MergeMutectStats --output $sampleID".mutect2.merge.vcf.stats" \$STATS
  """
}

/* 2.44: MuTect2 Concatenate VCFs
*/
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

/* 2.44: MuTect2 Contamination
*/
mutect2_contam
  .join(mutect2_merge)
  .join(mutect2_stats)
  .join(mutect2_f1r2_comb)
  .groupTuple()
  .map { it -> [it[0], it[1..5].flatten(), it[6], it[7], it[8]].flatten() }
  .set { mutect2_contam_merge }

process mutct2_contam_filter {

  label 'med_mem'

  publishDir path: "$params.outDir/samples/$sampleID/mutect2", mode: "copy", overwrite: true
  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*raw.vcf', overwrite: true
  publishDir path: "$params.outDir/samples/", mode: "copy", pattern: '*issue.table', overwrite: true

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(mergevcf), file(statsvcf), file(readorient) from mutect2_contam_merge
  tuple file(fasta), file(fai), file(dict), file(gps), file(gpstbi), file(intlist) from Channel.value([params.fasta, params.fai, params.dict, params.gps, params.gpstbi, params.intlist])

  output:
  file('*snv_indel.pass.vcf') into mutect2_veping mode flatten
  file('*.raw.vcf') into mutect2_rawVcf
  file('*') into completedmutect2call

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  gatk --java-options -Xmx$taskmem \
    GetPileupSummaries \
    -I $tumourbam \
    -V $gps \
    -O $sampleID".getpileupsummaries.table" \
    -L $intlist

  gatk CalculateContamination \
    -I $sampleID".getpileupsummaries.table" \
    -O $sampleID".calculatecontamination.table"

  Rscript --vanilla ${workflow.projectDir}/bin/MuTect2_contamination.call.R $sampleID".calculatecontamination.table" $sampleID

  gatk IndexFeatureFile \
    --feature-file $mergevcf

  gatk --java-options -Xmx$taskmem \
    FilterMutectCalls \
    --reference $fasta \
    --contamination-table $sampleID".calculatecontamination.table" \
    --interval-padding 5 \
    --output $sampleID".mutect2.FilterMutectCalls.vcf" \
    --unique-alt-read-count 3 \
    --variant $mergevcf \
    --stats $statsvcf \
    --disable-sequence-dictionary-validation true \
    --ob-priors $readorient \
    -L $intlist

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=$sampleID \
    DP=14 \
    MD=2 \
    VCF=$sampleID".mutect2.FilterMutectCalls.vcf"
  """
}

/* 2.5: Manta output is a pre-req for Strelka2, so call both here
*/
process mntstr {

  label 'high_mem'

  publishDir path: "$params.outDir/samples/$sampleID/manta-strelka2", mode: "copy"
  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*[.strelka2.snv_indel.raw.vcf, .strelka2.snv_indel.pass.vcf]'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(bedgz), file(bedgztbi) from Channel.value([params.bedgz, params.bedgztbi])

  output:
  file("${sampleID}.strelka2.snv_indel.pass.vcf") into strelka2_veping
  file("${sampleID}.strelka2.raw.vcf") into strelka2_rawVcf
  file('*.txt') into log_mantastrelka

  script:
  """
  {
  if [[ ${params.seqlevel} == "exome" ]];then
    CR="--exome --callRegions $bedgz"
  else
    CR="--callRegions $bedgz"
  fi

  configManta.py \$CR --referenceFasta=$fasta --normalBam=$germlinebam --tumourBam=$tumourbam --runDir=manta

  manta/runWorkflow.py -m local -j ${task.cpus}

  configureStrelkaSomaticWorkflow.py \$CR --referenceFasta=$fasta --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz --normalBam=$germlinebam --tumorBam=$tumourbam --runDir=strelka2

  strelka2/runWorkflow.py -m local -j ${task.cpus}

  ##merge into raw snv_indel
  gatk MergeVcfs -I strelka2/results/variants/somatic.snvs.vcf.gz -I strelka2/results/variants/somatic.indels.vcf.gz -O tmp.strelka2.snv_indel.vcf

  ${workflow.projectDir}/bin/manta_strelka2_rename_filter.sh  tmp.strelka2.snv_indel.vcf tmp2.strelka2.snv_indel.vcf ${sampleID} ${germlineID}

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=tmp2.strelka2.snv_indel.vcf

  } 2>&1 | tee > ${sampleID}.manta-strelka2.log.txt
  """
}

/* 2.61: Lancet
*/
lancetbedding = lancet_bedding.flatten()
lanceting
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5]]}
  .combine(lancetbedding)
  .set { lancetsbedding }

process lancet_sg {

  label 'med_mem'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(bed) from lancetsbedding
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple val(sampleID), file('*.sort.lancet.vcf') into lancet_gt

  when:
  params.seqlevel == "exome"

  script:
  """
  SCATGATHN=\$(echo $bed | perl -ane '@s=split(/\\./);print\$s[2];')
  lancet \
    --num-threads ${task.cpus} \
    --ref $fasta \
    --bed $bed \
    --tumor $tumourbam \
    --normal $germlinebam | \
    perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
      \$_=~s/TUMOR/$sampleID/;
      \$_=~s/NORMAL/$germlineID/;
      print \$_;}
    else{print \$_;}' > $sampleID"."\${SCATGATHN}".lancet.vcf"

  picard SortVcf \
    I=$sampleID"."\${SCATGATHN}".lancet.vcf" \
    O=$sampleID"."\${SCATGATHN}".sort.lancet.vcf" \
    SD=$dict
  """
}

/* 2.62: Lancet Merge
*/
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
  picard MergeVcfs I=vcf.list O=$sampleID".lancet.merge.vcf"
  """
}

/* 2.63: Lancet Filter
*/
process lancet_filter {

  label 'med_mem'

  publishDir path: "$params.outDir/samples/$sampleID/lancet"
  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  tuple val(sampleID), file(mergevcf) from lancet_merge

  output:
  file('*.pass.vcf') into lancet_veping mode flatten
  file('*.raw.vcf') into lancet_rawVcf
  file('*') into completedlancetcall

  script:
  """
  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=$sampleID \
    DP=14 \
    MD=2 \
    VCF=$mergevcf
  """
}

/* 3.0: Annotate Vcfs
*/
ALLVCFS = lancet_veping
          .mix( mutect2_veping )
          .mix( strelka2_veping )

process vepann {

  label 'med_mem'

  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*.vcf'

  input:
  each file(vcf) from ALLVCFS
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  file('*.vcf') into runGRanges

  script:
  """
  VCFANNO=\$(echo $vcf | sed "s/.vcf/.vep.vcf/")
  GRCHVER=\$(ls ${params.grchvers})
  VEPVERS=\$(ls ${params.pcgrdir}/data/*/.vep/homo_sapiens/ | cut -d "_" -f2)
  vep --dir_cache ${params.pcgrdir}/data/\$GRCHVER/.vep \
    --offline \
    --assembly \$VEPVERS \
    --vcf_info_field ANN \
    --symbol \
    --species homo_sapiens \
    --check_existing \
    --cache \
    --fork ${task.cpus} \
    --af_1kg \
    --af_gnomad \
    --vcf \
    --input_file $vcf \
    --output_file \$VCFANNO \
    --format "vcf" \
    --fasta $fasta \
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

process vcfGRa {

  label 'med_mem'

  publishDir "$params.outDir/output/pdf", pattern: '*.pdf'
  publishDir "$params.outDir/output/vcf", pattern: '*.vcf'
  publishDir "$params.outDir/output/data", pattern: '*.[*RData,*tab]'

  input:
  file(grangesvcfs) from allvcfs
  val(germlineID) from vcfGRaID.unique()
  tuple file(conscallR), file(consfuncR) from vcfGRa_Scripts

  output:
  file('*.ALL.pcgr.all.tab.vcf') into vcfs_pcgr
  file('*') into completedvcfGRangesConsensus

  script:
  """
  Rscript --vanilla $conscallR \
    $consfuncR \
    $germlineID \
    "snv_indel.pass.vep.vcf" \
    \"${params.incOrder}\"
  """
}

vcfs_pcgr
  .flatten()
  .set { vcfs_pcgrd }

/* 3.2 Create VCF for PCGR from consensus
*/
process pcgrVcf {

  label 'low_mem'

  input:
  file(vcf) from vcfs_pcgrd

  output:
  tuple val(sampleID), file("${sampleID}.snv_indel.pass.pcgr.vcf") into snvpass_pcgr

  script:
  sampleID = "${vcf}".split("\\.")[0]
  """
  for VCF in *ALL.pcgr.all.tab.vcf; do
    NVCF=\$(echo \$VCF | sed 's/ALL.pcgr.all.tab.vcf/snv_indel.pass.pcgr.vcf/')
    cat ${workflow.projectDir}/bin/vcf42.head.txt > \$NVCF
    head -n1 \$VCF >> \$NVCF
    tail -n+2 \$VCF | sort -V >> \$NVCF
  done
  """
}

snvpass_pcgr
  .join(metas_pcgr)
  .join(facets_pcgr)
  .map { it -> tuple(it[0], it[1..-1]).flatten() }
  .set { pcgr_inputs }

/* 3.3 PCGR report
* take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR
*/
process pcgrreport {

  label 'low_mem'

  publishDir "$params.outDir/reports", mode: "copy", pattern: "*html"
  publishDir "$params.outDir/samples/$sampleID/pcgr", mode: "copy", pattern: "*[!.html]"

  input:
  tuple val(sampleID), file(vcf), val(meta), file(jointsegs), file(ploidpur) from pcgr_inputs

  output:
  file('*') into completedPCGR

  script:
  """
  {
  ##want META to allow spaces, remove non-alphanum
  META=\$(echo $meta | perl -ane '\$_=~s/[^a-zA-Z0-9_ \\n]//g; print \$_;' | sed 's/\\s */_/g')

  PLOIDY=""; PURITY="";
  if [[ \$(cut -f 1 $ploidpur) != "NA" ]]; then
    PLOIDY="--tumor_ploidy \$(cut -f 1 $ploidpur)"
  fi
  if [[ \$(cut -f 2 $ploidpur) != "NA" ]]; then
    PURITY="--tumor_purity \$(cut -f 2 $ploidpur)"
  fi

  CONFIG=\$(readlink -e ${params.pcgrdir}/data/*/pcgr_configuration_default.toml)
  VERS=\$(ls ${params.pcgrdir}/data)
  pcgr.py ${params.pcgrdir} \
    ./ \
    \$VERS \
    \$CONFIG \
    \$META \
    --input_vcf $vcf \
    --input_cna $jointsegs \$PLOIDY \$PURITY \
    --no-docker \
    --force_overwrite

  } 2>&1 | tee > ${sampleID}.pcgr.log.txt
  """
}

/* 4.0 Run multiQC to finalise report
*/
process mltiQC {

  label 'low_mem'

  publishDir path: "$params.outDir/reports", mode: "copy", pattern: "*html"

  input:
  file(fastps) from fastp_multiqc.collect()
  file(fastqcs) from fastqc_multiqc.collect()
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()
  file(mrkdups) from mrkdup_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  multiqc . -i ${params.runID}".somatic_n-of-1" --tag DNA -f -c ${params.multiqcConfig}
  """
}
