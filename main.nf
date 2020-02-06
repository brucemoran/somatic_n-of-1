#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '------------------------------------------------------------'
  log.info 'NEXTFLOW 19.10 BAM QC, TRIM, ALIGN, SOMATIC SNV, CNA, REPORT'
  log.info '------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run main.nf \
              --sampleCsv "data/sample.csv" \
              --refDir "refs" \
              --includeOrder "tumour_A,tumour_B" \
              --germline \
              --multiqcConfig "~/.multiqc/my.config"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --sampleCsv      STRING      CSV format, headers: type (either "germline" or "somatic"),sampleID,meta,/path/to/read1.fastq.gz,/path/to/read2.fastq.gz ;use meta for naming in PCGR, CPSR report'
  log.info '    --refDir      STRING      dir in which reference data and required indices held; recommended to run associated reference creation NextFlow, DNAseq_references'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --includeOrder      STRING      in final plots, use this ordering of samples (if multiple somatic samples); comma-separated, no spaces'
  log.info '    --germline      STRING      run HaplotypeCaller on germline sample and annotate with CPSR'
  log.info '    --multiqcConfig      STRING      config file for multiqc'
  log.info ''
  exit 1
}

/* -2: Global Variables
*/
params.outDir = "analysis"

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

params.exomeintlist = Channel.fromPath("$params.refDir/exome.bed.interval_list").getVal()
params.exomebed = Channel.fromPath("$params.refDir/exome.bed").getVal()
params.exomebedgz = Channel.fromPath("$params.refDir/exome.bed.gz").getVal()
params.exomebedgztbi = Channel.fromPath("$params.refDir/exome.bed.gz.tbi").getVal()

params.dbsnp = Channel.fromPath("$params.refDir/dbsnp*.gz").getVal()
params.dbsnptbi = Channel.fromPath("$params.refDir/dbsnp*.tbi").getVal()
params.omni = Channel.fromPath("$params.refDir/KG_omni*.gz").getVal()
params.otbi = Channel.fromPath("$params.refDir/KG_omni*.gz.tbi").getVal()
params.kgp1 = Channel.fromPath("$params.refDir/KG_phase1*.gz").getVal()
params.ktbi = Channel.fromPath("$params.refDir/KG_phase1*.gz.tbi").getVal()
params.hpmp = Channel.fromPath("$params.refDir/hapmap*.gz").getVal()
params.htbi = Channel.fromPath("$params.refDir/hapmap*.gz.tbi").getVal()

params.gps = Channel.fromPath("$params.refDir/af-only-gnomad.*.vcf.gz").getVal()
params.gpstbi = Channel.fromPath("$params.refDir/af-only-gnomad.*.vcf.gz.tbi").getVal()

//PCGR, CPSR version and base data dir
params.grchvers = Channel.fromPath("$params.refDir/pcgr/data").getVal()
Channel.fromPath("$params.refDir/pcgr/", type: 'dir').into { CPSR; PCGR }

/* -1 Echo PATHs
*/
vers = Channel.fromPath("$params.grchvers", type: 'dir')
vepp = Channel.fromPath("$params.grchvers", type: 'dir')

process grchvers {
  executor 'local'

  input:
  file(datDir) from vers

  output:
  stdout into (cpsr_grchvers, pcgr_grchvers)

  """
  GRCHVERS=\$(ls $datDir/)
  echo -n \$GRCHVERS
  """
}

process vepp {
  executor 'local'
  echo true

  input:
  file(datDir) from vepp

  output:
  file('.vep/') into vepanndir

  """
  GRCHVERS=\$(ls $datDir)
  ln -s \$(readlink -e ${params.refDir}/pcgr/$datDir/\$GRCHVERS/.vep) ./
  """
}

/* 0.00: Input using sample.csv
*/
sampleInputs = Channel.fromPath("$params.sampleCsv")

/* 0.01: Push metadata on to CPSR, PCGR
*/
process meta {
  executor 'local'
  input:
  file(txt) from sampleInputs

  output:
  file('inputs.txt') into (sampleCsvInput, pcgrMetaInput)

  script:
  """
  cat $txt > "inputs.txt"
  """
}
//split for inputs
sampleCsvInput.splitCsv( header: true )
              .map { row -> [row.type, row.sampleID, row.meta, file(row.read1), file(row.read2)] }
              .set { bbduking }

/* 0.0: Input trimming
*/
process bbduk {

  publishDir path: "$params.outDir/samples/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from bbduking

  output:
  file('*') into completed0_0
  tuple val(type), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming
  tuple val(type), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz'), file(read1), file(read2) into fastping
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) into fastqcing

  script:
  """
  {
    sh bbduk.sh ${params.quarter_javamem} \
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
      ref=/opt/miniconda/envs/somatic_exome_n-of-1/opt/bbmap-38.57-0/resources/adapters.fa \
      tpe \
      tbo \
      stats=$sampleID".bbduk.adapterstats.txt" \
      overwrite=T
  } 2>&1 | tee > $sampleID".bbduk.runstats.txt"
  """
}

/* 0.1: QC of per, post-bbduk
*/
process fastp {

  publishDir "$params.outDir/samples/$sampleID/fastp", mode: "copy", pattern: "*.html"

  input:
  tuple val(type), val(sampleID), val(meta), file(preread1), file(preread2), file(postread1), file(postread2) from fastping

  output:
  file('*.html') into completed0_1
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h $sampleID"_pre.fastp.html" -j $sampleID"_pre.fastp.json" --in1 $preread1 --in2 $preread2

  fastp -w ${task.cpus} -h $sampleID"_post.fastp.html" -j $sampleID"_post.fastp.json" --in1 $postread1 --in2 $postread2
  """
}

/* 0.1: QC of per, post-bbduk
*/
process fastqc {

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

  cache 'deep'

  publishDir path: "$params.outDir/samples/$sampleID/bwa", mode: "copy", pattern: "*.log.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(read1), file(read2) from bwa_memming
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])

  output:
  tuple val(type), val(sampleID), val(meta), file('*.bam'), file('*.bai') into (cramming, dup_marking)
  file('*.log.txt') into completedbwa

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:$type\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  {
  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    $fasta \
    $read1 $read2 | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"

  samtools index $sampleID".sort.bam"

  } 2>&1 | tee > $sampleID".bwa-mem.log.txt"
  """
}

/* 1.01: CRAM alignment and output
* WIP: output upload schema for ENA/EGA
*/
process cram {

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

/* 1.1: MarkDuplicates
*/
process mrkdup {

  publishDir path: "$params.outDir/samples/$sampleID/picard/markdup", mode: "copy", pattern: "*[!.metrics.txt]"
  publishDir path: "$params.outDir/samples/$sampleID/picard/metrics", mode: "copy", pattern: "*.metrics.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from dup_marking

  output:
  file('*md.metrics.txt') into mrkdup_multiqc
  tuple val(type), val(sampleID), val(meta), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

  script:
  """
  OUTBAM=\$(echo $bam | sed 's/bam/md.bam/')
  OUTMET=\$(echo $bam | sed 's/bam/md.metrics.txt/')
  {
    picard ${params.quarter_javamem} \
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
  } 2>&1 | tee > $sampleID".picard_markDuplicates.log.txt"
  """
}

/* 1.2: GATK4 BestPractices
* as per best practices of GATK4
*/
process gtkrcl {

  publishDir path: "$params.outDir/samples/$sampleID/gatk4/bestpractice", mode: "copy", pattern: "*.log.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from gatk4recaling
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( germfiltering, gmultimetricing)
  tuple val(type), val(sampleID), val(meta), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into gatk_germ

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
    -L $exomeintlist

  #ApplyBQSR
  OUTBAM=\$(echo $bam | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R $fasta \
    -I $bam \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L $exomeintlist

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee > $sampleID".GATK4_recal.log.txt"
  """
}

/* 1.21: GATK4 Germline
*/
process gatkgerm {

  publishDir "$params.outDir/samples/$sampleID/gatk4/HC_germline", mode: "copy", pattern: "*.log.txt"

  input:
  tuple val(type), val(sampleID), val(meta), file(bam), file(bai) from gatk_germ
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  tuple file(omni), file(otbi), file(kgp1), file(ktbi), file(hpmp), file(htbi) from Channel.value([params.omni, params.otbi, params.kgp1, params.ktbi, params.hpmp, params.htbi])

  output:
  tuple val(type), val(sampleID), val(meta), file('*.HC.vcf.gz'), file('*.HC.vcf.gz.tbi') into germ_vcf

  script:
  """
  {
    #HaplotypeCaller
    INPUTBAM=$bam
    OUTVCF=\$(echo \$INPUTBAM | sed 's/bam/hc.vcf/')
    gatk --java-options ${params.full_javamem} HaplotypeCaller \
      -R $fasta \
      -I \$INPUTBAM \
      --dont-use-soft-clipped-bases \
      --standard-min-confidence-threshold-for-calling 20 \
      --dbsnp $dbsnp \
      --native-pair-hmm-threads ${task.cpus} \
      -O $sampleID".HC.vcf" \
      --disable-sequence-dictionary-validation true \
      -L $exomeintlist

    bgzip $sampleID".HC.vcf"
    tabix $sampleID".HC.vcf.gz"

  } 2>&1 | tee $sampleID".GATK4_HaplotypeCaller-germline.log.txt"
  """
}

/* 1.22: CPSR annotation of GATK4 Germline
*/
process cpsrreport {

  publishDir "$params.outDir/reports", mode: "copy", pattern: "*.html"
  publishDir "$params.outDir/output/cpsr", mode: "copy", pattern: "*[!.html]"

  input:
  tuple val(type), val(sampleID), val(meta), file(vcf), file(tbi) from germ_vcf
  each file(cpsr_grch) from CPSR
  val(grchvers) from cpsr_grchvers

  output:
  file('*') into cpsr_vcfs

  when:
  type == "germline"

  script:
  """
  {
    CONFIG=\$(readlink -e $cpsr_grch/data/*/cpsr_configuration_default.toml)
    META=\$(echo $meta | sed 's/\\s */_/g')
    cpsr.py \
      --no-docker \
      --no_vcf_validate \
      $vcf \
      $cpsr_grch \
      ./ \
      $grchvers \
      0 \
      \$CONFIG \
      \$META

  } 2>&1 | tee $sampleID".cpsr.log.txt"
  """
}

/* filter germ channel for all processes subsequent
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

/* combine germline with somatic and unique those outputs
*/
somagerm.somafiltered
    .map { [it[1], it[2..-1]] }
    .tap { somatap }
somagerm.germfiltered
    .map { [it[1], it[2..-1]] }
    .tap { germtap }
somatap.combine(germtap).tap{ somagermtap }

process combinegs {

  input:
  tuple val(sampleID), file(somafiles), val(germlineID), file(germlinefiles) from somagermtap

  output:
  tuple val(sampleID), file("${somafiles[0]}"), file("${somafiles[1]}"), val(germlineID), file("${germlinefiles[0]}"), file("${germlinefiles[1]}") into ( mutect2somaticing, facetsomaing, msisensoring, mantastrelka2ing, lanceting )
  val(germlineID) into vcfGRaID

  script:
  """
  """
}

/* 2.0: PicardTools metrics suite for MultiQC HTML report
*/
process mltmet {

  publishDir "$params.outDir/samples/$sampleID/metrics"

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from gmultimetricing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  file('*') into all2_0
  file('*.txt') into multimetrics_multiqc

  script:
  """
  {
    picard CollectHsMetrics \
      I=$bam \
      O=$sampleID".hs_metrics.txt" \
      TMP_DIR=./ \
      R=$fasta \
      BAIT_INTERVALS=$exomeintlist  \
      TARGET_INTERVALS=$exomeintlist

    picard CollectAlignmentSummaryMetrics \
      I=$bam \
      O=$sampleID".AlignmentSummaryMetrics.txt" \
      TMP_DIR=./ \
      R=$fasta

    picard CollectMultipleMetrics \
      I=$bam \
      O=$sampleID".CollectMultipleMetrics.txt" \
      TMP_DIR=./ \
      R=$fasta

    picard CollectSequencingArtifactMetrics \
      I=$bam \
      O=$sampleID".artifact_metrics.txt" \
      TMP_DIR=./ \
      R=$fasta

    picard EstimateLibraryComplexity \
      I=$bam \
      O=$sampleID".est_lib_complex_metrics.txt" \
      TMP_DIR=./

    picard CollectInsertSizeMetrics \
      I=$bam \
      O=$sampleID".insert_size_metrics.txt" \
      H=$bam".histogram.pdf" \
      TMP_DIR=./

  } 2>&1 | tee > $sampleID".picard.metrics.log"
  """
}

/*2.1: SCNA with facets CSV snp-pileup
*/
process fctcsv {

  publishDir "$params.outDir/samples/$sampleID/facets"

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from facetsomaing
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])

  output:
  file('*.cncf-jointsegs.pcgr.tsv') into facets_consensusing
  file('*.pcgr.tsv') into facets_pcgr
  file('*') into facetsoutputR

  script:
  """
  CSVFILE=\$(echo $tumourbam | sed 's/bam/facets.r10.csv/')

  {
    snp-pileup \
      $dbsnp \
      -r 10 \
      -p \
      \$CSVFILE \
      $germlinebam \
      $tumourbam

    Rscript --vanilla  ${workflow.projectDir}/bin/facets_cna.call.R \$CSVFILE

    echo -e "Chromosome\\tStart\\tEnd\\tSegment_Mean" > $sampleID".cncf-jointsegs.pcgr.tsv"
    tail -n+2 $sampleID".fit_cncf-jointsegs.tsv" | awk '{print \$1"\\t"\$10"\\t"\$11"\\t"\$5}' >> $sampleID".cncf-jointsegs.pcgr.tsv"

    tail -n+2 $sampleID".fit_ploidy-purity.tab" > $sampleID".fit_ploidy-purity.pcgr.tsv"

  } 2>&1 | tee > $sampleID".facets_snpp_call.log.txt"
  """
}

/* 2.11: SCNA consensus from facets
*/
process fctcon {

  publishDir "$params.outDir/output/scna/facets"

  input:
  file(filesn) from facets_consensusing.collect()
  file(dict) from Channel.value(params.dict)

  output:
  file('*') into completed2_11

  script:
  """
  {
  OUTID=\$(basename ${workflow.launchDir})
  Rscript --vanilla  ${workflow.projectDir}/bin/facets_cna_consensus.call.R \
    $dict \
    \$OUTID \
    ${workflow.projectDir}/bin/facets_cna_consensus.func.R
  } 2>&1 | tee > "facets_cons.log.txt"
  """
}

/* 2.2: MuTect2
* NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error
*/
process mutct2 {

  publishDir path: "$params.outDir/samples/$sampleID/mutect2", mode: "copy"
  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mutect2somaticing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  tuple file(gps), file(gpstbi) from Channel.value([params.gps, params.gpstbi])

  output:
  file('*.pass.vcf') into mutect2_veping mode flatten
  file('*.raw.vcf') into mutect2_rawVcf
  file('*') into completedmutect2call
  tuple val(sampleID), file('*calculatecontamination.table') into contamination

  script:
  """
  {
    gatk --java-options ${params.full_javamem} \
      Mutect2 \
      --native-pair-hmm-threads ${task.cpus} \
      --reference $fasta \
      --input $germlinebam \
      --input $tumourbam \
      --normal-sample $germlineID \
      --tumor-sample $sampleID \
      --output $sampleID".md.recal.mutect2.vcf" \
      --disable-sequence-dictionary-validation true \
      -L $exomeintlist

    gatk --java-options ${params.full_javamem} \
      GetPileupSummaries \
      -I $tumourbam \
      -V $gps \
      -O $sampleID".getpileupsummaries.table" \
      -L $exomeintlist

    gatk CalculateContamination \
      -I $sampleID".getpileupsummaries.table" \
      -O $sampleID".calculatecontamination.table"

    gatk --java-options ${params.full_javamem} \
      FilterMutectCalls \
      --reference $fasta \
      --contamination-table $sampleID".calculatecontamination.table" \
      --interval-padding 5 \
      --output $sampleID".md.recal.mutect2.FilterMutectCalls.vcf" \
      --unique-alt-read-count 3 \
      --variant $sampleID".md.recal.mutect2.vcf" \
      --disable-sequence-dictionary-validation true \
      -L $exomeintlist

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=$sampleID".md.recal.mutect2.FilterMutectCalls.vcf"

  } 2>&1 | tee > $sampleID".GATK4_mutect2.log.txt"
  """

}

/* 2.31: MuTect2 Contamination
*/
process mutct2_contam {

  publishDir path: "$params.outDir/samples/", mode: "copy", pattern: '*issue.table'

  input:
  tuple val(sampleID), file(contable) from contamination

  output:
  file('*.table') into completedcontam

  """
  Rscript --vanilla ${workflow.projectDir}/bin/MuTect2_contamination.call.R $contable $sampleID
  """
}

/* 2.4: Manta output is a pre-req for Strelka2, so call both here
*/
process mntstr {

  publishDir path: "$params.outDir/samples/$sampleID/manta-strelka2"
  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(exomebedgz), file(exomebedgztbi) from Channel.value([params.exomebedgz, params.exomebedgztbi])

  output:
  val(sampleID) into completed2_4
  file('*.strelka2.snv_indel.pass.vcf') into strelka2_veping
  file('*.raw.vcf') into strelka2_rawVcf
  file('manta/*') into completedmantacall

  script:
  """
  {
    configManta.py \
      --exome \
      --referenceFasta=$fasta \
      --callRegions $exomebedgz \
      --normalBam=$germlinebam \
      --tumourBam=$tumourbam \
      --runDir=manta

    manta/runWorkflow.py -m local -j ${task.cpus}

    configureStrelkaSomaticWorkflow.py \
      --exome \
      --referenceFasta=$fasta \
      --callRegions $exomebedgz \
      --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz \
      --normalBam=$germlinebam \
      --tumorBam=$tumourbam \
      --runDir=strelka2

    strelka2/runWorkflow.py -m local -j ${task.cpus}

    TUMOURSNVVCF=\$(echo $tumourbam | sed 's/bam/strelka2.snv.vcf/')
    gunzip -c strelka2/results/variants/somatic.snvs.vcf.gz | \
    perl -ane 'if(\$F[0]=~m/^#/){if(\$_=~m/^#CHROM/){
        \$_=~s/NORMAL/$germlineID/;
        \$_=~s/TUMOR/$sampleID/;
        print \$_;next;}
        else{print \$_;next;}
      }
      else{print \$_;}' > \$TUMOURSNVVCF

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
     ID=$sampleID \
     DP=14 \
     MD=2 \
     VCF=\$TUMOURSNVVCF

    TUMOURINDELVCF=\$(echo $tumourbam | sed 's/bam/strelka2.indel.vcf/')
    gunzip -c strelka2/results/variants/somatic.indels.vcf.gz | \
    perl -ane 'if(\$F[0]=~m/^#/){if(\$_=~m/^#CHROM/){
        \$_=~s/NORMAL/$germlineID/;
        \$_=~s/TUMOR/$sampleID/;
        print \$_;next;}
        else{print \$_;next;}}
      else{print \$_;}' > \$TUMOURINDELVCF

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=\$TUMOURINDELVCF

    PASSSNV=\$(ls | grep strelka2.snv.pass.vcf)
    PASSIND=\$(ls | grep strelka2.indel.pass.vcf)

    gatk MergeVcfs \
      -I \$PASSSNV \
      -I \$PASSIND \
      -O $sampleID".strelka2.snv_indel.pass.vcf"

  } 2>&1 | tee > $sampleID".manta-strelka2.log.txt"
  """
}

/* 2.5: Lancet
*/
process lancet {

  publishDir path: "$params.outDir/samples/$sampleID/lancet"
  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from lanceting
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomebed) from Channel.value(params.exomebed)

  output:
  file('*.pass.vcf') into lancet_veping mode flatten
  file('*.raw.vcf') into lancet_rawVcf
  file('*') into completedlancetcall

  script:
  """
  TUMOURVCF=\$(echo $tumourbam | sed 's/bam/lancet.legacy.vcf/')
  {
    lancet \
      --num-threads ${task.cpus} \
      --ref $fasta \
      --bed $exomebed \
      --tumor $tumourbam \
      --normal $germlinebam | \
      perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
        \$_=~s/TUMOR/$sampleID/;
        \$_=~s/NORMAL/$germlineID/;
        print \$_;}
      else{print \$_;}' > \$TUMOURVCF

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=\$TUMOURVCF

  } 2>&1 | tee > $sampleID".lancet.log.txt"

  """
}

/* 3.0: Annotate Vcfs
*/
ALLVCFS = lancet_veping
          .mix( mutect2_veping )
          .mix( strelka2_veping )

process vepann {

  publishDir path: "$params.outDir/output/vcf", mode: "copy", pattern: '*.vcf'

  input:
  each file(vcf) from ALLVCFS
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(pcgr_grch_vep) from vepanndir

  output:
  file('*.vcf') into runGRanges

  script:
  """
  VCFANNO=\$(echo $vcf | sed "s/.vcf/.vep.vcf/")
  VEPVERS=\$(ls $pcgr_grch_vep/homo_sapiens/ | cut -d "_" -f2)
  vep --dir_cache $pcgr_grch_vep \
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
* each sample has 9 associated (raw,snv,indel per caller)
* NB increment if adding callers!
*/
ALLVVCFS = runGRanges
         .mix(lancet_rawVcf)
         .mix(mutect2_rawVcf)
         .mix(strelka2_rawVcf)

process vcfGRa {

  publishDir "$params.outDir/output/pdf", pattern: '*.pdf'
  publishDir "$params.outDir/output/vcf", pattern: '*.vcf'
  publishDir "$params.outDir/output/data", pattern: '*.[*RData,*tab]'

  input:
  file(grangesvcfs) from ALLVVCFS.collect()
  val(germlineID) from vcfGRaID

  output:
  file('*.ALL.pcgr.all.tab.vcf') into pcgrvcfs
  file('*') into completedvcfGRangesConsensus

  script:
  """
  Rscript --vanilla ${workflow.projectDir}/bin/variants_GRanges_consensus_plot.call.R \
    ${workflow.projectDir}/bin/variants_GRanges_consensus_plot.func.R \
    $germlineID \
    "snv_indel.pass.vep.vcf" \
    \"${params.includeOrder}\"
  """
}

process pcgrVcf {

  input:
  file(vcf) from pcgrvcfs

  output:
  file('*snv_indel.pass.pcgr.vcf') into pcgrinput

  script:
  """
  for VCF in *vcf; do
    NVCF=\$(echo \$VCF | sed 's/ALL.pcgr.all.tab.vcf/snv_indel.pass.pcgr.vcf/')
    cat ${workflow.projectDir}/bin/vcf42.head.txt > \$NVCF
    head -n1 \$VCF >> \$NVCF
    tail -n+2 \$VCF | sort -V >> \$NVCF
  done
  """
}

/* 3.2 PCGR report
* take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR
*/
process pcgrreport {

  publishDir "$params.outDir/reports", mode: "copy", pattern: "*html"
  publishDir "$params.outDir/output/pcgr", mode: "copy", pattern: "*[!.html]"

  input:
  file(pcgrinputs) from pcgrinput.collect()
  file(pcgrfacets) from facets_pcgr.collect()
  file(pcgr_grch) from PCGR
  file(pcgrmeta) from pcgrMetaInput
  val(grchvers) from pcgr_grchvers

  output:
  file('*') into completedPCGR

  script:
  """
  {
    for VCF in *vcf; do
      SAMPLEID=\$(echo \$VCF | cut -d "." -f 1)
      METAID=\$(grep \$SAMPLEID $pcgrmeta | cut -d "," -f 3 | sed 's/\\s */_/g')

      ##tuple up tumour ploidy (from facets)
      PLOIDY=\$(cut -f 1 \$SAMPLEID".fit_ploidy-purity.pcgr.tsv")
      PURITY=\$(cut -f 2 \$SAMPLEID".fit_ploidy-purity.pcgr.tsv")
      PP="--tumor_ploidy \$PLOIDY --tumor_purity \$PURITY"
      FACETI=\$SAMPLEID".cncf-jointsegs.pcgr.tsv"

      if [[ \$PLOIDY == "NA" || \$PURITY == "NA" ]]; then
        PP=""
      fi

      if [[ -e \$FACETI ]]; then
        FPP="--input_cna \$FACETI \$PP"
      else
        FPP="\$PP"
      fi

      CONFIG=\$(readlink -e pcgr/data/*/pcgr_configuration_default.toml)
      pcgr.py $pcgr_grch \
        ./ \
        $grchvers \
        \$CONFIG \
        \$METAID \
        --input_vcf \$VCF \$FPP \
        --no-docker \
        --force_overwrite
    done
  } 2>&1 | tee pcgr.log.txt
  """
}

/* 4.0 Run multiQC to finalise report
*/
process mltiQC {

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
  OUTID=\$(basename ${workflow.launchDir})
  if [[ -e ${params.multiqcConfig} ]]; then
    multiqc . -i \$OUTID --tag DNA -f -c ${params.multiqcConfig}
  else
    multiqc . -i \$OUTID --tag DNA -f
  fi
  """
}
