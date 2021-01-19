#!/usr/bin/env nextflow
def helpMessage() {
  log.info"""
  --------------------------------------------------------------------------
              HUMAN REFERENCE FOR SOMATIC_N-OF-1 NEXTFLOW PIPELINE
  --------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/somatic_n-of-1/download-references.nf

  Mandatory arguments:

    -profile        [str]     singularity,refs
    --assembly      [str]     GRCh37 or GRCh38 (default, NB GRCh37 GRIDSS will not  function as resources are unavailable and need to be made by user)
    --exomeTag      [str]     naming for exome outputs when supplied; tag is then used in somatic_n-of-1 and batch_somatic pipelines to select relevant exome data

    and either:
    --exomeBedURL   [str]     URL to exome bed file for intervals; NB assumes GRCh37
    or
    --exomeBedFile  [str]     locally downloaded exome bed file for intervals; NB assumes GRCh37

  General Optional Arguments:

    --cosmicUser    [str]     COSMIC login credentials, user email
    --cosmicPass    [str]     COSMIC login credentials, password
    """.stripIndet()
}

if (params.help) exit 0, helpMessage()

/* 0.0: Global Variables
*/
params.outDir = "${params.assembly}"

//base URL for GRCh37, 38
params.gsurl37 = "gs://gatk-legacy-bundles/b37"
params.gsurl38 = "gs://genomics-public-data/resources/broad/hg38/v0"

//lowercase assembly
params.assemblylc = "${params.assembly}".toLowerCase()

/* 1.0: Download GATK4 resource bundle fasta
*/

process fasta_gs {

  label 'gs'
  errorStrategy 'retry'
  maxRetries 3

  output:
  tuple val(faval), file('*.fasta.gz') into fasta

  script:
  faval = params.assemblylc == 'grch37' ? 'human_g1k_v37' : 'GRCh38_Verily_v1.genome'
  if( params.assemblylc == 'grch37' )
    """
    ##http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
    gsutil -q cp ${params.gsurl37}/human_g1k_v37.fasta.gz .
    """
  else
    """
    ##moved to Verily as gs bucket more reliable
    gsutil -q cp gs://genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa .
    gzip -c ./GRCh38_Verily_v1.genome.fa > ./GRCh38_Verily_v1.genome.fasta.gz
    rm ./GRCh38_Verily_v1.genome.fa
    """
}

process fasta_process {

  echo true
  label 'low_mem'

  input:
  tuple val(faval), file(fagz) from fasta

  output:
  tuple file('*.noChr.fasta'), file('*.noChr.fasta.fai') into (fasta_bwa, fasta_seqza, fasta_msi, fasta_dict, fasta_2bit, fasta_exome_biall, fasta_wgs_biall)

  script:
  def fa = "${fagz}".split("\\.gz")[0]
  """
  gunzip -c \$(readlink ${fagz}) > ${fa}
  cat ${fa} | sed 's/>chr/>/g' > ./${faval}.noChr.fasta
  rm ${fa}
  samtools faidx ./${faval}.noChr.fasta
  """
}

/* 1.1: Dictionary for fasta
*/
process dict_pr {

  label 'low_mem'
  publishDir path: "$params.outDir/bwa", mode: "copy"

  input:
  tuple file(fa), file(fai) from fasta_dict

  output:
  file('*.dict') into dict_win
  tuple file(fa), file(fai), file('*.dict') into (fasta_dict_exome, fasta_dict_wgs, fasta_dict_gensiz, fasta_dict_gridss)

  """
  DICTO=\$(echo $fa | sed 's/fasta/dict/')
  picard CreateSequenceDictionary \
    R=$fa \
    O=\$DICTO
  """
}

/* 1.2: Download GATK4 resource bundle dbsnp, Mills
*/
process dbsnp_dl {

  label 'gs'
  errorStrategy 'retry'
  maxRetries 3

  output:
  file('*.vcf') into vcf_tabix
  file('KG_phase1.snps.high_confidence.*.vcf') into ascatloci

  when:
  !params.nodbsnp

  script:
  if( params.assemblylc == 'grch37' )
    """
    gsutil -q cp ${params.gsurl37}/1000G_phase1.snps.high_confidence.b37.vcf.gz ./KG_phase1.snps.high_confidence.b37.vcf.gz
    gsutil -q cp ${params.gsurl37}/dbsnp_138.b37.vcf.gz ./dbsnp_138.b37.vcf.gz
    gsutil -q cp ${params.gsurl37}/hapmap_3.3.b37.vcf.gz ./hapmap_3.3.b37.vcf.gz
    gsutil -q cp ${params.gsurl37}/1000G_omni2.5.b37.vcf.gz ./KG_omni2.5.b37.vcf.gz
    gsutil -q cp ${params.gsurl37}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz ./Mills_KG_gold.indels.b37.vcf.gz

    gunzip -cd dbsnp_138.b37.vcf.gz | sed 's/chr//g' > dbsnp_138.b37.vcf
    gunzip -cd hapmap_3.3.b37.vcf.gz | sed 's/chr//g' > hapmap_3.3.b37.sites.vcf
    gunzip -cd KG_omni2.5.b37.vcf.gz | sed 's/chr//g' > KG_omni2.5.b37.vcf
    gunzip -cd KG_phase1.snps.high_confidence.b37.vcf.gz | sed 's/chr//g' > KG_phase1.snps.high_confidence.b37.vcf
    gunzip -cd Mills_KG_gold.indels.b37.vcf.gz | sed 's/chr//g' > Mills_KG_gold.indels.b37.vcf
    """
  else
    """
    gsutil -q cp gs://genomics-public-data/cwl-examples/gdc-dnaseq-cwl/input/dbsnp_144.hg38.vcf.gz ./dbsnp_144.hg38.vcf.gz
    gsutil -q cp ${params.gsurl38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz ./KG_phase1.snps.high_confidence.hg38.vcf.gz
    gsutil -q cp ${params.gsurl38}/Homo_sapiens_assembly38.dbsnp138.vcf ./Homo_sapiens_assembly38.dbsnp138.vcf
    gsutil -q cp ${params.gsurl38}/hapmap_3.3.hg38.vcf.gz ./hapmap_3.3.hg38.vcf.gz
    gsutil -q cp ${params.gsurl38}/1000G_omni2.5.hg38.vcf.gz ./KG_omni2.5.hg38.vcf.gz
    gsutil -q cp ${params.gsurl38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ./Mills_KG_gold.indels.hg38.vcf.gz

    gunzip -cd dbsnp_144.hg38.vcf.gz | sed 's/chr//g' > dbsnp_144.hg38.vcf
    gunzip -cd hapmap_3.3.hg38.vcf.gz | sed 's/chr//g' > hapmap_3.3.hg38.vcf
    gunzip -cd KG_omni2.5.hg38.vcf.gz | sed 's/chr//g' > KG_omni2.5.hg38.vcf
    gunzip -cd KG_phase1.snps.high_confidence.hg38.vcf.gz | sed 's/chr//g' > KG_phase1.snps.high_confidence.hg38.vcf
    gunzip -cd Mills_KG_gold.indels.hg38.vcf.gz | sed 's/chr//g' > Mills_KG_gold.indels.hg38.vcf
    """
}

/* 1.3: KG ASCAT loci
*/
process ascat_loci {

  label 'low_mem'
  publishDir path: "$params.outDir", mode: "copy"

  input:
  file(vcf) from ascatloci

  output:
  file('*loci') into complete_ascat

  when:
  params.legacy

  script:
  """
  LOCIFILE=\$(echo $vcf | sed 's/vcf/maf0.3.loci/')
  cat $vcf | \
  perl -ane '@s=split(/[=\\;]/,\$F[7]);if(\$s[3]>0.3){print "\$F[0]\\t\$F[1]\\n";}' > \$LOCIFILE
  """
}

/* 2.0: Fasta processing
*/
process bwa_index {

  label 'med_mem'
  publishDir path: "$params.outDir/bwa", mode: "copy"

  input:
  tuple file(fa), file(fai) from fasta_bwa

  output:
  file('*') into complete_bwa

  script:
  """
  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
  bwa index -a bwtsw $fa
  """
}

/* 2.1: Dict processing
*/
process dict_pr2 {

  label 'low_mem'
  publishDir path: "$params.outDir", mode: "copy"

  input:
  file(win_dict) from dict_win

  output:
  file('*') into complete_dict

  when:
  params.legacy

  script:
  """
  perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' $win_dict > seq.dict.chr-size

  bedtools makewindows -g seq.dict.chr-size -w 35000000 | perl -ane 'if(\$F[1]==0){\$F[1]++;};print "\$F[0]:\$F[1]-\$F[2]\n";' > 35MB-window.bed
  """
}

/* 3.0: Exome bed file and liftOver
*/
if(params.exomeBedFile && params.exomeBedURL){
  Channel.from("Please only specify --exomeBedURL or --exomeBedFile!\nN.B. that subsquent runs using -resume can be used to add further -exomeBedURL or --exomeBedFile").println { it }
}
if(params.exomeBedFile && params.exomeBedURL){
  Channel.from("Only one of --exomeBedFile or --exomeBedURL should be used").println { it }
  exit 147
}

//set exomeTag
if(!params.exomeTag){
  if(params.exomeBedURL) {
    params.exomeTag = "${params.exomeBedURL}".split("\\.")[0]
  }
  if(params.exomeBedFile) {
    params.exomeTag = "${params.exomeBedFile}".split("\\.")[0]
  }
}

if(params.exomeBedURL){
  process exome_url {

    label 'low_mem'
    publishDir path: "$params.outDir/exome/$params.exomeTag", mode: "copy"

    output:
    tuple file("${params.exomeTag}.url.bed"), file("README.${params.exomeTag}.url.bed") into exome_bed

    script:
    """
    ##download URL
    echo "Exome bed used here is from:" > README.${params.exomeTag}.url.bed
    echo ${params.exomeBedURL} >> README.${params.exomeTag}.url.bed

    wget ${params.exomeBedURL}
    if [[ ${params.exomeBedURL} =~ zip\$ ]]; then
      unzip -p *.zip > ${params.exomeTag}.url.bed
    elif [[ ${params.exomeBedURL} =~ bed\$ ]]; then

      ##remove any non-chr, coord lines in top of file
      CHR=\$(tail -n1 ${params.exomeBedURL} | perl -ane 'print \$F[0];')
      if [[ \$CHR =~ "chr" ]]; then
        perl -ane 'if(\$F[0]=~m/^chr/){print \$_;}' ${params.exomeBedURL} >  ${params.exomeTag}.url.bed
      else
        perl -ane 'if(\$F[0]=~m/^[0-9MXY]/){print \$_;}' ${params.exomeBedURL} >  ${params.exomeTag}.url.bed
      fi

    else
      echo "No ZIP or BED files resulting from ${params.exomeBedURL}"
      echo "Please try another URL with ZIP or BED file resulting"
      exit 147
    fi
    """
  }
}

if(params.exomeBedFile){
  Channel.fromPath("${params.exomeBedFile}").set { exomebed_file }
  process exome_file {

    label 'low_mem'
    publishDir path: "$params.outDir/exome/$params.exomeTag", mode: "copy"

    input:
    file(exomeBedFile) from exomebed_file

    output:
    tuple file("${params.exomeTag}.file.bed"), file("README.${params.exomeTag}.file.bed") into exome_bed

    script:
    """
    ##use file as input
    echo "Exome bed used here is from:" > README.${params.exomeTag}.file.bed
    echo $exomeBedFile >> README.${params.exomeTag}.file.bed

    if [[ $exomeBedFile =~ bed\$ ]]; then

      ##remove any non-chr, coord lines in top of file
      CHR=\$(tail -n1 $exomeBedFile | perl -ane 'print \$F[0];')
      if [[ \$CHR =~ "chr" ]]; then
        perl -ane 'if(\$F[0]=~m/^chr/){print \$_;}' $exomeBedFile >  ${params.exomeTag}.file.bed
      else
        perl -ane 'if(\$F[0]=~m/^[0-9MXY]/){print \$_;}' $exomeBedFile >  ${params.exomeTag}.file.bed
      fi

    else
      echo "BED file $exomeBedFile is not a BED file, please retry"
      exit 147
    fi
    """
  }
}

process lift_over {

  label 'low_mem'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple file(exomebed), file(readme) from exome_bed

  output:
  file('*.lift.bed') into exome_bed_liftd

  script:
  if( params.assembly == "GRCh37" )
    """
    cp $exomebed ${params.exomeTag}.lift.bed
    """
  else
    """
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    liftOver $exomebed hg19ToHg38.over.chain.gz ${params.exomeTag}.lift.bed unmapped
    """
}

/* 3.11: Parse bed for exome
*/
process exome_bed_pr {

  label 'low_mem'
  publishDir path: "$params.outDir/exome/$params.exomeTag", mode: "copy", pattern: "*[.interval_list,.bed]"

  input:
  tuple file(fa), file(fai), file(dict) from fasta_dict_exome
  file(exomelift) from exome_bed_liftd

  output:
  file("${params.exomeTag}.bed.interval_list") into complete_exome
  file("${params.exomeTag}.bed") into (exome_tabix, exome_biallgz)

  script:
  """
  ##must test if all chr in fasta are in exome, else manta cries
  ##must test if all regions are greater than length zero or strelka cries
  ##must test if all seq.dict chrs are in bed and only they or BedToIntervalList cries
  perl -ane 'if(\$F[1] == \$F[2]){\$F[2]++;} if(\$F[0] !~m/^chrM/){print join("\\t", @F[0..\$#F]) . "\\n";}' $exomelift | grep -v chrM | sed 's/chr//g' > tmp.bed

   grep @SQ $dict | cut -f2 | sed 's/SN://' | while read CHR; do
   TESTCHR=\$(awk -v chrs=\$CHR '\$1 == chrs' tmp.bed | wc -l)
   if [[ \$TESTCHR != 0 ]];then
    awk -v chrs=\$CHR '\$1 == chrs' tmp.bed
   fi
  done >> tmp.dict.bed

  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=tmp.dict.bed O=${params.exomeTag}.interval_list SD=$dict

  ##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
  perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' ${params.exomeTag}.interval_list > ${params.exomeTag}.bed.interval_list

  ##output BED
  grep -v "@" ${params.exomeTag}.bed.interval_list | cut -f 1,2,3,5 > ${params.exomeTag}.bed
  """
}

/* 3.12: create bed for WGS
*/
process wgs_bed {

  label 'low_mem'
  publishDir path: "$params.outDir/wgs", mode: "copy"

  input:
  tuple file(fa), file(fai), file(dict) from fasta_dict_wgs

  output:
  file('wgs.bed.interval_list') into complete_wgs
  file('wgs.bed') into (wgs_tabix, wgs_fasta_biallgz)

  script:
  """
  ##WGS intervals = 1-LN for each chr
  grep @SQ $dict | cut -f 2,3 | perl -ane '\$chr=\$F[0];\$chr=~s/SN://;\$end=\$F[1];\$end=~s/LN://;print "\$chr\\t0\\t\$end\\n";' > tmp.wgs.dict.bed

  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=tmp.wgs.dict.bed O=wgs.bed.interval_list SD=$dict

  ##output BED
  grep -v "@" wgs.bed.interval_list | cut -f 1,2,3,5 > wgs.bed
  """
}

/* 3.2: Tabix those requiring tabixing
*/
wgs_tabix.concat(exome_tabix).set { bint_tabix }
process tabix_files {

  label 'low_mem'
  publishDir path: "$params.outDir/exome/$params.exomeTag", mode: "copy", pattern: "{${params.exomeTag}}*"
  publishDir path: "$params.outDir/wgs", mode: "copy", pattern: "{wgs}*"

  input:
  file(bed) from bint_tabix

  output:
  tuple file("${bed}.gz"), file("${bed}.gz.tbi") into complete_tabix

  script:
  """
  ##tabix
  bgzip $bed
  tabix $bed".gz"
  """
}


process gnomads {

  label 'gs'

  output:
  file('af-only-gnomad.*') into ( exome_biall_gnomad, wgs_biall_gnomad )

  script:
  if( params.assembly == "GRCh37" )
  """
  gsutil -q cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ./
  """
  else
  """
  gsutil -q cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz ./
  """
}

/* 3.31: Create Mutect2 af-only-gnomad file
*/
process exome_biall {

  label 'low_mem'
  publishDir path: "$params.outDir/exome/$params.exomeTag", mode: "copy", pattern: "![exome.biall.bed, af-only-gnomad.exomerh.hg38.noChr.vcf]"

  input:
  file(exomebed) from exome_biallgz
  file(gnomad) from exome_biall_gnomad
  tuple file(fasta), file(fai) from fasta_exome_biall

  output:
  tuple file("af-only-gnomad.${params.exomeTag}.hg*.noChr.vcf.gz"), file("af-only-gnomad.${params.exomeTag}.hg*.noChr.vcf.gz.tbi") into exome_biallelicgz
  file('exome.biall.bed') into pcgrtoml_exome

  script:
  if( params.assembly == "GRCh37" )
    """
    cut -f 1,2,3 $exomebed > exome.biall.bed
    bgzip ${gnomad}
    tabix ${gnomad}.gz
    gunzip -c ${gnomad}.gz |
    bcftools view -R exome.biall.bed ${gnomad}.gz | bcftools sort -T '.' > af-only-gnomad.exomerh.hg19.noChr.vcf
    perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.exomerh.hg19.noChr.vcf $fai > af-only-gnomad.${params.exomeTag}.hg19.noChr.vcf
    bgzip af-only-gnomad.${params.exomeTag}.hg19.noChr.vcf
    tabix af-only-gnomad.${params.exomeTag}.hg19.noChr.vcf.gz
    """
  else
    """
    cut -f 1,2,3 $exomebed > exome.biall.bed
    gunzip -c ${gnomad} | sed 's/chr//' | bgzip > af-only-gnomad.hg38.noChr.vcf.gz
    tabix af-only-gnomad.hg38.noChr.vcf.gz
    bcftools view -R exome.biall.bed af-only-gnomad.hg38.noChr.vcf.gz | bcftools sort -T '.' > af-only-gnomad.exomerh.hg38.noChr.vcf
    perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.exomerh.hg38.noChr.vcf $fai > af-only-gnomad.${params.exomeTag}.hg38.noChr.vcf
    bgzip af-only-gnomad.${params.exomeTag}.hg38.noChr.vcf
    tabix af-only-gnomad.${params.exomeTag}.hg38.noChr.vcf.gz
    """
}

/* 3.32: Create Mutect2 af-only-gnomad file
*/
process wgs_biall {

  label 'low_mem'
  publishDir path: "$params.outDir/wgs", mode: "copy"
  errorStrategy 'retry'
  maxRetries 3

  input:
  file(wgsbed) from wgs_fasta_biallgz
  file(gnomad) from wgs_biall_gnomad
  tuple file(fasta), file(fai) from fasta_wgs_biall

  output:
  tuple file('af-only-gnomad.wgs.*.noChr.vcf.gz'), file('af-only-gnomad.wgs.*.noChr.vcf.gz.tbi') into wgs_biallelicgz
  file('wgs.biall.bed') into pcgrtoml_wgs

  script:
  if( params.assembly == "GRCh37" )
    """
    cut -f 1,2,3 ${wgsbed} > wgs.biall.bed
    bgzip ${gnomad}
    tabix ${gnomad}.gz
    bcftools view -R wgs.biall.bed ${gnomad}.gz | bcftools sort -T '.' > af-only-gnomad.wgsh.hg19.noChr.vcf
    perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.wgsh.hg19.noChr.vcf $fai > af-only-gnomad.wgs.hg19.noChr.vcf
    bgzip af-only-gnomad.wgs.hg19.noChr.vcf
    tabix af-only-gnomad.wgs.hg19.noChr.vcf.gz
    """
  else
    """
    cut -f 1,2,3 $wgsbed > wgs.biall.bed
    gunzip -c ${gnomad} | sed 's/chr//' | bgzip > af-only-gnomad.hg38.noChr.vcf.gz
    tabix af-only-gnomad.hg38.noChr.vcf.gz
    bcftools view -R wgs.biall.bed af-only-gnomad.hg38.noChr.vcf.gz | bcftools sort -T '.' > af-only-gnomad.wgsh.hg38.noChr.vcf
    perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.wgsh.hg38.noChr.vcf $fai > af-only-gnomad.wgs.hg38.noChr.vcf
    bgzip af-only-gnomad.wgs.hg38.noChr.vcf
    tabix af-only-gnomad.wgs.hg38.noChr.vcf.gz
    """
}

/* 4.0 Index various VCFs
*/
process indexfeature_files {

  label 'low_mem'
  publishDir path: "$params.outDir/hc_dbs", mode: "copy", pattern: "{KG,Mills,hapmap}*"
  publishDir path: "$params.outDir/dbsnp", mode: "copy", pattern: "{dbsnp}*"
  input:
  file(tbtbx) from vcf_tabix.flatten()

  output:
  file('*') into indexfeatured

  script:
  """
  bgzip $tbtbx
  gatk IndexFeatureFile --input $tbtbx".gz"
  """
}

/* 5.0: Sequenza GC bias
*/
process seqnza {

  label 'low_mem'
  publishDir path: "$params.outDir", mode: "copy"

  input:
  set file(fa), file(fai) from fasta_seqza

  output:
  file('*') into sequenzaout

  when:
  params.legacy

  script:
  """
  GENOMEGC50GZ=\$(echo $fa | sed -r 's/.fasta/.gc50Base.txt.gz/')
  sequenza−utils.py GC-windows −w 50 $fa | gzip > \$GENOMEGC50GZ
  """
}

/* 6.0: MSIsensor microsatellites
*/
process msisen {

  label 'low_mem'
  publishDir "$params.outDir", mode: "copy"

  input:
  set file(fa), file(fai) from fasta_msi

  output:
  file('*') into completedmsisensor

  when:
  params.legacy

  script:
  """
  msisensor scan -d $fa -o msisensor_microsatellites.list
  """
}

/* 7.0: PCGR/CPSR data bundle
*/
process pcgr_vep {

  label 'low_mem'
  publishDir "$params.outDir/pcgr", mode: "copy", pattern: "data"
  errorStrategy 'retry'
  maxRetries 3

  input:
  file(exomebed) from pcgrtoml_exome
  file(wgsbed) from pcgrtoml_wgs

  output:
  file('data') into pcgrdata

  when:
  !params.nopcgr

  script:
  downloadURL = params.assembly == "GRCh37" ? "${params.pcgrURL37}" : "${params.pcgrURL38}"
  exometag = params.exomeTag == null ? "" : "${params.exomeTag}"
  """
  wget ${downloadURL}
  tar -xf *.tgz
  rm -rf *.tgz

  ##allows editting of TOML for WGS and exome if supplied
  sh ${workflow.projectDir}/bin/pcgr_edit_toml.sh \
    data/${params.assemblylc}/pcgr_configuration_default.toml \
    ${exomebed} ${exometag}
  mv pcgr_configuration*.toml data/${params.assemblylc}/

  ##build VEP cache using Singularity container 'vep_install' script
  ##however PCGR installs a version of vep cache, so test that matches required assembly, and only install if not
  ##required version is the VEP_INSTALL version in container

  ##variables for install and test
  VEP_INSTALL=\$(find /opt/miniconda/envs/somatic_n-of-1/share/*/vep_install)
  VEP_INSTVER=\$(echo \$VEP_INSTALL | perl -ane '@s=split(/\\//, \$F[0]); foreach \$k (@s){ if(\$k =~m/ensembl-vep/){@o=split(/[-.]/, \$k); print \$o[2];}}')
  VEP_PCGRVER=\$(cat data/${params.assemblylc}/RELEASE_NOTES | perl -ane 'if(\$F[0] eq "VEP"){@s=split(/\\./,\$F[5]); \$v=\$s[0]; \$v=~s/v//; print \$v;}')

  if [[ \$VEP_INSTVER != \$VEP_PCGRVER ]];then

    ##remove current and reinstall correct version
    \$VEP_INSTALL \
      --AUTO cf \
      --CACHE_VERSION \$VEP_INSTVER \
      --CACHEDIR "data/${params.assemblylc}/.vep/" \
      --SPECIES "homo_sapiens" \
      --ASSEMBLY ${params.assembly} \
      --NO_UPDATE \
      --NO_HTSLIB \
      --NO_BIOPERL \
      --NO_TEST
  fi
  """
}

/* 8.0: GenomeSize.xml for Pisces
*/
process gensizxml {

  label 'low_mem'
  publishDir "$params.outDir", mode: "copy"

  input:
  set file(fa), file(fai), file(dict) from fasta_dict_gensiz

  output:
  file('*') into complete_gensiz

  when:
  params.legacy

  script:
  """
  echo "<sequenceSizes genomeName=\"$dict\">" > GenomeSize.xml
  grep "@SQ" $dict | while read LINE; do
    CONTIGNAME=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[1]);print \$s[1];' | sed 's/chr//')
    TOTALBASES=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[2]);print \$s[1];')
    MD5SUM=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[3]);print \$s[1];')
    echo -e "\\t<chromosome fileName=\"$fa\" contigName=\"\$CONTIGNAME\" totalBases=\"\$TOTALBASES\" isCircular=\"false\" md5=\"\$MD5SUM\" ploidy=\"2\" knownBases=\"\$TOTALBASES\" type=\"Chromosome\" />" >> GenomeSize.xml
  done
  echo "</sequenceSizes>" >> GenomeSize.xml
  """
}

/* 9.0: Download hartwigmedical resource bundle
*/
process hartwigmed {

  label 'low_mem'
  publishDir path: "$params.outDir/gridss", mode: "copy"
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple file(fa), file(fai), file(dict) from fasta_dict_gridss

  output:
  file('dbs') into gridss_db
  file('refgenomes/human_virus') into gridss_hv
  file('gridss_blacklist.noChr.bed') into gridss_bl
  file('dbs/gridss/gridss.properties') into gridss_pr

  when:
  !params.exomeTag

  script:
  if( params.assembly == "GRCh37" )
    """
    curl -o gpl_ref_data_hg37.tar.gz "${params.hartwigGPLURL37}"
    tar -xf gpl_ref_data_hg37.tar.gz

    ##blacklist
    sed 's/chr//g' dbs/gridss/ENCFF001TDO.bed > gridss_blacklist.noChr.bed

    perl ${workflow.projectDir}/bin/exact_match_by_col.pl ${fai},0 gridss_blacklist.noChr.bed,0 > gridss_blacklist.noChr.1.bed
    mv gridss_blacklist.noChr.1.bed gridss_blacklist.noChr.bed
    """
  else
    """
    curl -o gpl_ref_data_hg38.tar.gz "${params.hartwigGPLURL38}"
    tar -xf gpl_ref_data_hg38.tar.gz
    rm -rf gpl_ref_data_hg38.tar.gz

    ##blacklist
    sed 's/chr//g' dbs/gridss/ENCFF001TDO.bed > gridss_blacklist.noChr.bed
    perl ${workflow.projectDir}/bin/exact_match_by_col.pl ${fai},0 gridss_blacklist.noChr.bed,0 > gridss_blacklist.noChr.1.bed
    mv gridss_blacklist.noChr.1.bed gridss_blacklist.noChr.bed
    """
}

process cosmic {

  publishDir "${params.outDir}/cosmic", mode: 'copy'

  when:
  params.cosmicUser && params.cosmicPass

  output:
  file 'cancer_gene_census.*'

  script:
  """
  ##README to show date of download
  echo \$(date) > cosmic_cancer-gene-census_dl.README.txt

  ##https://cancer.sanger.ac.uk/cosmic/help/file_download
  curl -H "Authorization: Basic \$(echo ${params.cosmicUser}:${params.cosmicPass} | base64)" https://cancer.sanger.ac.uk/cosmic/file_download/${params.assembly}/cosmic/${params.cosmicVers}/cancer_gene_census.csv > url.txt

  URL=\$(cut -d '"' -f4 url.txt)
  curl -o cancer_gene_census.csv \$URL

  ##parse into BED format (many thanks for comma in 'Name' field NOT)
  tail -n+2 cancer_gene_census.csv | \\
  perl -ane '@s=split(/\\,/);
             foreach \$k (@s){
               if(\$k=~m/\\d+:\\d+-\\d+/){
                   @p=split(/[:-]/, \$k);
                   if(\$s[-2]=~/^ENS/){
                     \$ens=\$s[-2];
                   }
                   if(\$s[-3]=~/^ENS/){
                     \$ens=\$s[-3];
                   }
               print "\$p[0]\\t\$p[1]\\t\$p[2]\\t\$s[0];\$ens\\n"; next;
             }};' | sort -V > cancer_gene_census.bed

  """
}
