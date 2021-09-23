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
    --assembly      [str]     GRCh37 or GRCh38
    --exomeTag      [str]     naming for exome outputs when supplied;
                              tag is then used in somatic_n-of-1 and
                              batch_somatic pipelines to select relevant exome
                              data (default: the first element when exome file
                              name is split on '.')

    one of:
    --exomeBedURL   [str]     URL of exome BED file (ZIP or BED)
    --exomeBedFile  [str]     path to local exome BED file

  General Optional Arguments:

    --exomeAssembly [str]     assembly the exome BED file used in it's creation
                              (default: GRCh37)
    --localPCGRdata [str]     local copy of PCGR data bundle to use
    --cosmicUser    [str]     COSMIC login credentials, user email
    --cosmicPass    [str]     COSMIC login credentials, password
    --microbiome    [bool]    Builds reference for Pathseq microbiome analysis (NB: massive output (100's Gb))
    """.stripIndet()
}

/*
================================================================================
                          0. INITIAL SETUP AND VARIABLES
================================================================================
*/

if (params.help) exit 0, helpMessage()

/* 0.0: Global Variables
*/
params.outDir = "${params.assembly}"

//base URL for GRCh37, 38
params.gsurl37 = "gs://gatk-legacy-bundles/b37"
params.gsurl38 = "gs://genomics-public-data/resources/broad/hg38/v0"

//lowercase assembly
params.assemblylc = "${params.assembly}".toLowerCase()

//specify single exome
if(params.exomeBedFile && params.exomeBedURL){
  Channel.from("Please only specify --exomeBedURL or --exomeBedFile!\nN.B. that subsquent runs using -resume can be used to add further -exomeBedURL or --exomeBedFile").println { it }
}
if(params.exomeBedFile && params.exomeBedURL){
  Channel.from("Only one of --exomeBedFile or --exomeBedURL should be used").println { it }
  exit 147
}

//make exomeTag if not extant but exomeBedFile or URL are specified
if(!params.exomeTag){
  if(params.exomeBedURL) {
    params.exomeTag = "${params.exomeBedURL}".split("\\.")[0]
  }
  if(params.exomeBedFile) {
    params.exomeTag = "${params.exomeBedFile}".split("\\.")[0]
  }
}

//param to specify if PCGR data exists yet
params.pcgrdata = false
if(!file("$params.outDir/pcgr/data/${params.assemblylc}").exists()){
  params.pcgrdata = true
}

/*
================================================================================
                          1. FASTA REFERENCES
================================================================================
*/

// 1.0: Download GATK4 resource bundle fasta if not extant
if(!file("$params.outDir/bwa").exists()){

  process fasta_gs {

    label 'gs'
    errorStrategy 'retry'
    maxRetries 3

    output:
    tuple val(faval), file('*.fasta.gz') into fasta

    script:
    def faval = params.assemblylc == 'grch37' ? 'human_g1k_v37' : 'GRCh38_Verily_v1.genome'
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
    tuple file('*.noChr.fasta'), file('*.noChr.fasta.fai') into (fasta_bwa, fasta_seqza, fasta_msi, fasta_dict, fasta_2bit, fasta_exome_biall, fasta_wgs_biall, fasta_pathseq)

    script:
    def fa = "${fagz}".split("\\.gz")[0]
    """
    gunzip -c \$(readlink ${fagz}) > ${fa}
    cat ${fa} | sed 's/>chr/>/g' > ./${faval}.noChr.fasta
    rm ${fa}
    samtools faidx ./${faval}.noChr.fasta
    """
  }

  process dict_pr {

    label 'low_mem'
    publishDir path: "${params.outDir}/bwa", mode: "copy"

    input:
    tuple file(fa), file(fai) from fasta_dict

    output:
    file(dict) into dict_win
    tuple file(fa), file(fai), file(dict) into ( fasta_dict_exome, fasta_dict_wgs, fasta_dict_gensiz )
    file(fai) into fai_gridss

    script:
    def dict = "${fa}".replace('fasta', 'dict')
    """
    picard CreateSequenceDictionary \
      R=${fa} \
      O=${dict}
    """
  }

  process bwa_index {

    label 'med_mem'
    publishDir path: "${params.outDir}/bwa", mode: "copy"

    input:
    tuple file(fa), file(fai) from fasta_bwa

    output:
    file('*') into complete_bwa

    script:
    """
    ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
    bwa index -a bwtsw ${fa}
    """
  }
}

// 1.1: send fasta etc. refs
if(file("$params.outDir/bwa").exists()){

  Channel.fromPath("${params.outDir}/bwa").set { bwa_chan }
  process send_dict_pr {

    label 'low_mem'

    input:
    file(bwa_dir) from bwa_chan

    output:
    file('bwa/*.dict') into dict_win
    tuple file('bwa/*noChr.fasta'), file('bwa/*.fai'), file('bwa/*.dict') into (fasta_dict_exome, fasta_dict_wgs, fasta_dict_gensiz)
    tuple file('bwa/*.noChr.fasta'), file('bwa/*.noChr.fasta.fai') into (fasta_bwa, fasta_seqza, fasta_msi, fasta_dict, fasta_2bit, fasta_exome_biall, fasta_wgs_biall, fasta_pathseq)
    file('bwa/*.noChr.fasta.fai') into fai_gridss

    script:
    """
    """
  }
}

/*
================================================================================
                          2. VCF REFERENCES
================================================================================
*/

// 2.0: Download dbsnp and other GATK4 reference bundle VCFs
if(!file("$params.outDir/dbsnp").exists()){

  process vcf_dl {

    label 'gs'
    errorStrategy 'retry'
    maxRetries 3

    output:
    file('*.vcf') into vcf_tabix

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

  process index_feature_files {

    label 'low_mem'
    publishDir path: "${params.outDir}/hc_dbs", mode: "copy", pattern: "{KG,Mills,hapmap}*"
    publishDir path: "${params.outDir}/dbsnp", mode: "copy", pattern: "dbsnp*"

    input:
    file(tbtbx) from vcf_tabix.flatten()

    output:
    file('*') into indexfeatured

    script:
    """
    bgzip ${tbtbx}
    gatk IndexFeatureFile --input ${tbtbx}".gz"
    """
  }
}

// 2.1: Download gnomad
if(!file("$params.outDir/gnomad").exists()){

  process gnomad_dl {

    label 'gs'
    publishDir path: "${params.outDir}/gnomad", mode: "copy"

    output:
    file('af-only-gnomad.*') into ( exome_biall_gnomad, wgs_biall_gnomad )

    script:
    if( params.assembly == "GRCh37" )
      """
      gsutil -q cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ./
      bgzip af-only-gnomad.raw.sites.vcf
      tabix af-only-gnomad.raw.sites.vcf.gz
      """
    else
      """
      gsutil -q cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz ./
      """
  }
}

// 2.2: else send gnomad
if(file("$params.outDir/gnomad").exists()){

  Channel.fromPath("${params.outDir}/gnomad").set { gnomad_chan }
  process send_gnomads {

    input:
    file(gnomad_dir) from gnomad_chan

    output:
    file('gnomad/af-only-gnomad.*') into ( exome_biall_gnomad, wgs_biall_gnomad )

    script:
    """
    """
  }
}

/*
================================================================================
                          3. EXOME PROCESS AND PARSING
================================================================================
*/

if(!file("$params.outDir/exome/$params.exomeTag").exists()){
  if(params.exomeBedURL){
    process exome_url {

      label 'low_mem'

      output:
      tuple file("${params.exomeTag}.url.bed"), file("README.${params.exomeTag}.url.bed") into exome_parse_bed

      script:
      """
      ##download URL
      echo "Exome bed used here is from:" > README.${params.exomeTag}.url.bed
      echo ${params.exomeBedURL} >> README.${params.exomeTag}.url.bed

      wget ${params.exomeBedURL}
      if [[ ${params.exomeBedURL} =~ zip\$ ]]; then
        7za x -so *.zip  > ${params.exomeTag}.url.bed
      fi
      if [[ ${params.exomeBedURL} =~ bed\$ ]]; then
        mv *.bed ${params.exomeTag}.url.bed
      fi
      if [[ ! ${params.exomeBedURL} =~ zip\$ || ${params.exomeBedURL} =~ bed\$ ]]; then
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
      publishDir path: "${params.outDir}/exome/${params.exomeTag}", mode: "copy"

      input:
      file(exome_bed_file) from exomebed_file

      output:
      tuple file("${params.exomeTag}.file.bed"), file("README.${params.exomeTag}.file.bed") into exome_parse_bed

      script:
      """
      ##use file as input
      if [[ ${exome_bed_file} =~ bed\$ ]]; then
        echo "Exome bed used here is from:" > README.${params.exomeTag}.file.bed
        echo ${exome_bed_file} >> README.${params.exomeTag}.file.bed
        mv ${exome_bed_file} ${params.exomeTag}.file.bed
      else
        echo "BED file ${exome_bed_file} is not a BED file, please retry"
        exit 147
      fi
      """
    }
  }

  process exome_parse {

    input:
    tuple file(exome_bed), file(readme) from exome_parse_bed

    output:
    tuple file("${params.exomeTag}.lo.bed"), file(readme) into exome_liftover

    script:
    """
    ##remove any non-chr, coord lines in top of file
    CHR=\$(tail -n1 ${exome_bed} | perl -ane 'print \$F[0];')
    if [[ \$CHR =~ "chr" ]]; then
      perl -ane 'if(\$F[0]=~m/^chr/){print \$_;}' ${exome_bed} > ${params.exomeTag}.lo.bed
    else
      perl -ane 'if(\$F[0]=~m/^[0-9MXY]/){print \$_;}' ${exome_bed} > ${params.exomeTag}.lo.bed
    fi
    """
  }

  process lift_over {

    label 'low_mem'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple file(exome_bed), file(readme) from exome_liftover

    output:
    tuple file("${params.exomeTag}.lift.bed"), file(readme) into exome_bed_liftd

    script:
    def hgTohg = params.exomeAssembly == "GRCh38" ? "hg38ToHg19" : "hg19ToHg38"
    def hg = params.exomeAssembly == "GRCh38" ? "hg38" : "hg19"

    if( params.assembly == params.exomeAssembly )
      """
      cp ${exome_bed} ${params.exomeTag}.lift.bed
      """
    else
      """
      wget http://hgdownload.cse.ucsc.edu/goldenPath/${hg}/liftOver/${hgTohg}.over.chain.gz
      liftOver ${exome_bed} ${hgTohg}.over.chain.gz ${params.exomeTag}.lift.bed unmapped

      echo -r"Liftover using:\\nwget http://hgdownload.cse.ucsc.edu/goldenPath/${hg}/liftOver/${hgTohg}.over.chain.gz
      liftOver ${exome_bed} ${hgTohg}.over.chain.gz ${params.exomeTag}.lift.bed unmapped" >> ${readme}
      """
  }

  process exome_bed_pr {

    label 'low_mem'
    publishDir path: "${params.outDir}/exome/${params.exomeTag}", mode: "copy", overwrite: "true"

    input:
    tuple file(fa), file(fai), file(dict) from fasta_dict_exome
    tuple file(exome_lift), file(readme) from exome_bed_liftd

    output:
    file("${params.exomeTag}.bed.interval_list") into complete_exome
    file("${params.exomeTag}.bed") into ( exome_biallgz, pcgrtoml_exome )
    tuple file("${params.exomeTag}.bed.gz"), file("${params.exomeTag}.bed.gz.tbi") into gztbi_exome

    script:
    """
    ##must test if all chr in fasta are in exome, else manta cries
    ##must test if all regions are greater than length zero or strelka cries
    ##must test if all seq.dict chrs are in bed and only they or BedToIntervalList cries
    perl -ane 'if(\$F[1] == \$F[2]){\$F[2]++;} if(\$F[0] !~m/^chrM/){print join("\\t", @F[0..\$#F]) . "\\n";}' ${exome_lift} | grep -v chrM | sed 's/chr//g' > tmp.bed

     grep @SQ ${dict} | cut -f2 | sed 's/SN://' | while read CHR; do
     TESTCHR=\$(awk -v chrs=\$CHR '\$1 == chrs' tmp.bed | wc -l)
     if [[ \$TESTCHR != 0 ]];then
      awk -v chrs=\$CHR '\$1 == chrs' tmp.bed
     fi
    done >> tmp.dict.bed

    ##always make interval list so we are in line with fasta
    picard BedToIntervalList I=tmp.dict.bed O=${params.exomeTag}.interval_list SD=${dict}

    ##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
    perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' ${params.exomeTag}.interval_list > ${params.exomeTag}.bed.interval_list

    ##output BED
    grep -v "@" ${params.exomeTag}.bed.interval_list | cut -f 1,2,3,5 > ${params.exomeTag}.bed

    ##tabix
    bgzip -c ${params.exomeTag}.bed > ${params.exomeTag}.bed.gz
    tabix ${params.exomeTag}.bed.gz
    """
  }

  process exome_biall {

    label 'low_mem'
    publishDir path: "${params.outDir}/exome/${params.exomeTag}", mode: "copy"

    input:
    file(exome_bed) from exome_biallgz
    file(gnomad) from exome_biall_gnomad
    tuple file(fasta), file(fai) from fasta_exome_biall

    output:
    tuple file("af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf.gz"), file("af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf.gz.tbi") into exome_biallelicgz
    file("exome.biall.bed") into exome_pcgr_toml

    script:
    hg = params.assembly == "GRCh37" ? "hg19" : "hg38"
    if( params.assembly == "GRCh37" )
      """
      cut -f 1,2,3 ${exome_bed} > exome.biall.bed
      bgzip ${gnomad}
      tabix ${gnomad}.gz
      gunzip -c ${gnomad}.gz |
      bcftools view -R exome.biall.bed ${gnomad}.gz | bcftools sort -T '.' > af-only-gnomad.exomerh.${hg}.noChr.vcf
      perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.exomerh.hg19.noChr.vcf ${fai} > af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf
      bgzip af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf
      tabix af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf.gz
      """
    else
      """
      cut -f 1,2,3 ${exome_bed} > exome.biall.bed
      gunzip -c ${gnomad} | sed 's/chr//' | bgzip > af-only-gnomad.${hg}.noChr.vcf.gz
      tabix af-only-gnomad.${hg}.noChr.vcf.gz
      bcftools view -R exome.biall.bed af-only-gnomad.${hg}.noChr.vcf.gz | bcftools sort -T '.' > af-only-gnomad.exomerh.${hg}.noChr.vcf
      perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.exomerh.${hg}.noChr.vcf ${fai} > af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf
      bgzip af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf
      tabix af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf.gz
      """
  }

  //exome-specific toml
  process pcgrtoml_exome {

    label 'low_mem'
    publishDir "${params.outDir}/pcgr/data/${params.assemblylc}", mode: "copy"
    errorStrategy 'retry'
    maxRetries 3

    input:
    file(exomebed) from exome_pcgr_toml

    output:
    file("pcgr_configuration_${params.exomeTag}.toml") into exome_tomlo

    when:
    params.pcgrdata

    script:
    """
    sh ${workflow.projectDir}/bin/pcgr_edit_toml.sh \
      ${params.outDir}/pcgr/data/${params.assemblylc}/pcgr_configuration_default.toml \
      ${exomebed} ${exometag}
    """
  }
}

/*
================================================================================
                          4. WGS PROCESS AND PARSING
================================================================================
*/

if(!file("$params.outDir/wgs").exists()){

  process wgs_bed {

    label 'low_mem'
    publishDir path: "${params.outDir}/wgs", mode: "copy"

    input:
    tuple file(fa), file(fai), file(dict) from fasta_dict_wgs

    output:
    file('wgs.bed.interval_list') into complete_wgs
    file('wgs.bed') into (wgs_tabix, wgs_fasta_biallgz)
    tuple file('wgs.bed.gz'), file('wgs.bed.gz.tbi') into gztbi_wgs

    script:
    """
    ##WGS intervals = 1-LN for each chr
    grep @SQ ${dict} | cut -f 2,3 | perl -ane '\$chr=\$F[0];\$chr=~s/SN://;\$end=\$F[1];\$end=~s/LN://;print "\$chr\\t0\\t\$end\\n";' > tmp.wgs.dict.bed

    ##always make interval list so we are in line with fasta
    picard BedToIntervalList I=tmp.wgs.dict.bed O=wgs.bed.interval_list SD=${dict}

    ##output BED
    grep -v "@" wgs.bed.interval_list | cut -f 1,2,3,5 > wgs.bed

    ##tabix
    bgzip -c wgs.bed > wgs.bed.gz
    tabix wgs.bed.gz
    """
  }

  process wgs_biall {

    label 'low_mem'
    publishDir path: "${params.outDir}/wgs", mode: "copy"
    errorStrategy 'retry'
    maxRetries 3

    input:
    file(wgsbed) from wgs_fasta_biallgz
    file(gnomadgz) from wgs_biall_gnomad
    tuple file(fasta), file(fai) from fasta_wgs_biall

    output:
    tuple file("af-only-gnomad.wgs.${hg}.noChr.vcf.gz"), file("af-only-gnomad.wgs.${hg}.noChr.vcf.gz.tbi") into wgs_biallelicgz
    file('wgs.biall.bed') into pcgrtoml_wgs

    script:
    hg = params.assembly == "GRCh37" ? "hg19" : "hg38"
    """
    cut -f 1,2,3 ${wgsbed} > wgs.biall.bed
    gunzip -c ${gnomadgz} | sed 's/chr//' | bgzip > af-only-gnomad.${hg}.noChr.vcf.gz
    tabix af-only-gnomad.${hg}.noChr.vcf.gz

    bcftools view -R wgs.biall.bed af-only-gnomad.${hg}.noChr.vcf.gz | bcftools sort -T '.' > af-only-gnomad.wgsh.${hg}.noChr.vcf
    perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.wgsh.${hg}.noChr.vcf ${fai} > af-only-gnomad.wgs.${hg}.noChr.vcf
    bgzip af-only-gnomad.wgs.${hg}.noChr.vcf
    tabix af-only-gnomad.wgs.${hg}.noChr.vcf.gz
    """
  }
}

if(file("$params.outDir/wgs").exists()){
  if(!file("$params.outDir/pcgr").exists()){
    Channel.fromPath("$params.outDir/wgs").set { wgs_chan }
    process send_pcgrtoml_wgs {

      label 'low_mem'

      input:
      file(wgs_dir) from wgs_chan

      output:
      file('wgs.biall.bed') into pcgrtoml_wgs

      script:
      """
      cut -f 1,2,3 wgs/wgs.bed > wgs.biall.bed
      """
    }
  }
}

/*
================================================================================
                          5. PCGR VEP
================================================================================
*/
if(!file("$params.outDir/pcgr").exists()){

  if(params.localPCGRdata){
    Channel.fromPath("${params.localPCGRdata}").set{ pcgr_data }
  }

  if(!params.localPCGRdata){
    process pcgr_data {

      label 'low_mem'
      errorStrategy 'retry'
      maxRetries 3

      output:
      file("*tgz") into pcgr_data

      script:
      downloadURL = params.assembly == "GRCh37" ? "${params.pcgrURL37}" : "${params.pcgrURL38}"
      """
      wget ${downloadURL}
      """
    }
  }

  //use data
  process pcgr_vep {

    label 'low_mem'
    publishDir "${params.outDir}/pcgr", mode: "copy", pattern: "data"
    errorStrategy 'retry'
    maxRetries 3

    input:
    file(exomebed) from pcgrtoml_exome
    file(wgsbed) from pcgrtoml_wgs
    file(tgz) from pcgr_data

    output:
    file('data') into pcgrdata

    script:
    exometag = params.exomeTag == null ? "" : "${params.exomeTag}"
    """
    tar -xf ${tgz} -C ./

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
}

/*
================================================================================
                          6. HARTWIG MEDICAL (GRIDSS)
================================================================================
*/

if(!file("$params.outDir/gridss").exists()){

  process hartwigmed {

    label 'low_mem'
    publishDir path: "${params.outDir}/gridss", mode: "copy"

    input:
    file(fai) from fai_gridss

    output:
    file('dbs') into gridss_db
    file('refgenomes/human_virus') into gridss_hv
    file('gridss_blacklist.noChr.bed') into gridss_bl
    file('dbs/gridss/gridss.properties') into gridss_pr

    script:
    if( params.assembly == "GRCh37" )
      """
      wget --content-disposition https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download?path=%2FHMFTools-Resources%2FGRIDSS-Purple-Linx-Docker

      7z x GRIDSS-Purple-Linx-Docker.zip
      mv GRIDSS-Purple-Linx-Docker/gpl_ref_data_37.gz gpl_ref_data_hg37.tar.gz
      tar -xf gpl_ref_data_hg37.tar.gz
      rm -rf GRIDSS-Purple-Linx-Docker.zip GRIDSS-Purple-Linx-Docker gpl_ref_data_hg37.tar.gz

      ##blacklist
      sed 's/chr//g' dbs/gridss/ENCFF001TDO.bed > gridss_blacklist.noChr.bed

      perl ${workflow.projectDir}/bin/exact_match_by_col.pl ${fai},0 gridss_blacklist.noChr.bed,0 > gridss_blacklist.noChr.1.bed
      mv gridss_blacklist.noChr.1.bed gridss_blacklist.noChr.bed
      """
    else
      """
      wget --content-disposition https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download?path=%2FHMFTools-Resources%2FGRIDSS-Purple-Linx-Docker

      7z x GRIDSS-Purple-Linx-Docker.zip
      mv GRIDSS-Purple-Linx-Docker/gpl_ref_data_38.gz gpl_ref_data_hg38.tar.gz
      tar -xf gpl_ref_data_hg38.tar.gz
      rm -rf GRIDSS-Purple-Linx-Docker.zip GRIDSS-Purple-Linx-Docker gpl_ref_data_hg38.tar.gz

      ##blacklist
      sed 's/chr//g' dbs/gridss/ENCFF001TDO.bed > gridss_blacklist.noChr.bed
      perl ${workflow.projectDir}/bin/exact_match_by_col.pl ${fai},0 gridss_blacklist.noChr.bed,0 > gridss_blacklist.noChr.1.bed
      mv gridss_blacklist.noChr.1.bed gridss_blacklist.noChr.bed
      """
  }
}

/*
================================================================================
                          7. COSMIC
================================================================================
*/

if(!file("$params.outDir/cosmic/cancer_gene_census.bed").exists()){

  process cosmic_cgc {

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
}

if(!file("$params.outDir/cosmic/CosmicFusionExport.tsv").exists()){

  process cosmic_fusion {

    publishDir "${params.outDir}/cosmic", mode: 'copy'

    when:
    params.cosmicUser && params.cosmicPass

    output:
    file 'CosmicFusionExport.*'

    script:
    """
    ##README to show date of download
    echo \$(date) > CosmicFusionExport-cosmic_dl.README.txt

    ##https://cancer.sanger.ac.uk/cosmic/help/file_download
    curl -H "Authorization: Basic \$(echo ${params.cosmicUser}:${params.cosmicPass} | base64)" https://cancer.sanger.ac.uk/cosmic/file_download/${params.assembly}/cosmic/${params.cosmicVers}/CosmicFusionExport.tsv.gz > url.txt

    curl -o CosmicFusionExport.tsv.gz \$(cut -d '"' -f4 url.txt)
    gunzip CosmicFusionExport.tsv.gz
    """
  }
}

if(!file("$params.outDir/cosmic/CosmicStructExport.tsv").exists()){

  process cosmic_struct {

    publishDir "${params.outDir}/cosmic", mode: 'copy'

    when:
    params.cosmicUser && params.cosmicPass

    output:
    file 'CosmicStructExport.*'

    script:
    """
    ##README to show date of download
    echo \$(date) > CosmicStructExport-cosmic_dl.README.txt

    ##https://cancer.sanger.ac.uk/cosmic/help/file_download
    curl -H "Authorization: Basic \$(echo ${params.cosmicUser}:${params.cosmicPass} | base64)" https://cancer.sanger.ac.uk/cosmic/file_download/${params.assembly}/cosmic/${params.cosmicVers}/CosmicStructExport.tsv.gz > url.txt

    curl -o CosmicStructExport.tsv.gz \$(cut -d '"' -f4 url.txt)
    gunzip CosmicStructExport.tsv.gz
    """
  }
}
/*
================================================================================
                          8. PATHSEQ
================================================================================
*/

if(!file("$params.outDir/pathseq").exists()){

  process Pathseq_HUMDL {

    publishDir "${params.outDir}/pathseq", mode: 'copy'

    when:
    params.microbiome

    output:
    file('*') into psdl_comp

    script:
    def csh = "gs://gcp-public-data--broad-references/hg38/v0/CrossSpeciesContamination/CrossSpeciesHuman"
    """
    gsutil -q cp ${csh}/pathseq_host.bfi ./
    gsutil -q cp ${csh}/pathseq_host.fa.img ./
    """
  }

  process Pathseq_MICDL {

    publishDir "${params.outDir}/pathseq", mode: 'copy'

    when:
    params.microbiome

    output:
    file('*') into psmicdl_comp

    script:
    def csc = "gs://gcp-public-data--broad-references/hg38/v0/CrossSpeciesContamination/CrossSpeciesContaminant"
    """
    gsutil -q cp ${csc}/pathseq_microbe_taxonomy.db ./
    gsutil -q cp ${csc}/pathseq_microbe.dict ./
    gsutil -q cp ${csc}/pathseq_microbe.fa.img ./
    """
  }
}



/*
================================================================================
                          Inf. LEGACY PROCESSES
================================================================================
*/

if(params.legacy){
  // Sequenza GC bias
  process seqnza {

    label 'low_mem'
    publishDir path: "${params.outDir}", mode: "copy"

    input:
    set file(fa), file(fai) from fasta_seqza

    output:
    file('*') into sequenzaout

    script:
    """
    GENOMEGC50GZ=\$(echo ${fa} | sed -r 's/.fasta/.gc50Base.txt.gz/')
    sequenza−utils.py GC-windows −w 50 ${fa} | gzip > \$GENOMEGC50GZ
    """
  }

  // MSIsensor microsatellites
  process msisen {

    label 'low_mem'
    publishDir "${params.outDir}", mode: "copy"

    input:
    set file(fa), file(fai) from fasta_msi

    output:
    file('*') into completedmsisensor

    script:
    """
    msisensor scan -d $fa -o msisensor_microsatellites.list
    """
  }

  // GenomeSize.xml for Pisces
  process gensizxml {

    label 'low_mem'
    publishDir "${params.outDir}", mode: "copy"

    input:
    set file(fa), file(fai), file(dict) from fasta_dict_gensiz

    output:
    file('*') into complete_gensiz

    script:
    """
    echo "<sequenceSizes genomeName=\"${dict}\">" > GenomeSize.xml
    grep "@SQ" ${dict} | while read LINE; do
      CONTIGNAME=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[1]);print \$s[1];' | sed 's/chr//')
      TOTALBASES=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[2]);print \$s[1];')
      MD5SUM=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[3]);print \$s[1];')
      echo -e "\\t<chromosome fileName=\"${fa}\" contigName=\"\$CONTIGNAME\" totalBases=\"\$TOTALBASES\" isCircular=\"false\" md5=\"\$MD5SUM\" ploidy=\"2\" knownBases=\"\$TOTALBASES\" type=\"Chromosome\" />" >> GenomeSize.xml
    done
    echo "</sequenceSizes>" >> GenomeSize.xml
    """
  }

  // Dict processing
  process dict_pr2 {

    label 'low_mem'
    publishDir path: "${params.outDir}", mode: "copy"

    input:
    file(win_dict) from dict_win

    output:
    file('*') into complete_dict

    script:
    """
    perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' ${win_dict} > seq.dict.chr-size

    bedtools makewindows -g seq.dict.chr-size -w 35000000 | perl -ane 'if(\$F[1]==0){\$F[1]++;};print "\$F[0]:\$F[1]-\$F[2]\n";' > 35MB-window.bed
    """
  }

  // KG ASCAT loci
  process ascat_loci {

    label 'low_mem'
    publishDir path: "${params.outDir}", mode: "copy"

    input:
    file(vcf) from ascatloci

    output:
    file('*loci') into complete_ascat

    script:
    """
    LOCIFILE=\$(echo ${vcf} | sed 's/vcf/maf0.3.loci/')
    cat ${vcf} | \
    perl -ane '@s=split(/[=\\;]/,\$F[7]);if(\$s[3]>0.3){print "\$F[0]\\t\$F[1]\\n";}' > \$LOCIFILE
    """
  }
}
