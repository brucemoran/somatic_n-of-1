#! /usr/bin/env bash
SAMPLEID=$4
GERMLINE=$5
export SAMPLEID GERMLINE

gunzip -c $1 | \
perl -ane 'if($F[0]=~m/^#/){if($_=~m/^#CHROM/){
    $_=~s/NORMAL/$ENV{'GERMLINE'}/;
    $_=~s/TUMOR/$ENV{'SAMPLEID'}/;
    print $_;next;
  }
    else{print $_;next;}
  }
else{print $_;}' > $2

perl $3 \
 ID=$SAMPLEID \
 DP=14 \
 MD=2 \
 VCF=$2
