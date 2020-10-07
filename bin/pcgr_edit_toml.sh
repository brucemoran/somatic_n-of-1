##calculate exome size in MB
toml=$1
exomebed=$2
exomeTag=$3

##size of exome
if [[ ! $exomebed == "" ]]; then
  bedtools merge -i $exomebed > exome.biall.merge.bed
  EMB=$(echo -n $(( $(awk '{s+=$3-$2}END{print s}' exome.biall.merge.bed) / 1000000 )))
  export EMB;

  ##perl to parse standard toml config and output ours
  perl -ane 'if($F[0]=~m/^\[mutational_burden/) {
      print "[mutational_burden]\ntmb_intermediate_limit = 10\ntarget_size_mb = $ENV{'EMB'}\n";
    }
    if($F[0]=~m/^vcf2maf/) {
      print "vcf2maf = false\n";
    }
    else { print $_; }' $toml > pcgr_configuration_${exomeTag}.toml
fi

perl -ane 'if($F[0]=~m/^\[mutational_burden/) {
    print "[mutational_burden]\ntmb_intermediate_limit = 10\n;"
  }
  if($F[0]=~m/^vcf2maf/) {
    print "vcf2maf = false\n";
  }
  else { print $_; }' $toml > pcgr_configuration_wgs.toml
