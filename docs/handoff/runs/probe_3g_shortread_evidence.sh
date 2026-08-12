#!/bin/bash
# Same measurement, but each truth locus gets three matched controls on the same
# chromosome (+-1Mb, +2Mb). Without a background rate the clip/disc counts mean nothing.
export PATH=/home/cos689/software/bin:$PATH
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
B=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT005/SMHT005-3G/illuminaNovaseq_bulkWgs/seq_data
U=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_UW/SMHT005/SMHT005-3G/illuminaNovaseq_bulkWgs/seq_data
CRAMS="$B/SMAFI5HKGSER.cram $B/SMAFI6WNA1XM.cram $U/SMAFI4852JH2.cram $U/SMAFI4D3DFD4.cram $U/SMAFIH9MAGML.cram"

probe () {  # chrm pos label
  local chrm=$1 pos=$2 lab=$3
  local lo=$((pos-300)) hi=$((pos+300))
  for c in $CRAMS; do samtools view -T $REF $c $chrm:$lo-$hi 2>/dev/null; done | awk -v P=$pos -v L="$lab" -v C="$chrm:$pos" '
    function boundary(cig, pos_,   c2,op,ln,ref) { ref=0; c2=cig
      while (match(c2,/^[0-9]+[MIDNSHP=X]/)) { op=substr(c2,RLENGTH,1); ln=substr(c2,1,RLENGTH-1)+0
        if (op=="M"||op=="D"||op=="N"||op=="="||op=="X") ref+=ln; c2=substr(c2,RLENGTH+1) }
      return pos_+ref }
    { tot++; cig=$6; seq=$10
      if (($7!="=" && $7!="*") || $9>1000 || $9<-1000) disc++
      if (match(cig,/^[0-9]+S/)) { l=substr(cig,1,RLENGTH-1)+0; d=$4-P
        if (l>=10) { if (d<=15&&d>=-15) c15++; if (d<=50&&d>=-50) c50++ } }
      if (match(cig,/[0-9]+S$/)) { l=substr(cig,RSTART,RLENGTH-1)+0; d=boundary(cig,$4)-P
        if (l>=10) { if (d<=15&&d>=-15) c15++; if (d<=50&&d>=-50) c50++ } } }
    END { printf "%-10s %-20s %6d %6d %6d %6d\n", L, C, tot, c15+0, c50+0, disc+0 }'
}

V=/n/data1/hms/dbmi/park-smaht_dac/DATA/Globus_DataShare/MEI/p25/filtered_calls/v1.0.0/SMHT005-SMaHT-sMEI-v1.0.0.vcf
printf "%-10s %-20s %6s %6s %6s %6s\n" SET LOCUS DEPTH CLIP15 CLIP50 DISC
grep -v '^#' $V | awk -F'\t' '$8 ~ /TISSUE=3G/ {print $1"\t"$2}' | while read chrm pos; do
  probe $chrm $pos TRUTH
  probe $chrm $((pos+1000000)) control
  probe $chrm $((pos-1000000)) control
  probe $chrm $((pos+2000000)) control
done
