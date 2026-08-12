#!/bin/bash
# Does the HapMap-mixture truth set carry short-read signal, and at which VAF does it die?
# Stratifies tier1 L1 truth loci by predicted VAF, counts clip/discordant evidence against
# matched same-chromosome controls. This is the check that SMHT005-3G failed (plan sec 7.2).
export PATH=/home/cos689/software/bin:$PATH
BAM=${BAM:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/HapMap_Mixture/illuminaNovaseq_bulkWgs/seq_data/SMAFI1LJ41ED.bam}
D=/n/data1/hms/dbmi/park/corinne/hapmap_mei_calls/hapmap_mix_calls/truth_set_10_18_2024
N=${N:-20}

probe () {
  samtools view "$BAM" $1:$(($2-300))-$(($2+300)) 2>/dev/null | awk -v P=$2 '
    function bnd(cig,p,  c2,op,ln,r){r=0;c2=cig; while(match(c2,/^[0-9]+[MIDNSHP=X]/)){op=substr(c2,RLENGTH,1);ln=substr(c2,1,RLENGTH-1)+0; if(op~/[MDN=X]/)r+=ln; c2=substr(c2,RLENGTH+1)} return p+r}
    { tot++; cig=$6
      if (($7!="=" && $7!="*") || $9>1000 || $9<-1000) disc++
      if (match(cig,/^[0-9]+S/)) { l=substr(cig,1,RLENGTH-1)+0; d=$4-P; if(l>=10&&d<=15&&d>=-15) clip++ }
      if (match(cig,/[0-9]+S$/)) { l=substr(cig,RSTART,RLENGTH-1)+0; d=bnd(cig,$4)-P; if(l>=10&&d<=15&&d>=-15) clip++ } }
    END { printf "%d %d %d\n", tot, clip+0, disc+0 }'
}

printf "%-10s %-5s %8s %8s %8s %8s\n" VAF_BIN n "depth" "clip_T" "clip_C" "disc_T/C"
for V in 0.0025 0.0100 0.0500 0.0700; do
  ct=0; cc=0; dt=0; dc=0; dep=0; k=0
  while IFS=$'\t' read id vaf; do
    printf -v r "%.4f" "$vaf"; [ "$r" = "$V" ] || continue
    k=$((k+1)); [ $k -gt $N ] && break
    chrm=${id%:*}; pos=${id#*:}
    read t c d <<< "$(probe $chrm $pos)";       dep=$((dep+t)); ct=$((ct+c)); dt=$((dt+d))
    read t2 c2 d2 <<< "$(probe $chrm $((pos+1000000)))"; cc=$((cc+c2)); dc=$((dc+d2))
  done < <(tail -n +2 $D/tier1_L1_VAFs.tsv)
  [ $k -eq 0 ] && continue
  awk -v v=$V -v k=$k -v dep=$dep -v ct=$ct -v cc=$cc -v dt=$dt -v dc=$dc \
    'BEGIN{printf "%-10s %-5d %8.0f %8.2f %8.2f %8.2f/%.2f\n", v, k, dep/k, ct/k, cc/k, dt/k, dc/k}'
done
