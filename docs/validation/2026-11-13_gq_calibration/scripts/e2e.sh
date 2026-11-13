#!/bin/bash
# Rerun the shared genotype + gVCF steps (same commands as the pipeline log) with the gq-calib
# worktree on a copy of a runs_b15 dataset's merged/ inputs, then compare to the offline result.
set -euo pipefail
WT=/n/scratch/users/c/cos689/xtea2_wt/gq_calib
SRC=/n/scratch/users/c/cos689/xtea2_testdata/gtval/runs_b15
OUT=/n/scratch/users/c/cos689/xtea2_review/gq_calib/e2e
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2
export PYTHONNOUSERSITE=1
for N in "$@"; do
  S=$SRC/$N/$N; D=$OUT/$N; rm -rf $D; mkdir -p $D/merged
  cp $S/bam_list.txt $D/
  cp $S/merged/candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt $S/merged/locus_type_map.tsv $D/merged/
  REF=$(grep -m1 -o -- '--ref [^ ]*' $S/${N}_xtea2.log | cut -d' ' -f2)
  RT=$(grep -m1 -o -- '--rtype [0-9]*' $S/${N}_xtea2.log | cut -d' ' -f2)
  I=$D/merged/candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt
  O=$D/merged/candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt
  python $WT/xtea/x_TEA_main.py --gntp_classify -i $I -o $O
  python $WT/xtea/x_TEA_main.py --gVCF -i $O -o $D/merged/ -b $D/bam_list.txt --ref $REF \
    --type-map $D/merged/locus_type_map.tsv --merged-vcf-out $D/$N.vcf --rtype $RT 2>&1 | grep -v '^\[W::' || true
  echo "$N: $(grep -vc '^#' $D/$N.vcf) records"; grep '^##xtea2_genotype_params' $D/$N.vcf
done
