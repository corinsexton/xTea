#!/bin/bash
#SBATCH -A park -p short -t 01:00:00 -c 1 --mem=8G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/germline_check/logs/%A_%a.out
# Re-run --gntp_classify and --gVCF on a copy of a gtval runs_b15 merged/ dir with the xTea2
# baseline (35d287f, arm base) or vaf-site-v2 (arm v2, with the --basic-info/--type-map args its
# shared tail passes, so the per-site VAF model runs). Array index 1-16: dataset x arm.
set -e
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval/runs_b15
OUT=/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/germline_check
DS=(HG002_HiSeq30x HG002_NovaSeq30x HG002_ds10x HG002_ds15x HG002_ds20x HG00438 HG02257 HG02622)
i=$((SLURM_ARRAY_TASK_ID-1)); n=${DS[$((i/2))]}
if [ $((i%2)) = 0 ]; then arm=base; CODE=$OUT/../base_code; else arm=v2; CODE=/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2; fi
src=$G/$n/$n; dst=$OUT/$arm/$n/$n; rm -rf $dst; mkdir -p $dst
for e in $src/* ; do b=$(basename $e); case $b in merged|$n.vcf) ;; *) ln -s $e $dst/$b;; esac; done
mkdir $dst/merged
for f in candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt locus_type_map.tsv; do cp $src/merged/$f $dst/merged/; done
ref=$(grep -o -- '--ref [^ ]*' $src/${n}_xtea2.log | head -1 | cut -d' ' -f2)
M=$dst/merged/candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene
cd $CODE
VA=""; [ $arm = v2 ] && VA="--basic-info $src/tmp_basic/basic_cov_is_rlth_info.txt --type-map $dst/merged/locus_type_map.tsv"
python xtea/x_TEA_main.py --gntp_classify -i $M.txt -o ${M}_gntp.txt $VA
python xtea/x_TEA_main.py --gVCF -i ${M}_gntp.txt -o $dst/merged/ -b $src/bam_list.txt --ref $ref --type-map $dst/merged/locus_type_map.tsv --merged-vcf-out $dst/$n.vcf --rtype 7
echo done $n $arm
