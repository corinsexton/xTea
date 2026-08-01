#!/bin/bash
# Stage 0 gate -- DOWNSTREAM_OVERHAUL_PLAN.md section 4:
#   "three-type calls byte-identical across the rename+fix; --types SVA and
#    --types L1,Alu,SVA agree on the SVA call set"
#
# Run after the five HG002 arms land. Login-node safe: this only diffs text files.
#
#   bash docs/handoff/runs/stage0_gate_check.sh
#
# THREE comparisons, and they are not interchangeable:
#
#   A. pre-fix 3-type  vs post-fix 3-type   MUST be identical
#        The rename is cosmetic and the fix cannot touch a three-type run, because
#        requested_types & {alu,sva} already equalled {Alu,SVA}. If this moves, the fix
#        changed something it had no business changing.
#
#   B. pre-fix SVA-only vs post-fix SVA-only  MUST DIFFER
#        This difference IS the bug. Pre-fix, --types SVA built only the SVA BED, so
#        mate_in_rep counted mate-in-SVA instead of mate-in-(Alu|SVA); union >= subset, so
#        the SVA cns consistency check saw a weaker signal. If these are identical, either
#        the fix is inert or the signal never mattered -- both are findings that must be
#        chased, not waved through.
#
#   C. post-fix SVA-only vs post-fix 3-type   the partition property
#        Compared at TWO levels, deliberately:
#          per-type (pre-merge)  -- the SVA pipeline's own calls. This is the real gate:
#                                   free of cross-family effects, so it must agree exactly.
#          merged VCF            -- what a user gets. May legitimately differ, because
#                                   locus_table reassigns borderline multi-family loci
#                                   (priority L1 -> SVA -> Alu), which only happens when
#                                   more than one family was requested. A difference here
#                                   is NOT automatically a failure; a difference at the
#                                   per-type level is.
set -uo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
SAMPLE=HG002.wgs
# per-type call list, pre-merge (post_filtering is the per-type final call set)
PT=candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt

post3=$SCRATCH/baseline0_hg002_L1-Alu-SVA/$SAMPLE
post1=$SCRATCH/baseline0_hg002_L1/$SAMPLE
postS=$SCRATCH/baseline0_hg002_SVA/$SAMPLE
pre3=$SCRATCH/prefix_hg002_L1-Alu-SVA/$SAMPLE
preS=$SCRATCH/prefix_hg002_SVA/$SAMPLE

rc=0
hr() { printf '\n%s\n' "------------------------------------------------------------"; }

count() { [ -f "$1" ] && grep -vc '^#' "$1" 2>/dev/null || echo "-"; }
exists() { [ -d "$1" ] || { echo "  MISSING run dir: $1"; return 1; }; }

echo "============================================================"
echo "STAGE 0 GATE -- HG002"
echo "============================================================"
for d in "$post3" "$post1" "$postS" "$pre3" "$preS"; do
    printf '  %-46s %s\n' "$(basename "$(dirname "$d")")" \
        "$([ -d "$d" ] && echo present || echo MISSING)"
done

hr
echo "CALL COUNTS (merged VCF, per family)"
printf '  %-26s %8s %8s %8s %8s\n' run LINE1 ALU SVA total
for pair in "post-3type:$post3" "post-L1:$post1" "post-SVA:$postS" \
            "pre-3type:$pre3" "pre-SVA:$preS"; do
    name=${pair%%:*}; d=${pair#*:}
    [ -d "$d" ] || { printf '  %-26s %8s\n' "$name" "MISSING"; continue; }
    printf '  %-26s %8s %8s %8s %8s\n' "$name" \
        "$(count "$d/merged/${SAMPLE}_LINE1.vcf")" \
        "$(count "$d/merged/${SAMPLE}_ALU.vcf")" \
        "$(count "$d/merged/${SAMPLE}_SVA.vcf")" \
        "$(count "$d/${SAMPLE}.vcf")"
done

hr
echo "A. pre-fix 3-type  vs  post-fix 3-type      EXPECT: identical"
if exists "$pre3" && exists "$post3"; then
    a_fail=0; a_ran=0
    for rel in "merged/${SAMPLE}_LINE1.vcf" "merged/${SAMPLE}_ALU.vcf" \
               "merged/${SAMPLE}_SVA.vcf" "${SAMPLE}.vcf" \
               "L1/$PT" "Alu/$PT" "SVA/$PT"; do
        if [ ! -f "$pre3/$rel" ] || [ ! -f "$post3/$rel" ]; then
            printf '  %-46s %s\n' "$rel" "SKIP (absent)"; continue
        fi
        a_ran=$((a_ran+1))
        if diff -q "$pre3/$rel" "$post3/$rel" >/dev/null; then
            printf '  %-46s %s\n' "$rel" "identical"
        else
            n=$(diff "$pre3/$rel" "$post3/$rel" | grep -c '^[<>]')
            printf '  %-46s %s\n' "$rel" "DIFFERS ($n lines)"; a_fail=1
        fi
    done
    # A gate that "passes" having compared nothing is worse than no gate -- absent files
    # must read as NOT RUN, never as agreement.
    if   [ "$a_ran" -eq 0 ]; then echo "  => A NOT RUN (0 files compared)"; rc=1
    elif [ "$a_fail" -eq 0 ]; then echo "  => A PASS ($a_ran files compared)"
    else echo "  => A FAIL"; rc=1; fi
else rc=1; fi

hr
echo "B. pre-fix SVA-only  vs  post-fix SVA-only  EXPECT: DIFFERENT (this is the bug)"
if exists "$preS" && exists "$postS"; then
    b_moved=0; b_ran=0
    for rel in "SVA/$PT" "merged/${SAMPLE}_SVA.vcf"; do
        if [ ! -f "$preS/$rel" ] || [ ! -f "$postS/$rel" ]; then
            printf '  %-46s %s\n' "$rel" "SKIP (absent)"; continue
        fi
        b_ran=$((b_ran+1))
        pre_n=$(count "$preS/$rel"); post_n=$(count "$postS/$rel")
        if diff -q "$preS/$rel" "$postS/$rel" >/dev/null; then
            printf '  %-46s identical  (pre %s / post %s)\n' "$rel" "$pre_n" "$post_n"
        else
            gained=$(comm -13 <(grep -v '^#' "$preS/$rel" | awk '{print $1"\t"$2}' | sort -u) \
                              <(grep -v '^#' "$postS/$rel" | awk '{print $1"\t"$2}' | sort -u) | wc -l)
            lost=$(comm -23 <(grep -v '^#' "$preS/$rel" | awk '{print $1"\t"$2}' | sort -u) \
                            <(grep -v '^#' "$postS/$rel" | awk '{print $1"\t"$2}' | sort -u) | wc -l)
            printf '  %-46s DIFFERS  pre %s -> post %s  (+%s / -%s positions)\n' \
                "$rel" "$pre_n" "$post_n" "$gained" "$lost"
            b_moved=1
        fi
    done
    if [ "$b_ran" -eq 0 ]; then
        echo "  => B NOT RUN (0 files compared)"; rc=1
    elif [ "$b_moved" -eq 1 ]; then
        echo "  => B PASS (the fix bites on a single-type run, as predicted)"
    else
        echo "  => B INCONCLUSIVE -- SVA-only output did not move. Either the fix is inert"
        echo "     on this sample or mate_in_rep never gated anything here. CHASE THIS;"
        echo "     do not record the fix as validated on the strength of A alone."
        rc=1
    fi
else rc=1; fi

hr
echo "C. post-fix SVA-only  vs  post-fix 3-type   the partition property"
if exists "$postS" && exists "$post3"; then
    echo "  -- per-type, pre-merge (THE GATE: must agree exactly) --"
    if [ -f "$postS/SVA/$PT" ] && [ -f "$post3/SVA/$PT" ]; then
        if diff -q "$postS/SVA/$PT" "$post3/SVA/$PT" >/dev/null; then
            echo "     SVA per-type call list: identical  => C PASS"
        else
            only_s=$(comm -23 <(awk '{print $1"\t"$2}' "$postS/SVA/$PT" | sort -u) \
                              <(awk '{print $1"\t"$2}' "$post3/SVA/$PT" | sort -u) | wc -l)
            only_3=$(comm -13 <(awk '{print $1"\t"$2}' "$postS/SVA/$PT" | sort -u) \
                              <(awk '{print $1"\t"$2}' "$post3/SVA/$PT" | sort -u) | wc -l)
            echo "     SVA per-type call list: DIFFERS (SVA-only has $only_s extra, 3-type has $only_3 extra)"
            echo "     => C FAIL -- --types does not partition cleanly"
            rc=1
        fi
    else echo "     => C NOT RUN (per-type list absent)"; rc=1; fi

    echo "  -- merged VCF (informational: family reassignment may legitimately differ) --"
    if [ -f "$postS/merged/${SAMPLE}_SVA.vcf" ] && [ -f "$post3/merged/${SAMPLE}_SVA.vcf" ]; then
        s_n=$(count "$postS/merged/${SAMPLE}_SVA.vcf")
        t_n=$(count "$post3/merged/${SAMPLE}_SVA.vcf")
        lost=$(comm -23 <(grep -v '^#' "$postS/merged/${SAMPLE}_SVA.vcf" | awk '{print $1"\t"$2}' | sort -u) \
                        <(grep -v '^#' "$post3/merged/${SAMPLE}_SVA.vcf" | awk '{print $1"\t"$2}' | sort -u) | wc -l)
        echo "     SVA-only $s_n vs 3-type $t_n; $lost SVA-only positions absent from the 3-type merged set"
        echo "     (expected direction: 3-type <= SVA-only, since L1/Alu can win a shared locus)"
    fi
else rc=1; fi

hr
echo "D. single-type L1 baseline exists (section 6.2 requirement)"
if [ -f "$post1/merged/${SAMPLE}_LINE1.vcf" ]; then
    l1_solo=$(count "$post1/merged/${SAMPLE}_LINE1.vcf")
    l1_in3=$(count "$post3/merged/${SAMPLE}_LINE1.vcf")
    echo "  L1-only $l1_solo  vs  L1-within-3-type $l1_in3"
    echo "  (mate_in_rep feeds only the SVA cns path, so L1 should be insensitive to it;"
    echo "   any difference here is cross-family family assignment, not the fix)"
else echo "  MISSING"; rc=1; fi

hr
echo "NOTE: every frozen VCF header carries ##reference=GRCh38_full_analysis_set_plus_decoy_hla.fa"
echo "regardless of the reference actually used (these runs used hg38_no_alt.fa). Pre-existing"
echo "hardcoding in x_gvcf.py, harmless to the comparisons above, but do not read that line as"
echo "provenance -- BASELINE-0/PROVENANCE.txt is the authority."

hr
[ "$rc" -eq 0 ] && echo "STAGE 0 GATE: PASS" || echo "STAGE 0 GATE: NOT CLEAN -- see above"
exit "$rc"
