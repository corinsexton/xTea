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
#   B. pre-fix SVA-only vs post-fix SVA-only  the MECHANISM must change
#        Pre-fix, --types SVA built only the SVA BED, so mate_in_rep counted mate-in-SVA
#        instead of mate-in-(Alu|SVA). What MUST change is that value (realign-out col 6),
#        and it must change monotonically upward, since union >= subset.
#        Whether CALLS move is a separate question and is sample-dependent: measured on
#        HG002, mate_in_rep changed at 1,434,938 positions -- including 38 of the 85 called
#        SVA sites -- and not one call moved, because the SVA cns consistency check was not
#        the binding constraint there. That is a PASS (output-neutral), not a failure.
#        The genuine failure is the mechanism NOT changing, which would mean the pre/post
#        distinction never took effect (wrong worktree, wrong flag) and would make A's
#        "identical" result meaningless.
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

# Per-type VCFs are named from the BAM BASENAME, not the sample name -- a 3-type HG002 run
# writes merged/HG002.novaseq.pcr-free.30x.dedup.grch38_LINE1.vcf, NOT merged/HG002.wgs_LINE1.vcf.
# Exact paths silently match nothing here, which in a gate reads as "no difference found".
# Resolve by glob; the top-level merged VCF IS sample-named.
fam_vcf() {  # $1 = run dir, $2 = family tag (LINE1|ALU|SVA)
    local m
    for m in "$1"/merged/*_"$2".vcf; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done
    return 0
}

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
        "$(count "$(fam_vcf "$d" LINE1)")" \
        "$(count "$(fam_vcf "$d" ALU)")" \
        "$(count "$(fam_vcf "$d" SVA)")" \
        "$(count "$d/${SAMPLE}.vcf")"
done

hr
echo "A. pre-fix 3-type  vs  post-fix 3-type      EXPECT: identical"
if exists "$pre3" && exists "$post3"; then
    a_fail=0; a_ran=0
    a_list=("${SAMPLE}.vcf" "L1/$PT" "Alu/$PT" "SVA/$PT")
    for fam in LINE1 ALU SVA; do
        v=$(fam_vcf "$post3" "$fam"); [ -n "$v" ] && a_list+=("merged/$(basename "$v")")
    done
    for rel in "${a_list[@]}"; do
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
echo "B. pre-fix SVA-only  vs  post-fix SVA-only  EXPECT: the MECHANISM to change"
if exists "$preS" && exists "$postS"; then
    b_moved=0; b_ran=0
    b_list=("SVA/$PT")
    bv=$(fam_vcf "$postS" SVA); [ -n "$bv" ] && b_list+=("merged/$(basename "$bv")")
    for rel in "${b_list[@]}"; do
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
    # Output movement is the plan's PREDICTED signature, not a law -- a correct fix can be
    # output-neutral on a given sample if no call sat on the margin. So test the MECHANISM
    # too, which is what must have changed: mate_in_rep is xtea-scan's 6th realign-out
    # column, and pre-fix it was mate-in-SVA while post-fix it is mate-in-(Alu|SVA).
    # Verdicts:
    #   mechanism moved + output moved    -> PASS, the fix bites here
    #   mechanism moved + output static   -> PASS (output-neutral on this sample)
    #   mechanism did NOT move            -> FAIL: wrong worktree, wrong flag, or inert fix
    echo "  -- mechanism: did mate_in_rep (realign-out col 6) itself change? --"
    b_mech=-1
    pre_t="$preS/scan_clip_realign.tsv"; post_t="$postS/scan_clip_realign.tsv"
    if [ -f "$pre_t" ] && [ -f "$post_t" ]; then
        changed=$(paste <(cut -f1,2,6 "$pre_t") <(cut -f6 "$post_t") | awk -F'\t' '$3!=$4' | wc -l)
        viol=$(paste <(cut -f1,2,6 "$pre_t") <(cut -f6 "$post_t") | awk -F'\t' '$4<$3' | wc -l)
        echo "     positions whose value changed: $changed"
        echo "     positions where post < pre:    $viol  (MUST be 0: union >= subset)"
        [ "$changed" -gt 0 ] && b_mech=1 || b_mech=0
        if [ "$viol" -ne 0 ]; then
            echo "     => MONOTONICITY VIOLATED -- the union count cannot be below the subset count"
            rc=1
        fi
    else
        echo "     SKIP (scan_clip_realign.tsv absent -- needs --keep-intermediates)"
    fi

    if [ "$b_ran" -eq 0 ] && [ "$b_mech" -eq -1 ]; then
        echo "  => B NOT RUN (nothing compared)"; rc=1
    elif [ "$b_mech" -eq 0 ]; then
        echo "  => B FAIL -- mate_in_rep is byte-identical between the two worktrees, so the"
        echo "     pre/post distinction did not take effect at all. Check that the pre-fix"
        echo "     worktree is really pre-fix and that each run used its own flag."
        rc=1
    elif [ "$b_moved" -eq 1 ]; then
        echo "  => B PASS (mechanism changed AND calls moved -- the fix bites on this sample)"
    else
        echo "  => B PASS, OUTPUT-NEUTRAL: mate_in_rep changed at many positions but no SVA"
        echo "     call moved. The fix is correct (monotone, union >= subset) and removes a"
        echo "     real type-set dependence; it simply flipped no decision on THIS sample,"
        echo "     because the SVA cns consistency check was not the binding constraint at"
        echo "     the called sites. Do NOT report this as 'the fix changed nothing'."
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
    sv=$(fam_vcf "$postS" SVA); tv=$(fam_vcf "$post3" SVA)
    if [ -n "$sv" ] && [ -n "$tv" ]; then
        s_n=$(count "$sv"); t_n=$(count "$tv")
        lost=$(comm -23 <(grep -v '^#' "$sv" | awk '{print $1"\t"$2}' | sort -u) \
                        <(grep -v '^#' "$tv" | awk '{print $1"\t"$2}' | sort -u) | wc -l)
        echo "     SVA-only $s_n vs 3-type $t_n; $lost SVA-only positions absent from the 3-type merged set"
        echo "     (expected direction: 3-type <= SVA-only, since L1/Alu can win a shared locus)"
    fi
else rc=1; fi

hr
echo "D. single-type L1 baseline exists (section 6.2 requirement)"
l1v=$(fam_vcf "$post1" LINE1)
if [ -n "$l1v" ]; then
    l1_solo=$(count "$l1v")
    l1_in3=$(count "$(fam_vcf "$post3" LINE1)")
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
