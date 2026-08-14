####
##04/23/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

#This module is designed for mosaic TE insertion calling from normal bulk tissue.
# Assume the insertions happen at low AF (say <1%).
# The data is sequenced at high coverage (say >200X), and low cutoff is set (say 2 disc and 2 clip)

####1. add blacklist, like gnomAD-SV
####2. filter with callset from other samples
####3. screenshot discordant and clipped reads only
####

# Explicit, not inherited from the star-imports below: get_AF_cutoff's warning writes to
# sys.stderr, and relying on `sys` leaking in through `from x_post_filter import *` means a
# tidy-up in an unrelated module turns a warning into a NameError at the worst moment.
import os
import sys

from x_post_filter import *
from x_somatic_calling import *

####
####this is inherited from XPostFilter class
class MosaicCaller(XPostFilter):
    def __init__(self, swfolder, n_jobs):
        XPostFilter.__init__(self, swfolder, n_jobs)
        self.brkpnt_slack=100

    ####This is assume the input is generated from high coverage data (with low cutoff)
    # def call_mosaic_high_coverage_output(self, sf_in):

    def run_call_mosaic(self, sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, sf_black_list, sf_new_out):
        xtea_parser = XTEARsltParser()
        xtea_parser.set_rep_support_type(self._two_side, self._one_half_side, self._one_side, self._other)
        l_old_rcd = xtea_parser.load_in_xTEA_rslt(sf_xtea_rslt)
        xannotation = self.construct_interval_tree(sf_rmsk, i_min_copy_len, self.b_rslt_with_chr)
        xtprt_filter = XTPRTFilter(self.swfolder, self.n_jobs)
        af_filter = AFConflictFilter(self.swfolder, self.n_jobs)
        l_types = af_filter.get_rep_type()
        m_cutoff = self.get_AF_cutoff(l_types, i_rep_type)

        ####Hard code here !!!!!!!!
        sf_new_out_bf_black_list=sf_new_out+".before_blacklist.txt"
        if (i_rep_type & 1) != 0:
            self.call_mosaic_L1_from_bulk(l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff,
                                          sf_new_out_bf_black_list)
        elif (i_rep_type & 2) != 0:
            self.call_mosaic_Alu_from_bulk(l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff,
                                     sf_new_out_bf_black_list)
        elif (i_rep_type & 4) != 0:
            self.call_mosaic_SVA_from_bulk(l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff,
                                     sf_new_out_bf_black_list)

        ##go through blacklist to filter out germline cases
        self.filter_by_blacklist(sf_new_out_bf_black_list, sf_black_list, self.brkpnt_slack, sf_new_out)
####
####
    def call_mosaic_L1_from_bulk(self, l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff, sf_new_out):
        with open(sf_new_out, "w") as fout_new:
            for rcd in l_old_rcd:
                # WITHDRAWN, deliberately not deleted. A rule dropping the raw subtype
                # one_side_coverage_conflict was added here on 2026-08-13 on a measured 91% FP
                # rate for Alu (30/33) and 100% for L1 (2/2). Both were counted on a call set
                # produced under the multi-BAM genotype-feature merge bug fixed in dca7363: the
                # subtype field itself is NOT corrupted by that bug (verified: zero subtype
                # flips between a buggy and a fixed cns_filter run over broad's SVA), but WHICH
                # calls survived to be counted was decided by AF gates reading corrupted
                # denominators, so the TP/FP split those rates came from is not trustworthy.
                # The rule may well be right -- re-derive it from a run made after dca7363
                # before putting it back.
                s_rep_supt_type = xtea_parser.get_ins_sub_type(rcd)
                ins_chrm = rcd[0]
                ins_pos = rcd[1]
                #boolean-only: i_pos was bound and never read. See §4 P1-EVALUATED.
                b_in_rep = xannotation.is_within_repeat_region_bool(ins_chrm, int(ins_pos))
                b_with_polyA = xtprt_filter.has_polyA_signal(rcd)
                if b_with_polyA is False:
                    continue
                if af_filter.is_qualified_mosaic_rcd(rcd[-1], m_cutoff) == False:
                    continue
                if s_rep_supt_type == self._two_side:  ####two sides
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) == True:
                        continue
                elif s_rep_supt_type == self._one_half_side:  ##one and half side
                    # check whether fall in repetitive region of the same type
                    if b_in_rep == True:
                        continue
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) is True:
                        continue
                elif s_rep_supt_type == self._one_side:###one side
                    if b_in_rep is True:
                        continue
                        # if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                else:  # for other type, just skip
                    continue
                # save the passed ones
                s_in_rep = "not_in_LINE1_copy"
                if b_in_rep is True:
                    s_in_rep = "Fall_in_LINE1_copy"
                s_pass_info = rcd[-1].rstrip() + "\t" + s_in_rep + "\n"
                fout_new.write(s_pass_info)
####
####
    def call_mosaic_Alu_from_bulk(self, l_old_rcd, xtea_parser, xtprt_filter,
                                  af_filter, xannotation, m_cutoff, sf_new_out):
        with open(sf_new_out, "w") as fout_new:
            for rcd in l_old_rcd:
                # one_side_coverage_conflict rule withdrawn here too -- see
                # call_mosaic_L1_from_bulk for why the evidence behind it does not hold.
                s_rep_supt_type = xtea_parser.get_ins_sub_type(rcd)
                ins_chrm = rcd[0]
                ins_pos = rcd[1]
                #boolean-only: i_pos was bound and never read. See §4 P1-EVALUATED.
                b_in_rep = xannotation.is_within_repeat_region_bool(ins_chrm, int(ins_pos))
                b_with_polyA = xtprt_filter.has_polyA_signal(rcd)
                if b_with_polyA is False:
                    continue
                if af_filter.is_qualified_mosaic_rcd(rcd[-1], m_cutoff) == False:
                    continue
                if s_rep_supt_type is self._two_side:  ####two sides
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) == True:
                        continue
                elif s_rep_supt_type is self._one_half_side:  ##one and half side
                    # check whether fall in repetitive region of the same type
                    if b_in_rep is True:
                        continue
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) is True:
                        continue
                elif s_rep_supt_type is self._one_side:  ###one side
                    if b_in_rep is True:
                        continue
                        # if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                else:  # for other type, just skip
                    continue
                # save the passed ones
                s_in_rep = "not_in_Alu_copy"
                if b_in_rep is True:
                    s_in_rep = "Fall_in_Alu_copy"
                s_pass_info = rcd[-1].rstrip() + "\t" + s_in_rep + "\n"
                fout_new.write(s_pass_info)

####
    def call_mosaic_SVA_from_bulk(self, l_old_rcd, xtea_parser, xtprt_filter,
                                  af_filter, xannotation, m_cutoff, sf_new_out):
        with open(sf_new_out, "w") as fout_new:
            for rcd in l_old_rcd:
                # one_side_coverage_conflict rule withdrawn here too -- see
                # call_mosaic_L1_from_bulk for why the evidence behind it does not hold.
                s_rep_supt_type = xtea_parser.get_ins_sub_type(rcd)
                ins_chrm = rcd[0]
                ins_pos = rcd[1]
                #boolean-only: i_pos was bound and never read. See §4 P1-EVALUATED.
                b_in_rep = xannotation.is_within_repeat_region_bool(ins_chrm, int(ins_pos))
                b_with_polyA = xtprt_filter.has_polyA_signal(rcd)
                if b_with_polyA is False:
                    continue
                if af_filter.is_qualified_mosaic_rcd(rcd[-1], m_cutoff) == False:
                    continue
                if s_rep_supt_type is self._two_side:  ####two sides
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) == True:
                        continue
                elif s_rep_supt_type is self._one_half_side:  ##one and half side
                    # check whether fall in repetitive region of the same type
                    if b_in_rep is True:
                        continue
                    # if xtprt_filter.is_polyA_dominant_two_side(rcd) is True:
                elif s_rep_supt_type is self._one_side:  ###one side
                    if b_in_rep is True:
                        continue
                        # if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                else:  # for other type, just skip
                    continue
                # save the passed ones
                s_in_rep = "not_in_SVA_copy"
                if b_in_rep is True:
                    s_in_rep = "Fall_in_SVA_copy"
                s_pass_info = rcd[-1].rstrip() + "\t" + s_in_rep + "\n"
                fout_new.write(s_pass_info)
####
####
    # (clip_upper, clip_lower, disc_upper). Master wrote clip_upper and disc_upper as the SAME
    # value -- is_ins_pass_mosaic_cutoff used one shared f_upper_af as the ceiling for BOTH the
    # clip-based allele fraction and the disc-read concordance ratio. Kept as a 3-tuple so the
    # two can move independently if a future calibration needs it; today they do not.
    #
    # THERE IS NO PER-FAMILY OVERRIDE, AND THE HISTORY IS WORTH KEEPING. An Alu-specific
    # disc_upper of 0.80 was added here on 2026-08-13 on the strength of a measurement that
    # true Alu positives carried f_disc_concod 0.421-0.754, entirely above 0.35, so that all
    # 156 of them were being rejected. That measurement was taken on output produced by the
    # multi-BAM genotype-feature merge bug fixed in dca7363: n_disc_pairs and n_concd_pairs are
    # the ratio's numerator and denominator and BOTH were inflated by one fixed unrelated
    # site's counts. Recomputed from the same run's per-BAM feature files with the merge done
    # correctly, true-positive f_disc_concod is:
    #
    #     L1  median 0.055 (was 0.295 under the bug)   0.9% of truth records exceed 0.35
    #     Alu median 0.052 (was 0.457)                 2.3% exceed 0.35
    #     SVA median 0.061 (was 0.459)                 8.4% exceed 0.35
    #
    # So the shared 0.35 ceiling was never the problem for any family -- it costs a couple of
    # percent of true positives, not all of them -- and 0.80 was pure over-correction for an
    # arithmetic artefact. Raising it again would only admit background. If a per-family
    # ceiling is ever wanted, derive it from a run made AFTER dca7363, not from these files.
    #
    # THE CLIP-AF FLOOR IS 0, LOWERED FROM 0.02 (2026-08-14). Master's 0.02 excluded the two
    # lowest strata of the truth set BY CONSTRUCTION: a 1-1.25% VAF insertion produces a clip
    # AF of about that, so the floor rejected it no matter how much coverage was thrown at it.
    # Measured on bcm (450x), truth-matching records reaching this filter, binned by the truth
    # locus's own VAF:
    #
    #     truth VAF   n     median clip AF   below the old 0.02 floor
    #     1-1.25%     283   0.0150           207/283  (73%)
    #     2-3.25%     150   0.0293            33/150  (22%)
    #     5-5.25%     351   0.0481             6/351   (2%)
    #     6-8.25%     401   0.0769             1/401   (0%)
    #
    # The detection was working; the gate was discarding it. That is why bcm's 1-1.25%
    # sensitivity (5.8%) was LOWER than shallower broad's (7.6%) -- extra depth cannot lift a
    # variant's allele fraction above a floor set higher than the variant.
    #
    # Swept 0.02 -> 0 on two arms, everything else held fixed (each point is a real
    # postF_mosaic re-run scored against the truth set):
    #
    #     arm    family   TP at 0.02 -> at 0     sensitivity        precision
    #     bcm    Alu       891 -> 1101           37.2% -> 46.0%     60.9% -> 62.8%  (RISES)
    #     bcm    L1         69 -> 72             16.4% -> 17.1%     66.3% -> 66.1%
    #     broad  Alu       751 -> 796            31.4% -> 33.2%     69.1% -> 69.2%
    #     broad  SVA        49 -> 53             27.4% -> 29.6%     65.3% -> 65.4%
    #     broad  L1         58 -> 59             13.7% -> 14.0%     84.1% -> 81.9%  (only loss)
    #
    # Precision RISING on the largest truth set (Alu, 2395 loci) is the load-bearing result:
    # what the floor was rejecting was enriched for real insertions, not background. The gain
    # also scales with depth (+213 TP on bcm vs +50 on broad), which noise-admission would not.
    #
    # ZERO RATHER THAN A SMALL POSITIVE VALUE, and the distinction is measured, not assumed:
    # the sweep SATURATES above zero -- at 0.01 on broad, 0.005 on bcm (the saturation point
    # falls as depth rises, since deeper coverage resolves lower allele fractions). Below each
    # arm's saturation point, 0.005/0.002/0.001/0.0 return call sets identical to within one
    # record. So zero is not "as permissive as possible"; it is the same answer every value
    # below saturation gives, chosen so no arm of a future depth ladder sits above its own
    # saturation point. What keeps it safe is that this is not the only low-AF gate:
    # has_polyA_signal and the subtype branches run alongside it, and nothing survives those
    # below ~0.005 in any arm measured. If that ever stops holding -- an input where records
    # with literally zero clip support reach here -- this is the constant to revisit first.
    MOSAIC_AF_BOUNDS = (0.35, 0.0, 0.35)

    def get_AF_cutoff(self, l_types, i_rep_type=1):
        """Subtype -> AF bounds, TOTAL over subtypes rather than only the enumerated ones.

        THE BUG THIS FIXES, found by running it (job 50434768, L1:postF_mosaic, exit 1):

            KeyError: 'orphan_or_sibling_transduction'   x_post_filter.py:1125

        `is_ins_pass_mosaic_cutoff` looks the record's own subtype (field 32) up in this dict,
        but the dict is built from `AFConflictFilter.get_rep_type()`, which enumerates 13
        subtypes and OMITS orphan_or_sibling_transduction -- even though
        global_values.ORPHAN_TRANSDUCTION defines it. On the HapMap chr1 mosaic run that
        subtype was 20 of 32 L1 records, so this is the dominant case, not a rare one; mosaic
        crashed on the first record it reached.

        Returning a defaulting map rather than adding the missing name to get_rep_type() is
        deliberate. get_rep_type() is shared with the GERMLINE af-conflict filter
        (x_post_filter.py:400, :1138); appending to it would change germline behaviour to fix a
        mosaic crash. This confines the change to the mosaic path.

        Defaulting is sound here precisely BECAUSE every branch assigns the same pair: a
        missing key can only mean "a subtype get_rep_type() does not enumerate", and the
        correct value is the one every enumerated subtype already gets. It is not a guess.
        Should the bounds ever become genuinely per-type, this default must be revisited --
        hence the warning, which names each unknown subtype once so the gap stays visible
        instead of being silently absorbed.
        """
        # i_rep_type is still taken (and still threaded through by run_call_mosaic) so a
        # per-family bound can be reinstated without touching call sites -- see the note on
        # MOSAIC_AF_BOUNDS for why the Alu one was withdrawn.
        bounds = self.MOSAIC_AF_BOUNDS

        class _DefaultingCutoff(dict):
            def __missing__(self, s_type):
                if not hasattr(self, "_warned"):
                    self._warned = set()
                if s_type not in self._warned:
                    self._warned.add(s_type)
                    sys.stderr.write(
                        "WARNING [mosaic]: insertion subtype '{0}' is not enumerated by "
                        "get_rep_type(); using the same AF bounds {1} that every enumerated "
                        "subtype gets.\n".format(s_type, bounds))
                return bounds

        m_cutoff = _DefaultingCutoff()
        for s_type in l_types:
            m_cutoff[s_type] = bounds
        return m_cutoff

####
    # "null" is this pipeline's established sentinel for "no blacklist given" -- it is the
    # DEFAULT for both entry points (run_xtea.py:1282, x_TEA_main.py:211) and run_xtea
    # normalises None/null/NULL onto it (:1767). Every other post-filter step tolerates it:
    # postF_hc ran with `--blacklist null` in the same chain and succeeded.
    #
    # This one did not. It passed the sentinel straight to load_sites(), which opened it:
    #
    #     FileNotFoundError: [Errno 2] No such file or directory: 'null'
    #
    # found by running it (job 50435713, broad WGS, L1:postF_mosaic, 4h07m in -- and only
    # reached because the KeyError fix above let it get that far). Every mosaic run that does
    # not explicitly name a blacklist dies on its very last step, after all the expensive work.
    #
    # Honouring the sentinel means passing the records through UNFILTERED, which is what "no
    # blacklist" means. It deliberately does not substitute a default blacklist: silently
    # filtering against a file the caller never asked for would change the call set on the
    # quiet, and for the HapMap benchmark it would also change what precision measures.
    _NO_BLACKLIST = ("null", "NULL", "", None)

    def filter_by_blacklist(self, sf_in, sf_black_list, i_slack, sf_out):
        somatic_caller = SomaticMEICaller(None)
        if sf_black_list in self._NO_BLACKLIST or not os.path.isfile(sf_black_list):
            if sf_black_list not in self._NO_BLACKLIST:
                sys.stderr.write("WARNING [mosaic]: blacklist '{0}' does not exist; passing "
                                 "records through unfiltered\n".format(sf_black_list))
            else:
                sys.stderr.write("[mosaic] no blacklist given ('{0}'); records pass through "
                                 "unfiltered\n".format(sf_black_list))
            m_MEIs = somatic_caller.load_sites(sf_in)
            somatic_caller.export_slcted_candidates(m_MEIs, sf_out)
            return
        m_MEIs = somatic_caller.load_sites(sf_in)
        m_pass = somatic_caller.filter_by_given_list(sf_black_list, m_MEIs, i_slack)
        somatic_caller.export_slcted_candidates(m_pass, sf_out)

####