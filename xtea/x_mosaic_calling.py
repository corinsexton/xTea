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
        m_cutoff = self.get_AF_cutoff(l_types)

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
    # (upper, lower) allele-fraction bounds. Master writes this as a four-branch if/elif chain
    # over the insertion subtype -- and every branch assigns the SAME pair. The per-type
    # distinction was never implemented, so the dict carries no per-type information at all.
    # Kept as one named constant rather than repeated four times, because the fix below only
    # makes sense once that is visible.
    MOSAIC_AF_BOUNDS = (0.35, 0.02)

    def get_AF_cutoff(self, l_types):
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
        class _DefaultingCutoff(dict):
            def __missing__(self, s_type):
                if not hasattr(self, "_warned"):
                    self._warned = set()
                if s_type not in self._warned:
                    self._warned.add(s_type)
                    sys.stderr.write(
                        "WARNING [mosaic]: insertion subtype '{0}' is not enumerated by "
                        "get_rep_type(); using the same AF bounds {1} that every enumerated "
                        "subtype gets.\n".format(s_type, MosaicCaller.MOSAIC_AF_BOUNDS))
                return MosaicCaller.MOSAIC_AF_BOUNDS

        m_cutoff = _DefaultingCutoff()
        for s_type in l_types:
            m_cutoff[s_type] = self.MOSAIC_AF_BOUNDS
        return m_cutoff

####
    def filter_by_blacklist(self, sf_in, sf_black_list, i_slack, sf_out):
        sf_tmp=None
        somatic_caller=SomaticMEICaller(sf_tmp)
        m_MEIs=somatic_caller.load_sites(sf_in)
        m_pass=somatic_caller.filter_by_given_list(sf_black_list, m_MEIs, i_slack)
        somatic_caller.export_slcted_candidates(m_pass, sf_out)

####