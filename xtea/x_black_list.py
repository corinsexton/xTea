##07/15/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
##
import os
from interval_index import GenomeIntervalIndex

class XBlackList():
    def __init__(self):
        self.m_regions={}
        self.b_with_chr = True #the bed file has "chr"
        self.m_interval_tree = None  # GenomeIntervalIndex, built by load_index_regions

    # This indicates whether the chromosomes is in format "chr1" or "1"
    def set_with_chr(self, b_with_chr):
        self.b_with_chr = b_with_chr

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        if self.b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif self.b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif self.b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    ####load and index bed files with the static sorted-array index
    def load_index_regions(self, sf_bed):
        if os.path.isfile(sf_bed)==False:
            return
        m_starts={}
        m_ends={}
        with open(sf_bed) as fin_bed:
            for line in fin_bed:
                fields=line.split()
                chrm=fields[0]
                istart=int(fields[1])
                iend=int(fields[2])
                m_starts.setdefault(chrm, []).append(istart)
                m_ends.setdefault(chrm, []).append(iend)
        index=GenomeIntervalIndex()
        for chrm in m_starts:
            index.set_chrom(chrm, m_starts[chrm], m_ends[chrm])
        self.m_interval_tree=index

####
    ####B8c FIXED (2026-07-25): this used to return (b_hit, start_pos), where start_pos was the
    ####first element of a SET iteration -- i.e. an arbitrary one of the overlapping blacklist
    ####regions. The author's own note said "here only check hit or not, doesn't matter how many
    ####records in all", and that was true: all five callers bound the position and never read
    ####it. But a dead binding next to an arbitrary value is a trap -- the moment someone starts
    ####using it, they inherit B8's bug silently.
    ####
    ####So the position is gone rather than merely fixed. Returning only the boolean, which is
    ####order-independent, makes the trap impossible to re-arm. See §4 B8c / P1-EVALUATED.
    def fall_in_region(self, ins_chrm1, ins_pos):
        if self.m_interval_tree is None:
            return False
        chrm = self._process_chrm_name(ins_chrm1)
        return self.m_interval_tree.contains(chrm, ins_pos)
####