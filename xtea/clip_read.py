##11/04/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu


import pysam
from x_annotation import *
from x_alignments import *
from x_intermediate_sites import *
import global_values
from x_polyA import *
from bwa_align import *

####Trim step 4 also removed class OneClipRead from here (22 lines): a left/right
####clip-position holder that was never instantiated -- grep for it across xtea/ and
####tests/ found only its own definition. Not part of the scan-only change as such;
####found while confirming what ClipReadInfo's deletion would orphan.
####
##Function: pool.map doesnt' accept function in the class
##So, use this wrapper to solve this
def unwrap_self_collect_clip_parts_l(arg, **kwarg):
    return LContigClipReadInfo.collect_clipped_parts_lr_by_chrm(*arg, **kwarg)

####

####Trim step 4 deleted class ClipReadInfo from here (1,062 lines): master's pysam
####clip-read collector. Its only caller was x_TEI_locator.TELocator, deleted in the
####same step -- the C++ xtea-scan pass does this work now. LContigClipReadInfo below
####stays: it is a different collector (clipped parts of a phased-assembly contig
####alignment), reached from x_contig.py via x_local_assembly.py, and untouched by the
####scan-only change.
####
class LContigClipReadInfo():
    def __init__(self, sf_bam, n_jobs, sf_ref, s_working_folder):
        self.sf_bam = sf_bam
        self.working_folder = s_working_folder
        self.n_jobs = n_jobs
        self.sf_reference = sf_ref

    # collect the clipped parts for contigs or long reads of given chromosome
    def collect_clipped_parts_lr_by_chrm(self, record):
        chrm = record[0]  #this is the chrm in the bam file
        sf_bam = record[1]
        sf_sites = record[2]
        islack=record[3]
        working_folder = record[4]

        m_sites = {}
        bam_info = BamInfo(sf_bam, self.sf_reference)
        b_with_chr = bam_info.is_chrm_contain_chr()
        with open(sf_sites) as fin_sites:  # first load in all the sites
            for line in fin_sites:
                fields = line.split()
                site_chrm1 = fields[0]
                site_pos = int(fields[1])
                site_chrm = bam_info.process_chrm_name(site_chrm1, b_with_chr)
                if site_chrm == chrm:
                    m_sites[site_pos]=1

        bamfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        sf_out = working_folder + chrm + global_values.LCLIP_FA_SUFFIX
        with open(sf_out, "w") as fout_clip:
            for algnmt in bamfile.fetch(chrm):
                if algnmt.is_unmapped == True:  # unmapped
                    continue
                l_cigar = algnmt.cigar
                if len(l_cigar) < 1:  # wrong alignment
                    continue
                if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                    continue
                if algnmt.is_supplementary or algnmt.is_secondary:###secondary and supplementary are not considered
                    continue
                query_name = algnmt.query_name
                query_seq = algnmt.query_sequence
                map_pos = algnmt.reference_start
                l_clip_pos=-1
                l_clip_lenth=0
                r_clip_pos=-1
                r_clip_lenth=0
                if l_cigar[0][0] == 4:  # left clipped
                    l_clip_pos=map_pos
                    ###check whether within the range of called clip position
                    hit_site=-1
                    for i in range(-1*islack, islack):
                        tmp_pos=l_clip_pos+i
                        if tmp_pos in m_sites:
                            hit_site=tmp_pos
                            break
                    if hit_site<0:
                        continue

                    l_clip_lenth=l_cigar[0][1]
                    clip_seq=query_seq[:l_clip_lenth]
                    shead=chrm+global_values.SEPERATOR+str(hit_site)+global_values.SEPERATOR+query_name+"L"
                    fout_clip.write(">"+shead+"\n")
                    fout_clip.write(clip_seq+"\n")
                if l_cigar[-1][0] == 4:  # right clipped
                    ##calculate the exact clip position
                    for (type, lenth) in l_cigar[:-1]:
                        if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                            continue
                        else:
                            map_pos += lenth
                    r_clip_pos=map_pos
                    hit_site = -1
                    for i in range(-1 * islack, islack):
                        tmp_pos = r_clip_pos + i
                        if tmp_pos in m_sites:
                            hit_site = tmp_pos
                            break
                    if hit_site < 0:
                        continue

                    r_clip_lenth = l_cigar[-1][1]
                    clip_seq = query_seq[-1 * r_clip_lenth:]
                    shead = chrm + global_values.SEPERATOR + str(hit_site) + global_values.SEPERATOR + query_name + "R"
                    fout_clip.write(">" + shead + "\n")
                    fout_clip.write(clip_seq + "\n")
        bamfile.close()
####
####
    # collect the clipped parts for contigs or long reads of given site
    def collect_clipped_parts_lr_by_region(self, record):
        sf_bam = record[0]
        chrm = record[1]# this is the chrm in the bam file
        ins_pos=record[2]
        iextend=record[3]#will extract reads in region chrm:(ins_pos+/- iextend)
        islack = record[4] #slack around the breakpoint
        working_folder = record[5]

        left_boundary=ins_pos-iextend
        if left_boundary<0:
            left_boundary=0
        right_boundary=ins_pos+iextend

        bamfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        s_site="{0}_{1}".format(chrm, ins_pos)
        sf_out = working_folder + s_site + global_values.LCLIP_FA_SUFFIX
        with open(sf_out, "w") as fout_clip:
            for algnmt in bamfile.fetch(chrm, left_boundary, right_boundary):#
                l_cigar = algnmt.cigar
                map_pos = algnmt.reference_start  # the mapping position
                query_seq = algnmt.query_sequence
                query_name=algnmt.query_name
                if query_seq == None: ###why this happen??????
                    continue
                seq_lenth = len(query_seq)
                lclip_len = 0
                if l_cigar[0][0] == 4:# left-clip
                    brkpnt = ins_pos - islack
                    lclip_len = l_cigar[0][1]
                    clip_pos = map_pos
                    if lclip_len > 100:
                        if brkpnt < clip_pos and brkpnt > (clip_pos - lclip_len):
                            #print "left-clip", algnmt.query_name
                            clip_seq = query_seq[:lclip_len]
                            shead = chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR + query_name + "L"
                            fout_clip.write(">" + shead + "\n")
                            fout_clip.write(clip_seq + "\n")

                if l_cigar[-1][0] == 4:# right clipped
                    brkpnt = ins_pos + islack
                    rclip_len = l_cigar[-1][1]
                    if rclip_len < 100:
                        continue
                    clip_pos = map_pos + seq_lenth - lclip_len - rclip_len
                    if brkpnt > clip_pos and brkpnt < (clip_pos + rclip_len):
                        #print "right-clip", algnmt.query_name
                        clip_seq = query_seq[-1 * rclip_len:]
                        shead = chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR + query_name + "R"
                        fout_clip.write(">" + shead + "\n")
                        fout_clip.write(clip_seq + "\n")

                #####also to consider the case that totally contained in the long read
                #####whole inserted ones
                #####
        bamfile.close()

    ####This function:
    ####1. given dictionary of clip position, ##in format {chrm: {map_pos: (left_cnt, right_cnt)}}
    ####2. write the clipped part into files
    ####iextnd: is the range [-iextnd, iextnd] where to collect the reads.
    ####islack: if the clip position is <= islack distance from insertion-position, then consider it clip at the pos
    def collect_clipped_parts(self, sf_sites, islack, sf_all_clip_fa):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references
        l_chrm_records = []
        xchrom = XChromosome()
        for chrm in references:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###Here filter out those aligned to decoy sequences
                print(("Skip chromosome {0}".format(chrm)))  ##
                continue
            l_chrm_records.append((chrm, self.sf_bam, sf_sites, islack, self.working_folder))
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_clip_parts_l, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()

        ##merge the clipped reads
        with open(sf_all_clip_fa, "w") as fout_all:
            for chrm in references:
                sf_clip_fa = self.working_folder + chrm + global_values.LCLIP_FA_SUFFIX
                if os.path.isfile(sf_clip_fa) == False:
                    continue
                with open(sf_clip_fa) as fin_clip:
                    for line in fin_clip:
                        fout_all.write(line)

####