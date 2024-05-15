##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

import os
from subprocess import *
import xtea.global_values
import numpy as np

class XIntermediateSites():
    
####
    def parse_sites_with_clip_cutoff_for_chrm(self, m_clip_pos_freq, cutoff_left_clip, cutoff_right_clip,
                                              cutoff_clip_mate_in_rep):
        print(f"CUTOFFS IN x_int_sites.py::: {cutoff_left_clip},{cutoff_right_clip},{cutoff_clip_mate_in_rep}")
        m_candidate_sites = {}
        for pos in m_clip_pos_freq:
            ####here need to check the nearby region
            nearby_left_freq = 0
            nearby_right_freq = 0
            nearby_mate_in_rep = 0
            for i in range(-1 * xtea.global_values.NEARBY_REGION, xtea.global_values.NEARBY_REGION):
                i_tmp_pos = pos + i
                if i_tmp_pos in m_clip_pos_freq:
                    nearby_left_freq += m_clip_pos_freq[i_tmp_pos][0]
                    nearby_right_freq += m_clip_pos_freq[i_tmp_pos][1]
                    nearby_mate_in_rep += (
                    m_clip_pos_freq[i_tmp_pos][2] + m_clip_pos_freq[i_tmp_pos][3] + m_clip_pos_freq[i_tmp_pos][4])

            b_candidate=False
            # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip \
            #         and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
            #     b_candidate=True
            if (nearby_left_freq >= cutoff_left_clip or nearby_right_freq >= cutoff_right_clip) \
                    and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
                b_candidate=True

            if b_candidate==True:
                # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip:
                i_left_cnt = m_clip_pos_freq[pos][0]
                i_right_cnt = m_clip_pos_freq[pos][1]
                i_mate_in_rep_cnt = m_clip_pos_freq[pos][2]
                m_candidate_sites[pos] = (i_left_cnt, i_right_cnt, i_mate_in_rep_cnt)
        return m_candidate_sites
####
    ####
    def parse_sites_with_clip_cutoff_for_chrm_with_polyA(self, m_clip_pos_freq, cutoff_left_clip, cutoff_right_clip,
                                              cutoff_clip_mate_in_rep, cutoff_polyA):
        m_candidate_sites = {}
        for pos in m_clip_pos_freq:
            ####here need to check the nearby region
            nearby_left_freq = 0
            nearby_right_freq = 0
            nearby_mate_in_rep = 0
            i_start=-1 * xtea.global_values.NEARBY_REGION
            i_end=xtea.global_values.NEARBY_REGION
            for i in range(i_start, i_end):
                i_tmp_pos = pos + i
                if i_tmp_pos in m_clip_pos_freq:
                    nearby_left_freq += m_clip_pos_freq[i_tmp_pos][0]
                    nearby_right_freq += m_clip_pos_freq[i_tmp_pos][1]
                    nearby_mate_in_rep += (#this is sum up: mate-in-rep + left-to-consensus + right-to-consensus
                    m_clip_pos_freq[i_tmp_pos][2] + m_clip_pos_freq[i_tmp_pos][3] + m_clip_pos_freq[i_tmp_pos][4])

            b_candidate=False
            # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip \
            #         and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
            #     b_candidate=True
            if (nearby_left_freq >= cutoff_left_clip or nearby_right_freq >= cutoff_right_clip) \
                    and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
                b_candidate=True

            if b_candidate==True:
                # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip:
                i_left_cnt = m_clip_pos_freq[pos][0]
                i_right_cnt = m_clip_pos_freq[pos][1]
                i_mate_in_rep_cnt = m_clip_pos_freq[pos][2]
                l_polyA=m_clip_pos_freq[pos][5]
                r_polyA=m_clip_pos_freq[pos][6]
                if (l_polyA+r_polyA)<cutoff_polyA:
                    continue
                m_candidate_sites[pos] = (i_left_cnt, i_right_cnt, i_mate_in_rep_cnt, l_polyA, r_polyA)
        return m_candidate_sites

    ###output the candidate list in a file
    def output_candidate_sites(self, m_candidate_list, sf_out):
        with open(sf_out, "w") as fout_candidate_sites:
            for chrm in m_candidate_list:
                if self.is_decoy_contig_chrms(chrm):  ####decoy and other contigs are not interested!!!!
                    continue
                for pos in m_candidate_list[chrm]:
                    lth = len(m_candidate_list[chrm][pos])
                    fout_candidate_sites.write(chrm + "\t" + str(pos) + "\t")
                    for i in range(lth):
                        s_feature = str(m_candidate_list[chrm][pos][i])
                        fout_candidate_sites.write(s_feature + "\t")
                    fout_candidate_sites.write("\n")

####
    def is_decoy_contig_chrms(self, chrm):
        fields = chrm.split("_")
        if len(fields) > 1:
            return True
        elif chrm == "hs37d5":
            return True

        # if this is not to call mitochondrial insertion, then filter out chrm related reads
        if xtea.global_values.GLOBAL_MITCHONDRION_SWITCH=="OFF":
            if chrm=="MT" or chrm=="chrMT" or chrm=="chrM":#doesn't consider the mitchrondrial DNA
                #print "[TEST]: global value is off"
                return True

        dot_fields = chrm.split(".")
        if len(dot_fields) > 1:
            return True
        else:
            return False
####

    def load_in_candidate_list(self, sf_candidate_list):
        m_list = {}
        with open(sf_candidate_list) as fin_candidate_sites:
            for line in fin_candidate_sites:
                fields = line.split()
                if len(fields)<3:
                    print((fields, "does not have enough fields"))
                    continue
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_list:
                    m_list[chrm] = {}
                if pos not in m_list[chrm]:
                    m_list[chrm][pos] = []
                for ivalue in fields[2:]:
                    m_list[chrm][pos].append(int(ivalue))
        return m_list

    ####
    def load_in_candidate_list_str_version(self, sf_candidate_list):
        m_list = {}
        with open(sf_candidate_list) as fin_candidate_sites:
            for line in fin_candidate_sites:
                fields = line.split()
                if len(fields)<3:
                    print((fields, "does not have enough fields"))
                    continue
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_list:
                    m_list[chrm] = {}
                if pos not in m_list[chrm]:
                    m_list[chrm][pos] = []
                for s_value in fields[2:]:
                    m_list[chrm][pos].append(s_value)
        return m_list

    ####In this version, we will calculate the standard derivation of the left and right clip cluster
    # In the previous step (call_TEI_candidate_sites), some sites close to each other may be introduced together
    # If there are more than 1 site close to each other, than use the peak site as a representative
    def call_peak_candidate_sites_with_std_derivation(self, m_candidate_sites, peak_window):
        m_peak_candidate_sites = {}
        for chrm in m_candidate_sites:
            l_pos = list(m_candidate_sites[chrm].keys())
            l_pos.sort()  ###sort the candidate sites
            pre_pos = -1
            set_cluster = set()
            for pos in l_pos:
                if pre_pos == -1:
                    pre_pos = pos
                    set_cluster.add(pre_pos)
                    continue

                if pos - pre_pos > peak_window:  # find the peak in the cluster
                    # find the representative position of this cluster
                    # also calc the standard derivation of the left and right cluster
                    max_clip = 0
                    tmp_candidate_pos = 0
                    l_lclip_pos=[]
                    l_rclip_pos=[]
                    for tmp_pos in set_cluster:
                        tmp_left_clip = int(m_candidate_sites[chrm][tmp_pos][0])  # left clip
                        tmp_right_clip = int(m_candidate_sites[chrm][tmp_pos][1])  # right clip
                        tmp_all_clip = tmp_left_clip + tmp_right_clip  # all the clip
                        if max_clip < tmp_all_clip:
                            tmp_candidate_pos = tmp_pos
                            max_clip = tmp_all_clip
                        l_tmp_lclip = [tmp_pos] * tmp_left_clip
                        l_tmp_rclip = [tmp_pos] * tmp_right_clip
                        l_lclip_pos.extend(l_tmp_lclip)
                        l_rclip_pos.extend(l_tmp_rclip)
                    ####
                    set_cluster.clear()

                    ####calculate the left and right cluster standard derivation
                    f_lclip_std=self.calc_std_derivation(l_lclip_pos)
                    f_rclip_std=self.calc_std_derivation(l_rclip_pos)
                    l_rclip_pos.extend(l_lclip_pos)
                    f_clip_std=self.calc_std_derivation(l_rclip_pos)

                    if chrm not in m_peak_candidate_sites:
                        m_peak_candidate_sites[chrm] = {}
                    if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip, f_lclip_std, f_rclip_std, f_clip_std]
                pre_pos = pos
                set_cluster.add(pre_pos)
####
            # push out the last group
            max_clip = 0
            tmp_candidate_pos = 0
            l_lclip_pos = []
            l_rclip_pos = []
            for tmp_pos in set_cluster:
                tmp_left_clip = int(m_candidate_sites[chrm][tmp_pos][0])
                tmp_right_clip = int(m_candidate_sites[chrm][tmp_pos][1])
                tmp_all_clip = tmp_left_clip + tmp_right_clip
                if max_clip < tmp_all_clip:
                    tmp_candidate_pos = tmp_pos
                    max_clip = tmp_all_clip
                l_tmp_lclip = [tmp_pos] * tmp_left_clip
                l_tmp_rclip = [tmp_pos] * tmp_right_clip
                l_lclip_pos.extend(l_tmp_lclip)
                l_rclip_pos.extend(l_tmp_rclip)
            ####calculate the left and right cluster standard derivation
            f_lclip_std = self.calc_std_derivation(l_lclip_pos)
            f_rclip_std = self.calc_std_derivation(l_rclip_pos)
            l_rclip_pos.extend(l_lclip_pos)
            f_clip_std = self.calc_std_derivation(l_rclip_pos)
            if chrm not in m_peak_candidate_sites:
                m_peak_candidate_sites[chrm] = {}
            if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                ##Here, use list in order to output the list (by call the output_candidate_sites function)
                m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip, f_lclip_std, f_rclip_std, f_clip_std]
        return m_peak_candidate_sites

####
    def merge_clip_disc(self, sf_disc_tmp, sf_clip, sf_out):
        with open(sf_out, "w") as fout_list:
            m_disc={}
            with open(sf_disc_tmp) as fin_disc:
                for line in fin_disc:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    s_left_disc=fields[2]
                    s_right_disc=fields[3]
                    if chrm not in m_disc:
                        m_disc[chrm]={}
                    m_disc[chrm][pos]=(s_left_disc, s_right_disc)#
            with open(sf_clip) as fin_clip:
                for line in fin_clip:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    if chrm not in m_disc:
                        print(("Error happen at merge clip and disc feature step: {0} not exist".format(chrm)))
                        continue
                    if pos not in m_disc[chrm]:
                        continue
                    s_left_disc=m_disc[chrm][pos][0]
                    s_right_disc=m_disc[chrm][pos][1]
                    fields.append(s_left_disc)
                    fields.append(s_right_disc)
                    fout_list.write("\t".join(fields) + "\n")
#
    ####This is to output all the candidates with all the clip, discord, barcode information in one single file
    def merge_clip_disc_barcode(self, sf_barcode_tmp, sf_disc, sf_out):
        with open(sf_out, "w") as fout_list:
            m_barcode={}
            with open(sf_barcode_tmp) as fin_barcode:
                for line in fin_barcode:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    s_nbarcode=fields[-1]
                    if chrm not in m_barcode:
                        m_barcode[chrm]={}
                    m_barcode[chrm][pos]=s_nbarcode
            with open(sf_disc) as fin_disc:
                for line in fin_disc:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    if chrm not in m_barcode:
                        print(("Error happen at merge clip, disc and barcode step: {0} not exist".format(chrm)))
                        continue
                    if pos not in m_barcode[chrm]:
                        continue
                    s_barcode=m_barcode[chrm][pos]
                    fields.append(s_barcode)
                    fout_list.write("\t".join(fields) + "\n")

    # In the previous step (call_TEI_candidate_sites), some sites close to each other may be introduced together
    # If there are more than 1 site close to each other, than use the peak site as a representative
    def call_peak_candidate_sites_all_features(self, m_candidate_sites, peak_window):
        m_peak_candidate_sites = {}
        for chrm in m_candidate_sites:
            l_pos = list(m_candidate_sites[chrm].keys())
            l_pos.sort()  ###sort the candidate sites
            pre_pos = -1
            set_cluster = set()
            for pos in l_pos:
                if pre_pos == -1:
                    pre_pos = pos
                    set_cluster.add(pre_pos)
                    continue

                if pos - pre_pos > peak_window:  # find the peak in the cluster
                    max_clip = 0
                    max_all=0
                    tmp_candidate_pos = 0
                    for tmp_pos in set_cluster:
                        tmp_clip = int(m_candidate_sites[chrm][tmp_pos][0])  # # ofclip related features
                        tmp_all = int(m_candidate_sites[chrm][tmp_pos][1])  # # of all features
                        if (max_clip < tmp_clip) or (max_clip==tmp_clip and max_all < tmp_all):
                            tmp_candidate_pos = tmp_pos
                            max_clip = tmp_clip
                            max_all=tmp_all
                    set_cluster.clear()
                    if chrm not in m_peak_candidate_sites:
                        m_peak_candidate_sites[chrm] = {}
                    if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
                pre_pos = pos
                set_cluster.add(pre_pos)
            # push out the last group
            max_clip = 0
            max_all = 0
            tmp_candidate_pos = 0
            for tmp_pos in set_cluster:
                tmp_clip = int(m_candidate_sites[chrm][tmp_pos][0])  # # ofclip related features
                tmp_all = int(m_candidate_sites[chrm][tmp_pos][1])  # # of all features
                if (max_clip < tmp_clip) or (max_clip == tmp_clip and max_all < tmp_all):
                    tmp_candidate_pos = tmp_pos
                    max_clip = tmp_clip
                    max_all = tmp_all
            set_cluster.clear()
            if chrm not in m_peak_candidate_sites:
                m_peak_candidate_sites[chrm] = {}
            if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
        return m_peak_candidate_sites

    ###this function is used to combine some sites are close to each other
    def combine_closing_sites(self, sf_input, iwindow, sf_out):
        m_original_sites={}
        with open(sf_input) as fin_sites:
            for line in fin_sites:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                cur_sum_clip=int(fields[2]) + int(fields[3]) + int(fields[4])
                cur_sum_all = cur_sum_clip + int(fields[5]) + int(fields[6])

                if chrm not in m_original_sites:
                    m_original_sites[chrm]={}
                m_original_sites[chrm][pos]=(cur_sum_clip, cur_sum_all)

        m_peak_candidate_sites=self.call_peak_candidate_sites_all_features(m_original_sites, iwindow)
        with open(sf_input) as fin_sites, open(sf_out, "w") as fout_sites:
            for line in fin_sites:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                if (chrm in m_peak_candidate_sites) and (pos in m_peak_candidate_sites[chrm]):
                    fout_sites.write(line)

    def are_sites_close(self, pos1, pos2, iwin):
        i_start=pos1-iwin
        i_end=pos1+iwin
        for i in range(i_start, i_end):
            if i==pos2:
                return True
        return False

    ####calculate the std derivation of a list of positions
    def calc_std_derivation(self, l_pos):
        if len(l_pos)<=0:
            return -1
        b = np.array([l_pos])
        f_std = round(np.std(b), 2)
        return f_std
####