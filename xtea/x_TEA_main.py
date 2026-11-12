##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com
##

####
import os
from shutil import copyfile
import global_values
import family_params
from x_local_assembly import *
from x_intermediate_sites import *
from x_reference import *
from x_clip_disc_filter import *
from x_somatic_calling import *
from optparse import OptionParser
from x_reads_collection import *
from x_mutation import *
from x_gene_annotation import *
from x_genotype_feature import *
from x_basic_info import *
from x_parameter import *
from x_post_filter import *
from x_mosaic_calling import *
from x_joint_calling import *
from x_igv import *
from x_gvcf import *
import x_genotype_lik
from x_orphan_transduction import *
#
####
##parse the options
def parse_option():
    parser = OptionParser()
    parser.add_option("-P", "--preprocess",
                      action="store_true", dest="preprocess", default=False,
                      help="Preprocessing stpes")
    parser.add_option("-N", "--filter_csn",
                      action="store_true", dest="filter_csn", default=False,
                      help="Filter out candidate sites from map position on consensus")
    parser.add_option("--resume",
                      action="store_true", dest="resume", default=False,
                      help="Resume the running, which will skip the step if output file already exists!")
    #transduction
    parser.add_option("--transduction",
                      action="store_true", dest="transduction", default=False,
                      help="Call transduction for sites")
    parser.add_option("--sibling",
                      action="store_true", dest="sibling", default=False,
                      help="Call sibling transduction for sites")
    parser.add_option("--spectrum",
                      action="store_true", dest="spectrum", default=False,
                      help="Spectrum analysis by tumor type")

    parser.add_option("-E", "--collect",
                      action="store_true", dest="collect", default=False,
                      help="Collect reads for candidate sites")
    parser.add_option("-I", "--mutation",
                      action="store_true", dest="mutation", default=False,
                      help="Call internal mutation")
    parser.add_option("-U", "--collect_Illumina",
                      action="store_true", dest="collect_illumina", default=False,
                      help="Collect reads for candidate sites from normal illumina alignment")
    parser.add_option("-F", "--filter_asm",
                      action="store_true", dest="filter_asm", default=False,
                      help="Filter out candidate sites from assembly")
    parser.add_option("-G", "--contig_realign",
                      action="store_true", dest="contig_realign", default=False,
                      help="Filter out candidate sites from assembly")
    parser.add_option("-T", "--trace",
                      action="store_true", dest="trace", default=False,
                      help="Trace the sources of TEIs")
    parser.add_option("-A", "--assembly",
                      action="store_true", dest="assembly", default=False,
                      help="Do local assembly for collected reads")
    parser.add_option("-L", "--local",
                      action="store_true", dest="local", default=False,
                      help="Assemble the TEIs on local machine")
    parser.add_option("-M", "--map",
                      action="store_true", dest="map", default=False,
                      help="map flank regions to the assembled contigs")
    parser.add_option("-V", "--visualization",
                      action="store_true", dest="visualization", default=False,
                      help="Show the heatmap figure of the selected regions")
    parser.add_option("-K", "--withflank",
                      action="store_true", dest="withflank", default=False,
                      help="Keep the flank regions with the repeat copies")
    parser.add_option("-J", "--joint",
                      action="store_true", dest="joint", default=False,
                      help="Joint calling")
    parser.add_option("--mit",
                      action="store_true", dest="mit", default=False,
                      help="Indicate call mitochondrion insertion")
    parser.add_option("--dna",
                      action="store_true", dest="dna", default=False,
                      help="Not RNA mediated insertion (no polyA)")
    parser.add_option("--cbs",
                      action="store_true", dest="cbs", default=False,
                      help="check by sample")#whether check by sample
    parser.add_option("--sva",
                      action="store_true", dest="sva", default=False,
                      help="For SVA insertion calling")
    parser.add_option("--gntp_feature",
                      action="store_true", dest="gntp_feature", default=False,
                      help="Collect genotyping features from bam")
    parser.add_option("--postF",
                      action="store_true", dest="postF", default=False,
                      help="Post filtering module")
    parser.add_option("--gntp_classify",
                      action="store_true", dest="gntp_classify", default=False,
                      help="Genotype every call of one sample (x_genotype_lik: beta-binomial "
                           "likelihood, biases fitted per sample by EM)")
    parser.add_option("--postFmosaic",
                      action="store_true", dest="postFmosaic", default=False,
                      help="Post filtering module for mosaic events")
    parser.add_option("--igv",
                      action="store_true", dest="igv", default=False,
                      help="Prepare screenshot command for given sites")
    parser.add_option("--case_control",
                      action="store_true", dest="case_control", default=False,
                      help="case control mode")
    parser.add_option("--tumor",
                      action="store_true", dest="tumor", default=False,
                      help="Working on tumor samples")
    # The germline-AF assumption is wrong for tumor AND for mosaic, in the same way and for the
    # same reason -- a subclonal/low-VAF insertion does not produce germline-strength read
    # ratios -- so the relaxation is named for what it means (somatic), not for one of the two
    # callers that needs it. --tumor still implies it, so existing tumor argv is unchanged.
    parser.add_option("--somatic",
                      action="store_true", dest="somatic", default=False,
                      help="Sample is somatic (tumor or mosaic): relax the germline "
                           "allele-fraction cutoffs in post-filtering")
    #convert to gVCF
    parser.add_option("--gVCF",
                      action="store_true", dest="gVCF", default=False,
                      help="Generate the gVCF from xTEA raw output")
####
    parser.add_option("--bed",
                      action="store_true", dest="bed", default=False,
                      help="Input annotation in bed format")
    parser.add_option("--flk_map",
                      action="store_true", dest="flk_map", default=False,
                      help="Map flanks to contigs")
    parser.add_option("--analysis",
                      action="store_true", dest="analysis", default=False,
                      help="Result analysis")
    parser.add_option("--flank", dest="flank", default=False,
                      help="flank regions")
    parser.add_option("--sv",
                      action="store_true", dest="sv", default=False,
                      help="Call promoted SVs")
    parser.add_option("--gene",
                      action="store_true", dest="gene", default=False,
                      help="Check whether the insertion falls in genes")
    parser.add_option("--somatic_hc",
                      action="store_true", dest="somatic_hc", default=False,
                      help="Get high confident somatic events")
    parser.add_option("--user",
                      action="store_true", dest="user_specific", default=False,
                      help="User specific parameters, by default automatically calc the parameters")
    parser.add_option("--single_sample",
                      action="store_true", dest="single_sample", default=False,
                      help="For single sample (like igv screenshot)")

#
####
    parser.add_option("-i", "--input", dest="input", default="",
                      help="input file ", metavar="FILE")
    parser.add_option("--input2", dest="input2", default="",
                      help="input file2 ", metavar="FILE")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file ", metavar="FILE")
    parser.add_option("-a", "--annotation", dest="annotation",
                      help="The annotation file ", metavar="FILE")
    # parser.add_option("-c", "--copies", dest="copies",
    #                   help="Repeat copies ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-d", "--barcode_bam", dest="barcode_bam",
                      help="Input barcode indexed bam file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string", default="./",
                      help="Working folder")
    parser.add_option("-n", "--cores", dest="cores", type="int", default=1,
                      help="number of cores")
    parser.add_option("-e", "--extend", dest="extend", type="int", default=0,
                      help="extend length")
    parser.add_option("-u", "--dup", dest="duplication",
                      help="duplication files", metavar="FILE")
    parser.add_option("--fflank", dest="fflank",
                      help="flank region file", metavar="FILE")
    parser.add_option("--flklen", dest="flklen", type="int",
                      help="flank region file")
    parser.add_option("--purity", dest="purity", type="float", default=0.45,#by default tumor purity set to 45%
                      help="Tumor purity")
    parser.add_option("--ref", dest="ref",
                      help="repeat consensus/copies", metavar="FILE")
    parser.add_option("--cr", dest="cliprep", type="int", default=1,
                      help="cutoff of minimum # of clipped parts fall in repeats")
    parser.add_option("--nd", dest="ndisc", type="int", default=5,
                      help="cutoff of minimum # of discordant pair")
    parser.add_option("--nb", dest="nbarcode", type="int", default=500,
                      help="cutoff of maximum # of molecure coverage")
    parser.add_option("--teilen", dest="teilen", type="int",
                      help="minimum length of the insertion for future analysis")
    parser.add_option("--iniclip", dest="iniclip", type="int", default=2,
                      help="initial minimum clip cutoff")
    parser.add_option("--af1", dest="af1", type="float", default=0.005,
                      help="minimal allel fraction")
    parser.add_option("--af2", dest="af2", type="float", default=0.45,
                      help="minimal allel fraction")
    parser.add_option("--rmsk_extnd", dest="rmsk_extnd", type="int", default=100,
                      help="Length of the left extended region when loading the repeatmasker output")
    parser.add_option("--rtype", dest="rep_type", type="int", default=1,
                      help="type of repeats: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-MIT, 32-MSTA")
    parser.add_option("--blacklist", dest="blacklist", default="null",
                      help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    parser.add_option("--type-map", dest="type_map", default=None,
                      help="xTea2 move A.2: chrm/pos/rtype sidecar (see locus_table.py) for --gVCF; "
                           "when set, writes one VCF per rtype present in the map instead of a single "
                           "global-rtype VCF", metavar="FILE")
    parser.add_option("--merged-vcf-out", dest="merged_vcf_out", default=None,
                      help="xTea2 move A.2: exact output path for the single merged VCF gVCF writes "
                           "in addition to the per-type ones when --type-map is set", metavar="FILE")
    # Trim step 4 removed the parse_option entries below along with the branches that
    # read them -- each had exactly one consumer, and that consumer was a deleted
    # branch. Kept as a list rather than silently dropped, because -C/-D in particular
    # appear in older run scripts and the useful answer to "why did -C stop working"
    # is here rather than in a commit message:
    #   barcode       (3 lines) -- -B/--barcode branch
    #   clip          (3 lines) -- -C branch
    #   cns           (2 lines) -- -C / --mosaic (sf_rep_cns)
    #   collect_clip  (3 lines) -- -Q/--collectclip branch
    #   cov           (2 lines) -- -C / --mosaic (max_cov_cutoff)
    #   cwfolder      (2 lines) -- -C / --mosaic / -Q (wfolder_pub_clip)
    #   discordant    (3 lines) -- -D branch
    #   force         (3 lines) -- -C / --mosaic (b_force)
    #   lclip         (2 lines) -- -C / --mosaic (cutoff_left_clip)
    #   mosaic        (3 lines) -- --mosaic branch
    #   rclip         (2 lines) -- -C / --mosaic (cutoff_right_clip)
    #   single        (3 lines) -- -C / --mosaic / -Q (b_se)
    #   siteclip      (2 lines) -- -C / --mosaic (site_clip_cutoff)
    #   somatic       (3 lines) -- -D branch (PEAK_WINDOW selection)
    # The scan-only replacements: --scan-min-clip (was --lc/--rc/--cr) and the
    # coverage-derived disc cutoff in run_xtea.py (was --nd). See run_xtea.py's
    # module docstring.
    (options, args) = parser.parse_args()
    return (options, args)
####
####
def automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False, b_tumor=False, f_purity=0.45):
    ####1. collect the basic information
    search_win = 500
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    rcd=x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)
    f_cov=rcd[0]

    ####2. based on the coverage, set the parameters
    xpar=Parameters()
    if b_tumor==True:
        f_cov=f_cov*f_purity
    par_rcd=xpar.get_par_by_cov(f_cov) #in format (iclip, idisc, i_clip-disc)
    print("Ave coverage is {0}: automatic parameters (clip, disc, clip-disc) with value ({1}, {2} ,{3})\n"
          .format(f_cov, par_rcd[0], par_rcd[1], par_rcd[2]))
    return par_rcd, rcd

####
def automatic_set_molecule_cutoff_for_10X_bam(sf_bam, sf_ref, s_working_folder, n_jobs):
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    search_win = 500
    f_cov=x_basic_info.calc_cov_from_bam(sf_bam, sf_ref, search_win)
    xpar=XParameters()
    f_molecule_cutoff=xpar.get_barcode_cov_cutoff(f_cov)
    return f_molecule_cutoff

def automatic_gnrt_parameters_case_control(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False):
    ####1. collect the basic information
    search_win = 500
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    rcd=x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)
    f_cov=rcd[0]

    ####2. based on the coverage, set the parameters
    xpar=CaseControlFilterPars()
    par_rcd=xpar.get_par_by_cov(f_cov) #in format (iclip, idisc, i_clip-disc)
    print("Ave coverage is {0}: automatic parameters (clip, disc, clip-disc) with value ({1}, {2} ,{3})\n"
          .format(f_cov, par_rcd[0], par_rcd[1], par_rcd[2]))
    return par_rcd, rcd

####
def adjust_cutoff_tumor(ncutoff=-1, i_adjust=1):
    if ncutoff-i_adjust>1:
        ncutoff=ncutoff-i_adjust
    return ncutoff

##main function
if __name__ == '__main__':
    (options, args) = parse_option()

    if options.mit:#if this to call mitochondrial insertion, then will not filter out chrM in "x_intermediate_sites.py"
        global_values.turn_on_mit()

    if options.dna:
        global_values.turn_off_rna_mediated()
    if options.cbs:
        global_values.turn_on_check_by_sample()
    # The family this invocation is running, resolved ONCE and passed explicitly to
    # every object that behaves differently per family -- replacing the process-wide
    # global_values.IS_CALL_SVA that used to be flipped here. `--sva` is the only family
    # signal the pre-Stage-2 steps have, and run_xtea.build_type_steps passes it to
    # cns_filter alone (matching the golden run_xTEA_pipeline.sh), so transduction and
    # sibling resolve GENERIC even when their --rtype says SVA. That is exactly what the
    # global did; see xtea/family_params.py for what it means and when to change it.
    x_family = family_params.from_sva_flag(options.sva)

    b_automatic=True
    if options.user_specific:
        b_automatic=False
    b_tumor=options.tumor #whether this is tumor sample
    f_purity=options.purity#tumor purity, by default 0.45
    b_resume=options.resume#resume the running, which will skip the step if output file already exist

    if options.preprocess:  #preprocess steps
        s_working_folder = options.wfolder
        sf_ref = options.reference
        sf_annotation = options.annotation
        sf_out_fa = options.output
        flank_lth = options.extend
        b_with_flank = options.withflank  # if not set then no flank region

        b_bed_fmt=options.bed
        x_annotation = XAnnotation(sf_annotation)
        i_lextnd=options.rmsk_extnd
        global_values.set_load_rmsk_left_extnd(i_lextnd)
        if b_bed_fmt==True:
            x_annotation.collect_seqs_of_TE_from_ref_bed_fmt(sf_ref, sf_out_fa, flank_lth)
        else:# this is for repeatmasker output
            b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
            x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
            x_annotation.load_rmsk_annotation()
            x_annotation.collect_seqs_of_TE_from_ref(sf_ref, flank_lth, b_with_flank, sf_out_fa)
            x_annotation.bwa_index_TE_seqs(sf_out_fa)

    elif options.flank:#preprocess the flank regions steps
        s_working_folder = options.wfolder
        sf_ref = options.reference
        sf_annotation = options.annotation
        sf_out_fa = options.output
        flank_lth = options.extend

        x_annotation = XAnnotation(sf_annotation)
        b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
        x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
        b_bed_fmt = options.bed
        if b_bed_fmt==True:
            print("load from bed")
            x_annotation.load_annotation_no_extnd_from_bed()
        else:
            print("load from rmsk")
            x_annotation.load_rmsk_annotation_no_extnd()

        x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa)  # only get the flank regions
        x_annotation.bwa_index_TE_seqs(sf_out_fa)

####
    # -C (clip) and -D (disc) are GONE -- trim step 4 committed xTea2 to the scan-only
    # path. The clip step is now the C++ xtea-scan single pass (scan_adapter) and the disc
    # step is a native pass over that same scan's --disc-pos stream (disc_adapter); both run
    # in run_xtea.py before this chain starts, so these two branches had no caller left.
    # They were the only users of TE_Multi_Locator's clip/disc entry points, i.e. of the
    # 1,315-line x_TEI_locator.py and the 1,062-line ClipReadInfo in clip_read.py, all
    # deleted with them. Recoverable from git, and an untouched master copy of the same code
    # lives at /home/cos689/xtea_exploration/xTea/xtea/ (what the cc3_master*_isolation.sh
    # attribution scripts already run against -- they never pointed at this checkout).
    #
    # NOTE for anyone reinstating a pysam fallback: don't. B9 (clipdisc spec section 4)
    # records that this code could not run under py3.10 at all -- x_annotation.py:582 did
    # `mid = (lo + hi) / 2` and indexed a list with the float, raising TypeError on every
    # call. The deleted branches were already dead in practice, not merely unused.
    ####
    ####
    elif options.filter_csn:  #filter out the FP by the pattern in the consensus repeat
        print("Working on \"clip-disc-filtering\" step!")
        sf_bam_list = options.bam  ###read in a bam list file
        s_working_folder = options.wfolder
        print("Current working folder is: {0}\n".format(s_working_folder))
        n_jobs = options.cores
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        sf_candidate_list = options.input#
        sf_raw_disc=options.input2#this is the raw disc file
        iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
        bin_size = 50000000  # block size for parallelization
        sf_cns = options.reference  ####repeat copies/cns here
        bmapped_cutoff = global_values.MIN_CLIP_MAPPED_RATIO
        sf_annotation = options.annotation
        i_concord_dist = 550  # this should be the 3*std_derivation, used to cluster disc reads on the consensus
        f_concord_ratio = 0.45
        sf_output = options.output
        sf_flank=options.fflank #this is the flanking region
        i_flank_lenth = options.flklen

        if b_resume == False or os.path.isfile(sf_output) == False:
            b_force=True
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force,
                                                       b_tumor, f_purity)
            ave_cov = basic_rcd[0]  # ave coverage
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard derivation
            print("Mean insert size is: {0}\n".format(mean_is))
            print("Standard derivation is: {0}\n".format(std_var))

            max_is = int(mean_is + 3 * std_var)
            if iextnd < max_is: #correct the bias
                iextnd = max_is
            if i_concord_dist < max_is: #correct the bias
                i_concord_dist = max_is
            global_values.set_read_length(rlth)
            global_values.set_insert_size(max_is)
            global_values.set_average_cov(ave_cov)
            print("Read length is: {0}\n".format(rlth))
            print("Maximum insert size is: {0}\n".format(max_is))
            print("Average coverage is: {0}\n".format(ave_cov))

            n_clip_cutoff = options.cliprep #this is the sum of left and right clipped reads
            n_disc_cutoff = options.ndisc  #each sample should have at least this number of discordant reads
            if b_automatic==True:
                n_clip_cutoff=rcd[0]
                n_disc_cutoff=rcd[1]
                # if b_tumor==True:
            print("Filter (on cns) cutoff: {0} and {1} are used!!!".format(n_clip_cutoff, n_disc_cutoff))

            # Phase C move B.1: emit temp_clip.sam/temp_disc.sam/low_mapq_clip_cnt.txt from a native
            # single-pass scan (docs/phase_c_unified_downstream.md sec 7), reusing THIS step's iextnd
            # so the collection window matches the cns math exactly. Opt-in via XTEA_SCAN_CNS_EMIT=<scan
            # bin>; on success sets XTEA_CNS_PRECOMPUTED_ALIGN so call_MEIs_consensus below reuses the
            # emitted SAMs instead of re-reading the BAM + re-running bwa.
            s_scan_bin = os.environ.get("XTEA_SCAN_CNS_EMIT", "")
            if s_scan_bin:
                import subprocess as _sp
                wf = s_working_folder if s_working_folder.endswith("/") else s_working_folder + "/"
                with open(sf_bam_list) as _fbl:
                    _bams = [ln.split()[0] for ln in _fbl if ln.strip()]
                if len(_bams) != 1:
                    print("[cns-emit] only single-BAM supported for now; using normal collect")
                else:
                    s_bwa = os.environ.get("XTEA_SCAN_BWA", global_values.BWA_PATH)
                    scan_cmd = [s_scan_bin, "--bam", _bams[0], "--reference", sf_ref,
                                "--cns-candidates", sf_candidate_list, "--cns-ref", sf_cns,
                                "--cns-clip-sam", wf + "temp_clip.sam",
                                "--cns-disc-sam", wf + "temp_disc.sam",
                                "--cns-lowq", wf + "low_mapq_clip_cnt.txt",
                                "--cns-extnd", str(iextnd), "--cns-sample", "0",
                                "--bwa", s_bwa, "--bwa-threads", str(n_jobs)]
                    print("[cns-emit] " + " ".join(scan_cmd))
                    _sp.run(scan_cmd, check=True)
                    os.environ["XTEA_CNS_PRECOMPUTED_ALIGN"] = "1"

            x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref,
                                          family=x_family)
            x_cd_filter.call_MEIs_consensus(sf_candidate_list, sf_raw_disc, iextnd, bin_size, sf_cns, sf_flank,
                                            i_flank_lenth, bmapped_cutoff, i_concord_dist, f_concord_ratio,
                                            n_clip_cutoff, n_disc_cutoff, sf_output)
    ####
####
    #why a new sub-module?
    #1. we need to collect clip and disc reads for all the candidates (filter out some existing and not qualified ones)
    #2. The old module doesn't work well!
    #3. we need to select the candidate
    elif options.transduction:#
        #need to re-collect all the clip, disc reads
        sf_bam_list = options.bam  ###read in a bam list file
        sf_candidate_list = options.input #this is the output from the "cns" step.
        sf_raw_disc = options.input2  # this is the raw disc file
        iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
        bin_size = 50000000  # block size for parallelization
        sf_cns = options.reference  ####repeat copies/cns here
        s_working_folder = options.wfolder
        print("Current working folder is: {0}\n".format(s_working_folder))
        n_jobs = options.cores
        sf_reference = options.ref  ###reference genome, some cram file require this file to open
        sf_flank = options.fflank  # this is the flanking region
        i_flank_lenth = options.flklen
        sf_output = options.output
        sf_rmsk = options.annotation
        i_rep_type = options.rep_type
####
        if b_resume == False or os.path.isfile(sf_output) == False:
            if os.path.isfile(sf_flank)==True:#for Alu and many others, there is no transduction
                b_force = False
                rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_reference, s_working_folder, n_jobs, b_force,
                                                           b_tumor, f_purity)
                ave_cov = basic_rcd[0]  # ave coverage
                rlth = basic_rcd[1]  # read length
                mean_is = basic_rcd[2]  # mean insert size
                std_var = basic_rcd[3]  # standard derivation
                print("Mean insert size is: {0}\n".format(mean_is))
                print("Standard derivation is: {0}\n".format(std_var))
                max_is = int(mean_is + 3 * std_var)
                if iextnd < max_is:  # correct the bias
                    iextnd = max_is
                i_concord_dist=550
                f_concord_ratio = 0.25
                if i_concord_dist < max_is:  # correct the bias
                    i_concord_dist = max_is
                global_values.set_read_length(rlth)
                global_values.set_insert_size(max_is)
                global_values.set_average_cov(ave_cov)
####
                n_clip_cutoff = options.cliprep  # this is the sum of left and right clipped reads
                n_disc_cutoff = options.ndisc  # each sample should have at least this number of discordant reads
                if b_automatic == True:
                    n_clip_cutoff = rcd[0]
                    n_disc_cutoff = rcd[1]

                xtransduction = XTransduction(s_working_folder, n_jobs, sf_reference,
                                              family=x_family)

                i_win=150 #if a site is close to an existing site, then will not be considered again
                sf_tmp_slct=sf_raw_disc+".slct"
                i_min_copy_len = 225
                xannotation = xtransduction.prep_annotation_interval_tree(sf_rmsk, i_min_copy_len)
                xtransduction.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_candidate_list, n_disc_cutoff, xannotation,
                                             i_rep_type, b_tumor, sf_tmp_slct)
                #now for the selected sites, re-evaluate each one
                x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_reference,
                                              family=x_family)
                i_max_cov=ave_cov*(global_values.MAX_COV_TIMES+1)
                sf_output_tmp=sf_output + global_values.TD_NON_SIBLING_SUFFIX
                xtransduction.call_candidate_transduction_v3(sf_tmp_slct, sf_candidate_list, x_cd_filter,
                                                             sf_flank, sf_cns, i_flank_lenth, iextnd, bin_size, n_clip_cutoff,
                                                             n_disc_cutoff, i_concord_dist, f_concord_ratio, xannotation,
                                                             sf_bam_list, i_rep_type, i_max_cov, ave_cov, sf_output_tmp)
####
                xorphan=XOrphanTransduction(s_working_folder, n_jobs, sf_reference,
                                            family=x_family)
                n_half_disc_cutoff=n_disc_cutoff/2
                i_search_win=2000
                sf_updated_cns=sf_output #this is the final updated

                #1.Call out the sibling transduction events from the current list
                sf_sibling_TD=sf_output+".sibling_transduction_from_existing_list"#
                xorphan.call_sibling_TD_from_existing_list(sf_output_tmp, sf_bam_list, iextnd, n_half_disc_cutoff,
                                                           i_search_win, xannotation, i_rep_type, i_max_cov,
                                                           sf_updated_cns, sf_sibling_TD)
#
                # #2. Call orphan "sibling" transdcution from non_existing list
                # #select the sites to exclude the already called out sites
                # xorpha.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_output_tmp, n_disc_cutoff, xannotation,
                #                                                i_rep_type, b_tumor, sf_tmp_slct2, b_with_original)
####
                #update high confident ones (in "cns" filter step, two results are generated)
                sf_ori_hc=sf_candidate_list+global_values.HIGH_CONFIDENT_SUFFIX
                sf_new_hc=sf_output+global_values.HIGH_CONFIDENT_SUFFIX
                xorphan.update_high_confident_callset(sf_ori_hc, sf_updated_cns, sf_new_hc)
####
            else:#rename the two files generated in previous step
                copyfile(sf_candidate_list, sf_output)
                sf_ori_hc = sf_candidate_list + global_values.HIGH_CONFIDENT_SUFFIX
                sf_new_hc=sf_output+global_values.HIGH_CONFIDENT_SUFFIX
                copyfile(sf_ori_hc, sf_new_hc)

####
    elif options.sibling:#sibling orphan transduction
        '''
        Todo: 09-29-2019: Add filtering modules:
        1. using background low mapq reads (multiple mapped reads) for filtering
        2. Set a upper-bound cutoff for discordant reads
        3. using blacklist for filtering
        '''
        sf_bam_list = options.bam  ###read in a bam list file
        sf_pre_step_out = options.input  # this is the output from the "cns" step.
        sf_raw_disc = options.input2  # this is the raw disc file
        iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
        bin_size = 50000000  # block size for parallelization
        sf_cns = options.reference  ####repeat copies/cns here
        s_working_folder = options.wfolder
        print("Current working folder is: {0}\n".format(s_working_folder))
        n_jobs = options.cores
        sf_reference = options.ref  ###reference genome, some cram file require this file to open
        sf_flank = options.fflank  # this is the flanking region
        i_flank_lenth = options.flklen
        sf_output = options.output
        sf_rmsk = options.annotation
        i_rep_type = options.rep_type
        sf_black_list = options.blacklist

        if os.path.isfile(sf_flank) == True:  #for Alu and many others, there is no transduction
            b_force = False
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_reference, s_working_folder, n_jobs, b_force,
                                                       b_tumor, f_purity)
            ave_cov = basic_rcd[0]  # ave coverage
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard derivation
            print("Mean insert size is: {0}\n".format(mean_is))
            print("Standard derivation is: {0}\n".format(std_var))
            max_is = int(mean_is + 3 * std_var)

            i_concord_dist = 550
            f_concord_ratio = 0.25
            if i_concord_dist < max_is:  # correct the bias
                i_concord_dist = max_is
            global_values.set_read_length(rlth)
            global_values.set_insert_size(mean_is) #here set mean inset size
            global_values.set_average_cov(ave_cov)

            n_clip_cutoff = options.cliprep  # this is the sum of left and right clipped reads
            n_disc_cutoff = options.ndisc  # each sample should have at least this number of discordant reads
            if b_automatic == True:
                n_clip_cutoff = rcd[0]
                n_disc_cutoff = rcd[1]

            xorphan = XOrphanTransduction(s_working_folder, n_jobs, sf_reference,
                                          family=x_family)
            i_min_copy_len = 225
            xorphan.set_boundary_extend(mean_is)
            xannotation = xorphan.prep_annotation_interval_tree(sf_rmsk, i_min_copy_len)
            # 2. Call orphan "sibling" transdcution from non_existing list
            b_with_original = False
            sf_tmp_slct2 = sf_raw_disc + ".slct2"
            sf_output_tmp = sf_pre_step_out + global_values.TD_NON_SIBLING_SUFFIX
            # select the sites to exclude the already called out sites, and filter out sites fall in black_list
            xorphan.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_output_tmp, n_disc_cutoff, xannotation,
                                                    i_rep_type, b_tumor, sf_tmp_slct2, b_with_original)

            #re-select transduction candidates based on disc-clip consistency (clip position encompass disc ones?)
            m_failed_ori_td=xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list, i_concord_dist,
                                                                           n_clip_cutoff, n_disc_cutoff, sf_black_list,
                                                                           sf_rmsk)
            sf_sibling_TD=sf_output
            if os.path.isfile(sf_black_list)==False:
                print("Blacklist file {0} does not exist!".format(sf_black_list))
            xorphan.call_novel_sibling_TD_from_raw_list(sf_tmp_slct2, sf_bam_list, i_concord_dist, n_clip_cutoff,
                                                        n_disc_cutoff, sf_black_list, sf_rmsk, sf_sibling_TD)

            ####append the newly called events to existing list
            xorphan.append_to_existing_list(sf_sibling_TD, sf_pre_step_out, m_failed_ori_td)
            xorphan.append_to_existing_list(sf_sibling_TD, sf_pre_step_out+global_values.HIGH_CONFIDENT_SUFFIX,
                                            m_failed_ori_td)
        else:
            if os.path.isfile(sf_pre_step_out)==True:#do td filtering only
                b_force = False
                rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_reference, s_working_folder, n_jobs, b_force,
                                                           b_tumor, f_purity)
                ave_cov = basic_rcd[0]  # ave coverage
                rlth = basic_rcd[1]  # read length
                mean_is = basic_rcd[2]  # mean insert size
                std_var = basic_rcd[3]  # standard derivation
                print("Mean insert size is: {0}\n".format(mean_is))
                print("Standard derivation is: {0}\n".format(std_var))
                max_is = int(mean_is + 3 * std_var)

                i_concord_dist = 550
                f_concord_ratio = 0.25
                if i_concord_dist < max_is:  # correct the bias
                    i_concord_dist = max_is
                global_values.set_read_length(rlth)
                global_values.set_insert_size(mean_is)  # here set mean inset size
                global_values.set_average_cov(ave_cov)
                n_clip_cutoff = options.cliprep  # this is the sum of left and right clipped reads
                n_disc_cutoff = options.ndisc  # each sample should have at least this number of discordant reads
                if b_automatic == True:
                    n_clip_cutoff = rcd[0]
                    n_disc_cutoff = rcd[1]
                xorphan = XOrphanTransduction(s_working_folder, n_jobs, sf_reference,
                                              family=x_family)
                # re-select transduction candidates based on disc-clip consistency (clip position encompass disc ones?)
                m_failed_ori_td = xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list,
                                                                                   i_concord_dist, n_clip_cutoff,
                                                                                   n_disc_cutoff, sf_black_list, sf_rmsk)
                xorphan.update_existing_list_only(sf_pre_step_out, m_failed_ori_td)
####
####
    ####this module for: 1) tumor case-control files;
    elif options.case_control:#case-control mode to call somatic events
        b_somatic_hc = options.somatic_hc
        sf_candidate_list = options.input  # this is the list called from case
        sf_output = options.output
        s_working_folder = options.wfolder
        if b_somatic_hc==False:#this is for all the raw call set
            sf_bam_list = options.bam  # this is the control bam file list
            sf_ref = options.ref  # reference genome
            n_jobs = options.cores
            nclip_cutoff = options.cliprep  # this is the sum of left and right clipped reads
            ndisc_cutoff = options.ndisc  # each sample should have at least this number of discordant reads
            sf_rep_cns = options.reference  ####repeat copies/cns here
            sf_flank = options.fflank  # this is the flanking region
            i_flk_len = options.flklen

            b_force=True
            rcd=None
            basic_rcd=None
            n_polyA_cutoff=0
            if b_automatic==True:
                rcd, basic_rcd=automatic_gnrt_parameters_case_control(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force)
                nclip_cutoff=rcd[0]
                ndisc_cutoff=rcd[1]
                n_polyA_cutoff=rcd[2]

            ccm=CaseControlMode(sf_ref, s_working_folder, n_jobs)
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            global_values.set_insert_size(int(mean_is))
            std_var = basic_rcd[3]  # standard derivation
            max_is = int(mean_is + 3 * std_var) + int(rlth)
            extnd = max_is
            bin_size = 50000000  # block size for parallelization
            print("clip,disc,polyA-cutoff is ({0}, {1}, {2})".format(nclip_cutoff, ndisc_cutoff, n_polyA_cutoff))
            n_polyA_cutoff=ndisc_cutoff #if both sides have more than cutoff polyA, then filter out
            ccm.call_somatic_TE_insertion(sf_bam_list, sf_candidate_list, extnd, nclip_cutoff, ndisc_cutoff,
                                          n_polyA_cutoff, sf_rep_cns, sf_flank, i_flk_len, bin_size, sf_output, b_tumor)
    ####
        else:#This is to parse out the high confident somatic ones, assume already have the raw somatic callset
            sf_raw_somatic=options.input2
            sf_ref=""
            n_jobs=1
            ccm = CaseControlMode(sf_ref, s_working_folder, n_jobs)
            ccm.parse_high_confident_somatic(sf_candidate_list, sf_raw_somatic, sf_output)
####

####
    # --mosaic (master's mosaic clip step) is gone with trim step 4: it called
    # TE_Multi_Locator.call_TEI_candidate_sites_from_multiple_alignmts_mosaic, part of the
    # deleted pysam clip collector. This does NOT block MOSAIC_PLAN.md: its M.-1 reference
    # run is explicitly against master's own checkout in the `xTea` env, and its M.1 replaces
    # this path with a mosaic_adapter.py reading xtea-scan's polyA tally rather than calling
    # master's slow Python mosaic clip function at all.
    # --barcode (10X barcode-coverage filter, the only TELocator user) is gone with trim
    # step 4. xTea2 does not support 10X input.
    elif options.postF:#post filtering step
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_xtea_rslt = options.input
        sf_rmsk = options.annotation
        sf_new_out = options.output
        i_rep_type=options.rep_type
        sf_black_list=options.blacklist

        i_min_copy_len=225 #when check whether fall in repeat region, require the minimum copy length
        b_pf_mosaic=options.postFmosaic
        if b_pf_mosaic is True:#for mosaic events
            xpf_mosic = MosaicCaller(s_working_folder, n_jobs)
            xpf_mosic.run_call_mosaic(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, sf_black_list, sf_new_out)
        else:
            xbasic=X_BasicInfo(s_working_folder)
            f_cov, f_ave_cov=xbasic.load_get_cov()
            xpost_filter = XPostFilter(s_working_folder, n_jobs)
            # --tumor and --somatic both select the relaxed (non-germline) allele-fraction
            # cutoffs; --tumor is kept as an alias so the existing tumor/case-control argv keeps
            # working unchanged. run_post_filtering's parameter is still spelled b_tumor because
            # it threads down into master's post_processing_* functions untouched.
            b_somatic = b_tumor or options.somatic
            xpost_filter.run_post_filtering(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, f_cov, sf_black_list,
                                            sf_new_out, b_somatic)
####
    ####
    elif options.gntp_feature:#generate the genotype features
        sf_bam_list = options.bam
        sf_ref = options.ref
        sf_candidate_list = options.input
        n_jobs = options.cores
        s_working_folder = options.wfolder
        sf_output = options.output

        x_gntper = XGenotyper(sf_ref, s_working_folder, n_jobs)
        extnd = 450
        x_gntper.call_genotype(sf_bam_list, sf_candidate_list, extnd, sf_output)

    elif options.gntp_classify:
        # One EM fit over every call in the input, so the input must hold ALL of one sample's
        # calls (the merged, all-family list in the shared tail). Appends GT, GQ, GTL, AF_SR,
        # AF_PE, VAF, VAF_CI as fields 50-56 and writes the fitted parameters to
        # <output>.genotype_params, which --gVCF puts in the VCF header.
        x_genotype_lik.genotype_raw_file(options.input, options.output)

####
    elif options.gVCF:
        sf_ref = options.ref
        sf_bam_list = options.bam  ###read in a bam list file
        sf_raw_rslt = options.input
        i_rep_type=options.rep_type##
        sf_prefix=options.output
        sf_bam=""
        s_sample_id="null"

        with open(sf_bam_list) as fin_blist:
            for line in fin_blist:
                fields=line.split()
                if len(fields)<=0:
                    break
                sf_bam=fields[0]
                bam_fields=sf_bam.split("/")
                if len(bam_fields)<=0:
                    break
                s_sample_id2=bam_fields[-1]
                tmp_fields=s_sample_id2.split(".")
                s_sample_id=".".join(tmp_fields[:-1])
                break
        if os.path.isfile(sf_bam)==True:
            gvcf=gVCF()
            if sf_prefix[-1]!="/":
                sf_prefix+="/"
            reptype = RepType()
            if options.type_map:
                # xTea2 move A.2: --rtype is the OR'd bitmask of every type present in the
                # unified locus set (e.g. 7 = L1|Alu|SVA); decompose it and write one VCF
                # per type, with each record's own type looked up from the sidecar (built
                # by locus_table.assign_and_write) instead of the single global --rtype.
                l_rtypes = [b for b in (1, 2, 4, 8, 16, 32, 64) if i_rep_type & b]
                m_sf_vcf = {rt: sf_prefix + s_sample_id + "_" + reptype.get_rep_type(rt) + ".vcf"
                            for rt in l_rtypes}
                gvcf.cvt_raw_rslt_to_gvcf_multi(s_sample_id, sf_bam, sf_raw_rslt, options.type_map,
                                                sf_ref, m_sf_vcf, options.merged_vcf_out)
            else:
                s_rep_type = reptype.get_rep_type(i_rep_type)
                sf_vcf=sf_prefix+s_sample_id+"_"+s_rep_type+".vcf"
                gvcf.cvt_raw_rslt_to_gvcf(s_sample_id, sf_bam, sf_raw_rslt, i_rep_type, sf_ref, sf_vcf)
        else:
            print("Wrong bam file: {0}".format(sf_bam))
####
    elif options.joint:
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_rslt_list = options.input
        sf_out=options.output
        islack=50
        if options.postFmosaic==True:
            mzc_joint=MosaicJointCalling(s_working_folder, n_jobs)
            mzc_joint.call_mosaic_from_multi_samples(sf_rslt_list, islack, sf_out)
####
    elif options.igv:#prepare the igv screenshot script for multiple bams and sites
        sf_sites = options.input #site list
        sf_bam_list=options.bam #bam list
        s_screenshot_folder = options.wfolder
        sf_out = options.output #output file
        sf_gnm=options.ref #"hg19" or "hg38"
        i_extnd = options.extend #"-e", "--extend"
        b_single_sample=options.single_sample

        x_igv = XIGV()#
        if b_single_sample==True:
            sf_sites=x_igv.gnrt_sites_single_sample(sf_sites, sf_bam_list)
        x_igv.prepare_igv_scripts_multi_bams(sf_sites, sf_bam_list, s_screenshot_folder, i_extnd, sf_gnm, sf_out)

####
    elif options.collect:  # collect the reads for each candidate site
        sf_ori_bam = options.bam
        sf_barcode_bam = options.barcode_bam
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_candidate_list = options.input
        i_cov_cutoff = options.nbarcode
        sf_annotation = options.annotation
        i_flank_lenth = options.flklen #extend the boundary of each annotated repeat
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        # collect reads for all the sites
        i_extend = 1500 #by default, collect barcode in [-1500, 1500] region
        xread_collection=XReadsCollection(s_working_folder, sf_ref)
        xread_collection.collect_phased_reads_all_TEIs(sf_ori_bam, sf_barcode_bam, sf_candidate_list, i_extend,
                                                       i_cov_cutoff, sf_annotation, i_flank_lenth, n_jobs)
####
    elif options.mutation:#call out the internal mutations by aligning the reads
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_cns = options.reference  ####repeat copies/cns here
        sf_sites=options.input
        sf_merged_vcf=options.output
        n_len_cutoff=options.teilen
        xmutation=XMutation(s_working_folder)
        xmutation.call_mutations_from_reads_algnmt(sf_sites, sf_cns, n_len_cutoff, n_jobs, sf_merged_vcf)

####
    elif options.gene:#
        #To run: --gene -a -i -o
        sf_gene_annotation=options.annotation
        sf_input=options.input
        sf_output=options.output
        gff=GFF3(sf_gene_annotation)
        iextnd=global_values.UP_DOWN_GENE
        i_user_extnd = options.extend
        if i_user_extnd>iextnd:
            iextnd=i_user_extnd
        gff.load_gene_annotation_with_extnd(iextnd)
        gff.index_gene_annotation_interval_tree()
        gff.annotate_results(sf_input, sf_output)

    # this step needs to retrieve the seqs of mate reads
    elif options.collect_illumina:  # this is only for normal illumina data
        sf_ori_bam = options.bam
        sf_candidate_list = options.input
        n_jobs = options.cores
        s_working_folder = options.wfolder

    elif options.assembly:  # assemble the reads for all the sites (this is for 10X module)
        b_local = options.local
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_candidate_list = options.input
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        xlasm = XLocalAssembly(s_working_folder, sf_ref)
        if b_local == True:
            xlasm.assemble_all_phased_TEIs_locally(sf_candidate_list, n_jobs)
        else:
            xlasm.assemble_all_phased_TEIs_slurm_cluster(sf_candidate_list)

####
####
    elif options.map:####align the asm to reference genome
        sf_ref = options.reference
        s_working_folder = options.wfolder
        n_jobs = int(options.cores)
        sf_sites = options.input
        sf_cns=options.ref####repeat consensus
        sf_final_list=options.output
        sf_final_seqs=sf_final_list+".fa"

        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.align_asm_contigs_to_reference(sf_sites, sf_ref, s_working_folder)
        xctg.call_MEIs_from_all_group_contigs(sf_sites, sf_ref, sf_cns, s_working_folder, sf_final_list, sf_final_seqs)

    elif options.flk_map:  #gnrt the flank regions and align the regions to the contigs
        sf_ref = options.ref
        s_working_folder = options.wfolder
        n_jobs = int(options.cores)
        sf_sites = options.input

        i_extend = 500
        xref = XReference()
        xref.gnrt_flank_region_for_sites(sf_sites, i_extend, n_jobs, sf_ref, s_working_folder)
        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.align_flanks_to_phased_contig(sf_sites)
####
    elif options.filter_asm:##
        s_working_folder = options.wfolder
        n_jobs = int(options.cores)
        sf_sites = options.input
        sf_keep_sites = options.output
        sf_repeat_copies = options.reference  ####repeat copies here

        flank_length = 500
        f_map_cutoff = 0.45 #25% of the bases are required to be mapped! This is a small ratio
        i_slack = 35
        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.validate_TEI_from_phased_asm_algnmt(sf_sites, flank_length, f_map_cutoff, i_slack, sf_repeat_copies,
                                                 sf_keep_sites)

    # -Q/--collectclip is gone with trim step 4: it was a thin wrapper over
    # TE_Multi_Locator.collect_all_clipped_from_multiple_alignmts, i.e. the pysam clip
    # collector this step deleted. Nothing in xTea2 invoked it.
