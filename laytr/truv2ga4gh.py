"""
Consolidate truvari tp/fn/fp results and annotate with GA4GH intermediates tags
"""
import os
import sys
import json
import logging
import argparse
from collections import defaultdict

import pysam
import truvari
import pandas as pd
from intervaltree import IntervalTree

def parse_args(args):
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(prog="truv2ga4gh", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True,
                        help="Truvari result directory")
    parser.add_argument("-o", "--output", required=True,
                        help="Output suffix")
    parser.add_argument("-w", "--with-refine", action="store_true",
                        help="Consolidate with `truvari refine` output")
    parser.add_argument("-b", "--bSample", default=0,
                        help="Sample name to annotate in truth VCF (default first)")
    parser.add_argument("-c", "--cSample", default=0,
                        help="Sample name to annotate in query VCF (default first)")

    args = parser.parse_args(args)
    return args
    
def edit_header(header):
    """
    Add INFO for new fields to vcf
    """
    header.add_line(('##INFO=<ID=IsRefined,Number=0,Type=Flag,'
                     'Description="True if variant was pulled from refinement">'))
    header.add_line(('##FORMAT=<ID=BD,Number=1,Type=String,'
                     'Description="Decision for call (TP/FP/FN/N)">'))
    header.add_line(('##FORMAT=<ID=BK,Number=1,Type=String,'
                     'Description="Sub-type for decision (match/mismatch type) lm variants had any counterpart to which it could compare">'))
    return header

def in_tree(tree, entry):
    """
    return True if entry within boundaries of tree interval
    """
    astart, aend = truvari.entry_boundaries(entry)
    # Filter these early so we don't have to keep checking overlaps
    if astart == aend - 1:
        return tree[entry.chrom].overlaps(astart)
    m_ovl = tree[entry.chrom].overlap(astart, aend)
    if len(m_ovl) != 1:
        return False
    m_ovl = list(m_ovl)[0]
    # Edge case - the variant spans the entire include region
    return astart >= m_ovl.begin and aend <= m_ovl.end

def build_tree(regions):
    """
    Build tree from regions
    """
    tree = defaultdict(IntervalTree)
    for _, i in regions.iterrows():
        tree[i['chrom']].addi(i['start'], i['end'] + 1)
    return tree
        
def get_truvari_filenames(in_dir):
    return {'tp-base': os.path.join(in_dir, 'tp-base.vcf.gz'),
            'tp-comp': os.path.join(in_dir, 'tp-comp.vcf.gz'),
            'fn': os.path.join(in_dir, 'fn.vcf.gz'),
            'fp': os.path.join(in_dir, 'fp.vcf.gz'),
            'params': os.path.join(in_dir, 'params.json')}

def pull_entries(in_vcf_fn, keep_tree=None, is_refined=False):
    """
    pull variants inside regions (or all if None)
    """
    in_vcf = pysam.VariantFile(in_vcf_fn)
    m_iter = in_vcf if keep_tree is None else truvari.region_filter(in_vcf, keep_tree, is_refined)
    for entry in m_iter:
        yield entry

def move_record(entry, out_vcf, sample):
    """
    Move the new entry (with the single specified sample)
    to be writeable to the out_vcf
    """
    n_sampname = out_vcf.header.samples[0]
    ret = out_vcf.new_record(contig=entry.contig, start=entry.start, stop=entry.stop, alleles=entry.alleles,
                              id=entry.id, qual=entry.qual, filter=entry.filter, info=entry.info)
    for k,v in entry.samples[sample].items():
        ret.samples[0][k] = v
    return ret

def parse_bench_dir(in_dir, t_vcf, q_vcf, tree=None, is_refined=False):
    """
    Pull relevant entries from relevant files
    """
    in_vcfs = get_truvari_filenames(in_dir)
    params = json.load(open(in_vcfs['params'], 'r'))
    bsamp = params['bSample']
    csamp = params['cSample']
    ops_to_do = [("tp-base", "TP", t_vcf, bsamp), ("fn", "FN", t_vcf, bsamp),
                 ("tp-comp", "TP", q_vcf, csamp), ("fp", "FP", q_vcf, csamp)]
    for filename, bdkey, out_vcf, samp in ops_to_do:
        for entry in pull_entries(in_vcfs[filename], tree, is_refined):
            n_entry = move_record(entry, out_vcf, samp)
            n_entry.info["IsRefined"] = is_refined
            n_entry.samples[0]["BD"] = bdkey
            if bdkey.startswith('T'):
                n_entry.samples[0]["BK"] = 'gm' if entry.info["GTMatch"] == 0 else 'am'
            else:
                n_entry.samples[0]["BK"] = '.' if entry.info["TruScore"] is None else 'lm'
            out_vcf.write(n_entry)

def check_bench_dir(dirname):
    """
    Make sure all the files are there
    """
    check_fail = False
    b_files = get_truvari_filenames(dirname)
    for k, v in b_files.items():
        if not os.path.exists(v):
            logging.error("%s bench file doesn't exist", k)
            check_fail = True
    return check_fail


def check_args(args):
    """
    Ensure everything we're looking for is available
    return True if check failed
    """
    check_fail = False
    if not os.path.isdir(args.input):
        logging.error("input is not a directory")
        check_fail = True
    else:
        check_fail |= check_bench_dir(args.input)

    if os.path.exists(args.output + '_truth.vcf.gz'):
        logging.error("%s_truth.vcf.gz already exists", args.output)
        check_fail = True
    if os.path.exists(args.output + '_truth.vcf.gz.tbi'):
        logging.error("%s_truth.vcf.gz.tbi already exists", args.output)
        check_fail = True
    if os.path.exists(args.output + '_query.vcf.gz'):
        logging.error("%s_query.vcf.gz already exists", args.output)
        check_fail = True
    if os.path.exists(args.output + '_query.vcf.gz.tbi'):
        logging.error("%s_query.vcf.gz.tbi already exists", args.output)
        check_fail = True
    
    if not args.with_refine:
        return check_fail

    pdir = os.path.join(args.input, 'phab_bench')
    if not os.path.exists(pdir):
        logging.error("phab_bench dir doesn't exist")
        check_fail = True
    else:
        check_fail |= check_bench_dir(pdir)
    
    return check_fail

def truv2ga4gh_main(args):
    """
    """
    args = parse_args(args)
    truvari.setup_logging()

    if check_args(args):
        logging.error("Unable to do conversion")
        sys.exit(1)

    logging.info("Consolidating VCFs")
    
    bench_vcfs = get_truvari_filenames(args.input)

    b_header = edit_header(pysam.VariantFile(bench_vcfs['tp-base'], 'r').header.copy())
    t_vcf_fn = args.output + '_truth.vcf'
    t_vcf = pysam.VariantFile(t_vcf_fn, 'w', header=b_header)

    c_header = edit_header(pysam.VariantFile(bench_vcfs['tp-comp'], 'r').header.copy())
    q_vcf_fn = args.output + '_query.vcf'
    q_vcf = pysam.VariantFile(q_vcf_fn, 'w', header=c_header)

    regions = None
    refine_tree = None
    if args.with_refine:
        regions = pd.read_csv(os.path.join(args.input, "refine.regions.txt"), sep='\t')
        refine_tree = build_tree(regions[regions['refined']])
        parse_bench_dir(os.path.join(args.input, 'phab_bench'), t_vcf, q_vcf,
                        tree=refine_tree, is_refined=True)
        
    parse_bench_dir(args.input, t_vcf, q_vcf, tree=refine_tree, is_refined=False)
    t_vcf.close()
    q_vcf.close()
    truvari.compress_index_vcf(t_vcf_fn)
    truvari.compress_index_vcf(q_vcf_fn)
    logging.info("Finished")

if __name__ == '__main__':
    truv2ga4gh_main(sys.argv[1:])
