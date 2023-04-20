"""
Consolidate truvari tp/fn/fp results and annotate with GA4GH intermediates tags
"""
import os
import sys
import logging
import argparse

import pysam
import truvari

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
    header.add_line(('##FORMAT=<ID=BD,Number=1,Type=String,'
                     'Description="Decision for call (TP/FP/FN/N)">'))
    header.add_line(('##FORMAT=<ID=BK,Number=1,Type=String,'
                     'Description="Sub-type for decision (match/mismatch type) lm variants had to which it could compare">'))
    return header

def consolidate_annotate(pos_fn, neg_fn, out_fn, bsamp=0, csamp=0):
    """
    Consolidate the two truvari results
    """
    f_key = 'FN' if 'base' in pos_fn else 'FP'
    pos = pysam.VariantFile(pos_fn)
    neg = pysam.VariantFile(neg_fn)
    n_header = edit_header(pos.header.copy())
    out = pysam.VariantFile(out_fn, 'w', header=n_header)
    ts = 0
    for entry in pos:
        ts += 1
        entry.translate(n_header)
        entry.samples[bsamp]["BD"] = "TP"
        entry.samples[bsamp]["BK"] = 'gm' if entry.info["GTMatch"] == 0 else 'am'
        out.write(entry)

    fs = 0
    for entry in neg:
        fs += 1
        entry.translate(n_header)
        entry.samples[csamp]["BD"] = f_key
        entry.samples[csamp]["BK"] = '.' if entry.info["TruScore"] is None else 'lm'
        out.write(entry)

    out.close()
    logging.info("Total of %d variants (%d:TP, %d:%s)", ts + fs, ts, fs, f_key)
    return

def truv2ga4gh_main(args):
    """
    """
    args = parse_args(args)
    truvari.setup_logging()
    logging.info("Consolidating Truth VCF")
    t_vcf = args.output + '_truth.vcf'
    consolidate_annotate(os.path.join(args.input, "tp-base.vcf.gz"),
                         os.path.join(args.input, "fn.vcf.gz"),
                         t_vcf, args.bSample, args.cSample)
    truvari.compress_index_vcf(t_vcf)
    
    logging.info("Consolidating Query VCF")
    q_vcf = args.output + '_query.vcf'
    consolidate_annotate(os.path.join(args.input, "tp-comp.vcf.gz"),
                         os.path.join(args.input, "fp.vcf.gz"),
                         q_vcf, args.bSample, args.cSample)
    truvari.compress_index_vcf(q_vcf)

if __name__ == '__main__':
    truv2ga4gh_main(sys.argv[1:])
