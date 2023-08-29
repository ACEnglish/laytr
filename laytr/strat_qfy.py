import argparse

import pysam
import truvari
import pandas as pd


def parse_args(args):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description="Calculate variant metrics")
    parser.add_argument("--base_vcf", required=True, help="Path to the base VCF file")
    parser.add_argument("--comp_vcf", required=True, help="Path to the comp VCF file")
    parser.add_argument("--input_bed", required=True, help="Path to the input BED file")
    parser.add_argument(
        "--sample_name",
        default=None,
        help="Sample name to parse in VCF. If not provided, the first sample will be used.",
    )
    parser.add_argument(
        "--output", default="output.csv", help="Output performance metrics csv file."
    )
    args = parser.parse_args(args)
    return args


def get_vcf_entry_state(entry, sample_name, vcf_type):
    """
    Returns BD and BK format fields from the first sample of a vcf entry
    """
    if sample_name == None:
        sample = entry.samples[0]
    else:
        sample = entry.samples[sample_name]

    bd = sample["BD"]
    try:
        bk = sample["BK"]
    except:
        bk = "not defined"

    ## TODO - incorporate bk
    if bd == "TP":
        if bk == "gm":
            if vcf_type == "comp":
                return "tp"
            elif vcf_type == "base":
                return "tpbase"
            else:
                print("Not defined TP")
                return None
        else:
            return None
    elif bd == "FP":
        return "fp"
    elif bd == "FN":
        return "fn"
    ## TODO - check this logic, not sure how to account for UNK
    # elif bd == "UNK":
    #     if vcf_type == "comp":
    #         return "fp"
    #     elif vcf_type == "base":
    #         return "fn"
    #     else:
    #         print("Not defined UNK")
    #         return None
    else:
        # print(f"Not defined BD - BD:{bd} BK:{bk}")
        return None


def calculate_var_states(region, base_vcf, comp_vcf, sample_name):
    """
    For variant entries within a region, sum the variant states
    returns a dictionary of {state: count}
    """
    chrom, start, end = region
    var_states = {"tpbase": 0, "tp": 0, "fn": 0, "fp": 0}
    for base_entry in base_vcf.fetch(chrom, int(start), int(end)):
        var_state = get_vcf_entry_state(base_entry, sample_name, vcf_type = "base")
        if var_state in var_states:
            var_states[var_state] += 1
    for comp_entry in comp_vcf.fetch(chrom, int(start), int(end)):
        var_state = get_vcf_entry_state(comp_entry, sample_name, vcf_type = "comp")
        if var_state in var_states:
            var_states[var_state] += 1
    return var_states


def strat_qfy_main(args):
    """
    Main
    """
    args = parse_args(args)

    df = pd.read_csv(args.input_bed, sep="\t")

    base_vcf = pysam.VariantFile(args.base_vcf)
    comp_vcf = pysam.VariantFile(args.comp_vcf)

    # Calculate the variant state counts for each region
    regions = df[["chrom", "start", "end"]].values.tolist()
    ## Note - potential place to introduce parallel processing
    region_var_states = [
        calculate_var_states(region, base_vcf, comp_vcf, args.sample_name)
        for region in regions
    ]

    # Calculate the performance metrics for each column
    out = {}
    ## Note - a second place for parallelization
    for column in df.columns[5:]:
        column_values = df[column].values.tolist()
        for i, value in enumerate(column_values):
            if value == 1:
                if column in out:
                    for k, v in region_var_states[i].items():
                        out[column][k] += v                    
                else:
                    out[column] = region_var_states[i]
        if column in out:
            out[column] = truvari.performance_metrics(**out[column])

    result = pd.DataFrame(out)

    result.to_csv(args.output, index=False)


if __name__ == "__main__":
    strat_qfy_main()
