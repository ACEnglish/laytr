import logging

import pysam
import truvari
import pandas as pd

def calculate_metrics_method1(var_states):
    """
    Loose Regional Comparison
    """
    TP = var_states['tp']
    FP = var_states['fp']
    FN = var_states['fn']
    # TN = var_states['tn']
    # NC = var_states['nc']
    # NCV = var_states['ncv']
    P = TP + FN
    N = TP + FP + FN
    # N = TN + FP

    metrics = {
        'sensitivity': TP / (TP + FN) if (TP + FN) != 0 else 0,
        'positive_predictive_value': TP / (TP + FP) if (TP + FP) != 0 else 0,
        'false_positive_rate': FP / N if N != 0 else 0,
       # 'no_call_rate': NC / (P + N) if (P + N) != 0 else 0,
       # 'no_call_variant_rate': NCV / P if P != 0 else 0
    }

    return metrics

def calculate_metrics_method2(var_states):
    """
    Allele Match Required
    """
    TP = var_states['tp']
    FP = var_states['fp']
    FN = var_states['fn']
    # NC = var_states['nc']
    # NCV = var_states['ncv']
    P = TP + FN
    N = TP + FP + FN

    metrics = {
        'sensitivity': TP / (TP + FN) if (TP + FN) != 0 else 0,
        'precision': TP / (TP + FP) if (TP + FP) != 0 else 0,
        'negative_predictive_value': FP / (FP + FN) if (FP + FN) != 0 else 0,
        'false_positive_rate': FP / N if N != 0 else 0,
        # 'no_call_rate': NC / (P + N) if (P + N) != 0 else 0,
        # 'no_call_variant_rate': NCV / P if P != 0 else 0
    }

    return metrics

def calculate_metrics_method3(var_states):
    """
    Genotype Match Required as well as Genotype Match and Local Phasing Required
    """
    TP = var_states['tp']
    FP = var_states['fp']
    FN = var_states['fn']
    # NC = var_states['nc']
    # NCV = var_states['ncv']
    P = TP + FN
    N = TP + FP + FN

    metrics = {
        'sensitivity': TP / (TP + FN) if (TP + FN) != 0 else 0,
        'precision': TP / (TP + FP) if (TP + FP) != 0 else 0,
        'negative_predictive_value': FP / (FP + FN) if (FP + FN) != 0 else 0,
        'false_positive_rate': FP / N if N != 0 else 0,
        # 'no_call_rate': NC / (P + N) if (P + N) != 0 else 0,
        # 'no_call_variant_rate': NCV / P if P != 0 else 0
    }

    return metrics

# TODO: Why are there multiple methods calculating very similar things?
# TODO: Try to figure out a way to "don't repeat yourself" (DRY principle)
metrics_functions = {
    'method1': calculate_metrics_method1,
    'method2': calculate_metrics_method2,
    'method3': calculate_metrics_method3,
}

def get_vcf_entry_state(entry):
    """
    Returns BD and BK format fields from the first sample of a vcf entry

    # TODO: by default, its fine to parse the first sample. Would parsing multi-sample ever be needed?
    # If so, this method should take a second parameter 'sample_name' which is available to the user via
    # the cli (which would be handled in `parse_args`)
    """
    bd = entry.samples[0]["BD"]
    bk = entry.samples[0]["BK"]
    if bd == 'TP':
        if bk == 'gm':
            return 'tp'
        else:
            return 'tpbase'
    elif bd == 'FP':
        return 'fp'
    elif bd == 'FN':
        return 'fn'
    else:
        return None

def calculate_var_states(region, truth_vcf, query_vcf):
    """
    For variant entries within a region, sum the variant states
    returns a dictionary of {state: count}
    """
    # TODO: this would have been `print`ed potentially millions of times
    # Put it in debug (if at all). I think it shouldn't be used
    # logging.debug(f"Calculating var states for region: {region} ")
    chrom, start, end = region
    var_states = {"tpbase": 0, "tp": 0, "fn": 0, "fp": 0}
    for truth_entry in truth_vcf.fetch(chrom, int(start), int(end)):
        for query_entry in query_vcf.fetch(chrom, int(start), int(end)):
            var_states[get_vcf_entry_state(truth_entry)] += 1
            var_states[get_vcf_entry_state(query_entry)] += 1
    return var_states

def calculate_metrics(var_states, metrics_func='truvari'):
    """
    Run a metrics calculation method
    """
    if metrics_func == 'truvari':
        return truvari.performance_metrics(**var_states)
    else:
        return metrics_functions[metrics_func](var_states)

def parse_args(args):
    """
    Given a list of command line arguments, parse them
    
    # TODO: This should use argparse and look like laytr.giabTR_report.parse_args
    """
    truth_vcf = "GRCh38_v4.2.1_HG2-verrkoV1.1-V0.7_dipcall-z2k_TRUTH.vcf.gz"  # replace with the path to the truth VCF file
    query_vcf = "GRCh38_v4.2.1_HG2-verrkoV1.1-V0.7_dipcall-z2k_QUERY.vcf.gz"  # replace with the path to the query VCF file
    # TODO: where does this file exist?
    # probably need to compress it and add it to the repository
    input_bed = "v3.1-GRCh38-subset-multiintersect.bed"  # replace with the path to the input TSV file

    # Let user choose metrics function, default is 'truvari'
    # TODO: interactive prompts get in the way of pipelining. Turn this into a parameter that has the `choices`
    metrics_func = input('Enter metrics function (truvari/method1/method2/method3/method4): ') or 'truvari'

    # Logging
    truvari.setup_logging()
    return truth_vcf, query_vcf, input_bed, metrics_func

def strat_qfy_main(args):
    """
    Main
    """
    truth_vcf, query_vcf, input_bed, metrics_func = parse_args(args)

    df = pd.read_csv(input_bed, sep='\t')

    truth_vcf = pysam.VariantFile(truth_vcf)
    query_vcf = pysam.VariantFile(query_vcf)

    # Calculate the variant state counts for each region
    regions = df[['chrom', 'start', 'end']].values.tolist()
    region_var_states = [calculate_var_states(region, truth_vcf, query_vcf) for region in regions]

    # Calculate the performance metrics for each column
    out = {}
    for column in df.columns[5:]:
        column_values = df[column].values.tolist()
        for i, value in enumerate(column_values):
            if value != 0:
                if column not in out:
                    out[column] = region_var_states[i]
                else:
                    for k, v in region_var_states[i].items():
                        out[column][k] += v
        if column in out:
            out[column] = calculate_metrics(out[column], metrics_func)

    result = pd.DataFrame(out)
    # TODO: output file should be an argument with a default of output.csv
    result.to_csv(f"output_{metrics_func}.csv", index=False)
