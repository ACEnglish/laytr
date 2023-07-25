import pandas as pd
import pysam
import truvari

# Method #1: Loose Regional Comparison
def calculate_metrics_method1(var_states):
    TP = var_states['tp']
    FP = var_states['fp']
    FN = var_states['fn']
#    TN = var_states['tn']
    # NC = var_states['nc']
    # NCV = var_states['ncv']
    P = TP + FN
    N = TP + FP + FN
 #   N = TN + FP

    metrics = {
        'sensitivity': TP / (TP + FN) if (TP + FN) != 0 else 0,
        'positive_predictive_value': TP / (TP + FP) if (TP + FP) != 0 else 0,
        'false_positive_rate': FP / N if N != 0 else 0,
       # 'no_call_rate': NC / (P + N) if (P + N) != 0 else 0,
       # 'no_call_variant_rate': NCV / P if P != 0 else 0
    }

    return metrics

# Method #2: Allele Match Required
def calculate_metrics_method2(var_states):
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

# Method #3: Genotype Match Required
def calculate_metrics_method3(var_states):
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

# Method #4: Genotype Match and Local Phasing Required
# Same as method #3
calculate_metrics_method4 = calculate_metrics_method3


metrics_functions = {
    'method1': calculate_metrics_method1,
    'method2': calculate_metrics_method2,
    'method3': calculate_metrics_method3,
    'method4': calculate_metrics_method4,
}

def get_vcf_entry_state(entry):
    bd = entry.samples.values()[0]["BD"]
    bk = entry.samples.values()[0]["BK"]
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
    print(f"Calculating var states for region: {region} ")
    chrom, start, end = region
    var_states = {"tpbase": 0, "tp": 0, "fn": 0, "fp": 0}
    for truth_entry in truth_vcf.fetch(chrom, int(start), int(end)):
        for query_entry in query_vcf.fetch(chrom, int(start), int(end)):
            var_states[get_vcf_entry_state(truth_entry)] += 1
            var_states[get_vcf_entry_state(query_entry)] += 1
    return var_states

def calculate_metrics(var_states, metrics_func='truvari'):
    if metrics_func == 'truvari':
        return truvari.performance_metrics(**var_states)
    else:
        return metrics_functions[metrics_func](var_states)

if __name__ == '__main__':
    truth_vcf = "GRCh38_v4.2.1_HG2-verrkoV1.1-V0.7_dipcall-z2k_TRUTH.vcf.gz"  # replace with the path to the truth VCF file
    query_vcf = "GRCh38_v4.2.1_HG2-verrkoV1.1-V0.7_dipcall-z2k_QUERY.vcf.gz"  # replace with the path to the query VCF file
    # truth_vcf = "ga4gh_truth.vcf.gz"  # replace with the path to the truth VCF file
    # query_vcf = "ga4gh_query.vcf.gz"  # replace with the path to the query VCF file
    input_bed = "v3.1-GRCh38-subset-multiintersect.bed"  # replace with the path to the input TSV file

    # Let user choose metrics function, default is 'truvari'
    metrics_func = input('Enter metrics function (truvari/method1/method2/method3/method4): ') or 'truvari'

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
    result.to_csv(f"output_{metrics_func}.csv", index=False)
