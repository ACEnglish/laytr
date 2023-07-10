import sys
import json
import laytr
import base64
import joblib
import logging
import truvari
import argparse
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from io import BytesIO
from matplotlib import colors, cm
from jinja2 import Environment, PackageLoader

def parse_args(args):
    """
    parse arguments
    """
    parser = argparse.ArgumentParser(prog="giabTRreport", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--regionsummary", type=str, required=True,
                        help="Truvari refine region summary ")
    parser.add_argument("-b", "--includebed", type=str, required=True,
                        help="GIAB TR benchmark")
    parser.add_argument("-t", "--trcatalog", type=str, required=True,
                        help="Adotto tandem repeat catalog")
    parser.add_argument("-s", "--som", type=str, required=True,
                        help="Laytr self-organizing map")
    parser.add_argument("-m", "--sommap", type=str, required=True,
                        help="Mapping of regions to neurons")
    parser.add_argument("-o", "--output", default="output.html",
                        help="Output html report")
    parser.add_argument("-d", "--dumpcsv", type=str,
                        help="Dump each report's CSV to filename prefix")
    args = parser.parse_args(args)
    return args

def inline_plot(plt):
    tmpfile = BytesIO()
    plt.figure.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
    
    return '<img src=\'data:image/png;base64,{}\'>'.format(encoded)

def make_hue_matrix(data, column):
    ret = np.empty((25, 25))

    for idx, row in data.iterrows():
        ret[idx] = row[column]
    return ret

def get_max_motif_length(row):
    data = json.loads(row['annos'])
    return max([anno['period'] for anno in data])

def region_stats(data):
    result = {}
    base = (data['out_tpbase'] != 0) | (data['out_fn'] != 0)
    baseP = int(base.sum())
    baseN = int((~base).sum())
    comp = (data['out_tp'] != 0) | (data['out_fp'] != 0)
    compP = int(comp.sum())
    compN = int((~comp).sum())

    state_cnts = data['state'].value_counts()
    result["TP"] = state_cnts["TP"] if "TP" in state_cnts else 0
    result["TN"] = state_cnts["TN"] if "TN" in state_cnts else 0

    fpfn = state_cnts["FN,FP"] if "FN,FP" in state_cnts else 0
    result["FP"] = fpfn + (state_cnts["FP"] if "FP" in state_cnts else 0)
    result["FN"] = fpfn + (state_cnts["FN"] if "FN" in state_cnts else 0)

    result["base P"] = baseP
    result["base N"] = baseN
    result["comp P"] = compP
    result["comp N"] = compN
    # precision
    result["PPV"] = result["TP"] / result["comp P"] if result["comp P"] != 0 else None
    # recall
    result["TPR"] = result["TP"] / result["base P"] if result["base P"] != 0 else None
    # specificity
    result["TNR"] = result["TN"] / result["base N"] if result["base N"] != 0 else None
    # negative predictive value
    result["NPV"] = result["TN"] / result["comp N"] if result["comp N"] != 0 else None
    # accuracy
    if result["base P"] + result["base N"] != 0:
        result["ACC"] = (result["TP"] + result["TN"]) / (result["base P"] + result["base N"])
    else:
        result["ACC"] = None
    if result["TPR"] is not None and result["TNR"] is not None:
        result["BA"] = (result["TPR"] + result["TNR"]) / 2
    else:
        result["BA"] = None

    if result["PPV"] and result["TPR"]:
        result["F1"] = 2 * ((result["PPV"] * result["TPR"]) / (result["PPV"] + result["TPR"]))
    else:
        result["F1"] = None
    return result

def get_expcon(x):
    if x['ad1'] == 0 and x['ad2'] == 0:
        return "REF"
    if x['ad1'] <= 0 and x['ad2'] <= 0:
        return "CON"
    if x['ad1'] >= 0 and x['ad2'] >= 0:
        return "EXP"
    return "MIX"

def get_maxadbin(x):
    sz = max(abs(x['ad1']), abs(x['ad2']))
    return truvari.get_sizebin(sz)

# Reports

def subset_report(data):
    """
    Main subsets table
    """
    full_report = region_stats(data)
    t1_report = region_stats(data[data["tier"] == "Tier1"])
    t2_report = region_stats(data[data["tier"] == "Tier2"])
    any_var_report = region_stats(data[data["var_flag"] != 0])
    gt5 = region_stats(data[(data['var_flag'] & 1 != 0) | (data['var_flag'] & 4 != 0)])

    sub_reports = pd.DataFrame([full_report, t1_report, t2_report, any_var_report, gt5],
                               index=pd.Index(["Full", "Tier1", "Tier2", "AnyVar", ">=5"], name="Subset"))
    return sub_reports

def make_recm_report(data):
    names = []
    rows = []
    for recm, sub in data.groupby('exp/con'):
        names.append(recm)
        rows.append(region_stats(sub))
    
    recm_report = pd.DataFrame(rows, index=pd.Index(names, name="ExpCon"))
    view = (recm_report.melt(ignore_index=False,
                            value_vars=["ACC", "BA", "F1"],
                            var_name="Metric",
                            value_name='Value')
                .reset_index())
    plt.figure()
    p1 = sb.barplot(data=view,
                   x="ExpCon",
                   y="Value",
                   hue="Metric",
                   order=["REF", "CON", "EXP", "MIX"])
    _ = p1.set(title="Performance by Allele Delta State", xlabel="Allele Deltas")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()

    view = (recm_report.melt(ignore_index=False,
                            value_vars=["PPV", "TPR", "TNR"],
                            var_name="Metric",
                            value_name='Value')
                .reset_index())
    plt.figure()
    p2 = sb.barplot(data=view,
                   x="ExpCon",
                   y="Value",
                   hue="Metric",
                   order=["REF", "CON", "EXP", "MIX"])
    _ = p2.set(title="Performance by Allele Delta State", xlabel="Allele Deltas")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()

    return recm_report, p1, p2

def make_perf_by_max_entropy(data):
    ba_by_ent_min = []
    ba_by_ent = []
    for i in range(0, 101, 2):
        i /= 100
        sub = region_stats(data[data['entropy'] <= i])
        ba_by_ent_min.append(i)
        ba_by_ent.append(sub)
    max_entropy = pd.DataFrame(ba_by_ent, index=pd.Index(ba_by_ent_min, name="Max Entropy"))

    plt.figure()
    p1 = sb.lineplot(data=max_entropy[["ACC", "BA", "F1"]].dropna())
    _ = p1.set(title="Performance by region entropy", ylabel="Metric")
    plt.tight_layout()

    plt.figure()
    p2 = sb.lineplot(data=max_entropy[["PPV", "TPR", "TNR"]].dropna())
    _ = p2.set(title="Performance by region entropy", ylabel="Metric")
    plt.tight_layout()
    
    return max_entropy, p1, p2

def make_perf_by_entropy_bin(data):
    names = []
    rows = []
    for idx, sub in data.groupby('entropy_bin'):
        names.append(idx)
        rows.append(region_stats(sub))
    return pd.DataFrame(rows, index=pd.Index(names, name="Entropy Bin"))

def make_perf_by_mxszbin(data):
    names = []
    rows = []
    for idx, sub in data.groupby('mx_szbin'):
        names.append(idx)
        rows.append(region_stats(sub))
    mxszbin_report = pd.DataFrame(rows, index=pd.Series(names).astype(truvari.SZBINTYPE)).loc[truvari.SZBINS]

    plt.figure()
    p1 = sb.lineplot(data=mxszbin_report[["ACC", "BA", "F1"]])
    plt.xticks(rotation=45, ha='right')
    _ = p1.set(title="Performance by max allele_delta length", xlabel="Size Bin", ylabel="Metric")
    plt.tight_layout()

    plt.figure()
    p2 = sb.lineplot(data=mxszbin_report[["PPV", "TPR", "TNR"]])
    plt.xticks(rotation=45, ha='right')
    _ = p2.set(title="Performance by max allele_delta length", xlabel="Size Bin", ylabel="Metric")
    plt.tight_layout()
    return mxszbin_report, p1, p2

def make_gene_report(data):
    no_gene = region_stats(data[data['gene_flag'] == 0])
    any_gene = region_stats(data[data['gene_flag'] != 0])
    prot_gene = region_stats(data[data['biotype'].str.contains("protein")])

    gene_report = pd.DataFrame([no_gene, any_gene, prot_gene],
                               index=pd.Index(["Intergenic", "Genic", "Protein Coding Gene"], name="Subset"))
    return gene_report

def make_interspersed_report(data):
    no_inter = region_stats(data[data['interspersed'] == '.'])
    has_inter = region_stats(data[data['interspersed'] != '.'])

    inter_report = pd.DataFrame([no_inter, has_inter], 
                               index=pd.Index(["No Interspersed", "Interspersed"], name="Subset"))
    return inter_report

def make_repeat_cpx_report(data):
    # Regions with a single isolated annotation are typically the easiest
    simple_overlap = (data["ovl_flag"] == 1) & (data["n_subregions"] == 1)

    # Regions with a single parent/nested annotations are next
    parent_overlap = (data["ovl_flag"].isin([2, 4, 6])) & (data["n_subregions"] == 1)

    # We'll consider the rest of the overlaps as being 'complex'
    cpx_overlap = ~simple_overlap & ~parent_overlap

    # Regions with multiple annotations or multiple sub-regions are typically difficult
    multi_anno = data["n_annos"] >= 2
    multi_subreg = data["n_subregions"] >= 2

    # The above subsets are boolean masks. We can use list comprehension to 
    # run each subset through `region_stats` and put them all together 
    results = [region_stats(data[_]) for _ in [simple_overlap,
                                               parent_overlap,
                                               cpx_overlap,
                                               multi_anno,
                                               multi_subreg]]
    cpx_report = pd.DataFrame(results, 
                               index=pd.Index(["Simple Overlap", "Parent Overlap", "Complex Overlap",
                                                "Multi Anno", "Multi SubReg"],
                               name="Subset"))
    return cpx_report


def make_max_motif_report(data):
    names = []
    rows = []
    for mbin, subset in data.groupby("max_motif_bin"):
        names.append(mbin)
        rows.append(region_stats(subset))
    motif_report = pd.DataFrame(rows, index=pd.Index(names, name="Motif Length Bin"))

    plt.figure()
    p = sb.lineplot(data=motif_report[["ACC", "BA", "F1"]])
    _ = p.set(title="Performance by maximum motif length", ylabel="Value")
    plt.tight_layout()
    return motif_report, p

def make_hompct_report(data):
    ba_by_hompct_min = []
    ba_by_hompct = []
    for i in range(0, 101, 5):
        sub = region_stats(data[data['hompct'] <= i])
        ba_by_hompct_min.append(i)
        ba_by_hompct.append(sub)
    max_hompct = pd.DataFrame(ba_by_hompct, index=pd.Index(ba_by_hompct_min, name="Max HomPct"))

    plt.figure()
    p1 = sb.lineplot(data=max_hompct[["ACC", "BA", "F1"]].dropna())
    _ = p1.set(title="Performance by region hompct", ylabel="Metric")
    plt.tight_layout()

    plt.figure()
    p2 = sb.lineplot(data=max_hompct[["PPV", "TPR", "TNR"]].dropna())
    _ = p2.set(title="Performance by region hompct", ylabel="Metric")
    plt.tight_layout()
    return max_hompct, p1, p2

def make_som_report(data, som):
    index = []
    rows = []
    for neuron, subset in data.groupby(["neuron_x", "neuron_y"]):
        index.append(neuron)
        rows.append(region_stats(subset))
    neuron_report = pd.DataFrame(rows, index=pd.Index(index))

    ba_matrix = make_hue_matrix(neuron_report, "BA")
    acc_matrix = make_hue_matrix(neuron_report, "ACC")
    f1_matrix = make_hue_matrix(neuron_report, "F1")
    
    plt.figure()
    norm = colors.Normalize(vmin=0.50, vmax=1.00)
    p1 = laytr.make_hex_plot(som, hue=ba_matrix, hue_label="Balanced Accuracy", color_map=cm.RdBu, color_norm=norm)
    plt.tight_layout()

    plt.figure()
    p2 = sb.histplot(data=neuron_report, x="BA")
    _ = p2.set(title="Balanced Accuracy per-Neuron", ylabel="Neuron Count")
    plt.tight_layout()

    plt.figure()
    norm = colors.Normalize(vmin=0.90, vmax=1.00)
    p3 = laytr.make_hex_plot(som, hue=acc_matrix, hue_label="Accuracy", color_map=cm.RdBu, color_norm=norm)
    plt.tight_layout()

    plt.figure()
    p4 = sb.histplot(data=neuron_report, x="ACC")
    _ = p4.set(title="Accuracy per-Neuron", ylabel="Neuron Count")
    plt.tight_layout()

    plt.figure()
    p5 = laytr.make_hex_plot(som, hue=f1_matrix, hue_label="F1", color_map=cm.RdBu)
    plt.tight_layout()

    plt.figure()
    p6 = sb.histplot(data=neuron_report, x="F1")
    _ = p6.set(title="F1 per-Neuron", ylabel="Neuron Count")
    plt.tight_layout()

    return p1, p2, p3, p4, p5, p6

def load_data(reg_fn, bed_fn, trc_fn, som_fn):
    """
    """
    regions = pd.read_csv(reg_fn, sep='\t')
    bed = pd.read_csv(bed_fn, sep='\t', names=["chrom", "start", "end", "tier", "replicates", "var_flag",
                                               "entropy", "ad1", "ad2"])   

    reg_header = ("chrom start end ovl_flag up_buff dn_buff hompct n_filtered "
                  "n_annos n_subregions mu_purity pct_annotated interspersed "
                  "patho codis gene_flag biotype annos").split(' ')
    catalog = (pd.read_csv(trc_fn, sep='\t', header=None, names=reg_header)
               .set_index(["chrom", "start", "end"]))

    bed.set_index(["chrom", "start", "end"], inplace=True)
    regions.set_index(["chrom", "start", "end"], inplace=True)
    
    som_map = joblib.load(som_fn)
    som_map_df = pd.DataFrame(som_map['map'],
                          index=pd.MultiIndex.from_tuples(som_map['index'], names=["chrom", "start", "end"]),
                          columns=['neuron_x', 'neuron_y'])

    data = regions.join(bed).join(catalog).join(som_map_df)

    data['exp/con'] = data.apply(get_expcon, axis=1)
    data['mx_szbin'] = data.apply(get_maxadbin, axis=1)
    data['entropy_bin'] = pd.cut(data['entropy'], bins=[0, 0.25, 0.50, 0.60, 0.70, 0.80, 0.90, 1], right=True)
    data['max_motif'] = data.apply(get_max_motif_length, axis=1)
    data['max_motif_bin'] = pd.cut(data["max_motif"], bins=[0, 5, 10, 20, 50, sys.maxsize], right=True,
                                   labels=["<5", "[5,10)", "[10, 20)", "[20, 50)", ">=50"])
    return data

def giabTRreport_main(args):
    """
    Args 
    """
    sb.set() # style the graphs
    args = parse_args(args)
    truvari.setup_logging()
    logging.info("Loading data")
    data = load_data(args.regionsummary, args.includebed, args.trcatalog, args.sommap)
    
    logging.info("Subsets report")
    sub_report = subset_report(data)

    logging.info("RECM report")
    recm_report, recm_p1, recm_p2 = make_recm_report(data)

    logging.info("Entropy report")
    entropy_report, entropy_p1, entropy_p2 = make_perf_by_max_entropy(data)

    logging.info("Sizebin report")
    sizebin_report, sizebin_p1, sizebin_p2 = make_perf_by_mxszbin(data)
    
    logging.info("Gene report")
    genetr_report = make_gene_report(data)
    
    logging.info("Interspersed report")
    inter_report = make_interspersed_report(data)
    
    logging.info("Complexity report")
    cpx_report = make_repeat_cpx_report(data)

    logging.info("Motif report")
    motif_report, motif_p1 = make_max_motif_report(data)
    
    logging.info("HomPct report")
    hompct_report, hompct_p1, hompct_p2 = make_hompct_report(data)

    logging.info("SOMs")
    som_1, som_2, som_3, som_4, som_5, som_6 = make_som_report(data, joblib.load(args.som))

    # Configure the Jinja environment to load templates from the current directory
    env = Environment(loader=PackageLoader('laytr', 'templates'))
    template = env.get_template('giabtr.html')

    # Render the template with the data and plots
    html_output = template.render(subset_report=sub_report.to_html(),
                                  subset_csv=sub_report.to_csv(),
                                  recm_report=recm_report.to_html(),
                                  recm_csv=recm_report.to_csv(),
                                  recm_p1=inline_plot(recm_p1),
                                  recm_p2=inline_plot(recm_p2),
                                  entropy_csv=entropy_report.to_csv(),
                                  entropy_p1=inline_plot(entropy_p1),
                                  entropy_p2=inline_plot(entropy_p2),
                                  sizebin_report=sizebin_report.to_html(),
                                  sizebin_csv=sizebin_report.to_csv(),
                                  sizebin_p1=inline_plot(sizebin_p1),
                                  sizebin_p2=inline_plot(sizebin_p2),
                                  genetr_report=genetr_report.to_html(),
                                  gene_csv=genetr_report.to_csv(),
                                  inter_report=inter_report.to_html(),
                                  inter_csv=inter_report.to_csv(),
                                  cpx_report=cpx_report.to_html(),
                                  repeat_csv=cpx_report.to_csv(),
                                  motif_report=motif_report.to_html(),
                                  motif_csv=motif_report.to_csv(),
                                  motif_p1=inline_plot(motif_p1),
                                  hompct_csv=hompct_report.to_csv(),
                                  hompct_p1=inline_plot(hompct_p1),
                                  hompct_p2=inline_plot(hompct_p2),
                                  som_1=inline_plot(som_1),
                                  som_2=inline_plot(som_2),
                                  som_3=inline_plot(som_3),
                                  som_4=inline_plot(som_4),
                                  som_5=inline_plot(som_5),
                                  som_6=inline_plot(som_6),
                                  )

    # Write the HTML output to a file
    with open(args.output, 'w') as f:
        f.write(html_output)
    
    if args.dumpcsv:
        logging.info("Dumping CSVs")
        for suffix, report in [("sub", sub_report), ("recm", recm_report), ("entropy", entropy_report),
                               ("sizebin", sizebin_report), ("gene", genetr_report), ("inter", inter_report),
                               ("cpx", cpx_report), ("motif", motif_report)]:
            report.to_csv(f"{args.dumpcsv}_{suffix}.csv")
        data.to_csv(f"{args.dumpcsv}_data.csv")
