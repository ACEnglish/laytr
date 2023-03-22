import logging
import joblib
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


pd.set_option("display.max_rows", None)
sb.set()
def get_scores(df, state="state", features=["cat2","svtype","szbin"], min_obs=10):
    """
    Actively filter 0 in either - has to be something
    """
    tot_counts = []
    for feat in features:
        cnt = df[feat].value_counts()
        cnt.name = "count"
        cnt = pd.DataFrame(cnt)
        cnt['feat'] = feat
        tot_counts.append(cnt)
    tot_counts = pd.concat(tot_counts)
    tot_counts = tot_counts.reset_index().set_index(['feat', 'index'])
    
    cont_obs = pd.crosstab(df[state], [df[i] for i in features])
    view = pd.DataFrame(cont_obs.sum(axis=0))
    cont_obs = cont_obs.T[view[0] >= min_obs].T
    
    chi, pval, dof, exp = chi2_contingency(cont_obs)
    score = (np.sign(cont_obs.values[1] - exp[1])) * \
            (cont_obs.values[1] - exp[1] )**2 / exp[1]
    
    df_rank = pd.DataFrame(list(zip(list(cont_obs.columns), 
                                     list(cont_obs.values[1]), list(exp[1]), list(score))), 
                            columns=["Values", "observed", "expected", "score"])

    df_rank.sort_values(["score"], ascending=True, inplace=True)
    df_rank.reset_index(inplace=True, drop=True)
    df_rank.index.names = ['rank']
    df_rank[features] = pd.DataFrame(df_rank["Values"].tolist(), index=df_rank.index)
    df_rank = df_rank.drop("Values", axis=1)
    return df_rank, tot_counts

def permutation_test(a_values, b_values, n_samps=10000, tailed="two"):
    """
    a_values - b_values mean difference from permutations' mean differences
    tail can be "two", "left", or "right"
    """
    delta = a_values.mean() - b_values.mean()
    if tailed == "two":
        delta = abs(delta)
    logging.debug("Case %.2f (%d)", a_values.mean(), len(a_values))
    logging.debug("Control %.2f (%d)", b_values.mean(), len(b_values))
    logging.debug("delta %.2f", delta)
    def run_t(pooled, sizeA, sizeB):
        np.random.shuffle(pooled)
        starA = pooled[:sizeA]
        starB = pooled[-sizeB:]
        return starA.mean() - starB.mean()
    pooled = np.hstack([a_values, b_values])
    estimates = np.array([run_t(pooled, a_values.size, b_values.size) for _ in range(n_samps)])
    if tailed == "left":
        diffCount = len(np.where(estimates <= delta)[0])
    elif tailed == "right":
        diffCount = len(np.where(estimates >= delta)[0])
    else:
        estimates = np.abs(estimates)
        diffCount = len(np.where(estimates >= delta)[0])
    pct_extreme = float(diffCount) / float(n_samps)
    one, fif, nine = np.quantile(estimates, [0.01, 0.5, 0.99])
    logging.debug("Quantiles: 1% {:.10f} - 50% {:.10f} - 99% {:.10f}".format(one, fif, nine))
    logging.debug("Pct at least as extreme %.2f", pct_extreme)
    return pct_extreme, delta, one, fif, nine

def surb_score(m_data, state="state", features=["cat2", "svtype", "szbin"], tail="left", min_obs=10, alpha=0.01):
    """
    create the surb score data-frame that's plottable
    """
    rank, counts = get_scores(m_data, state, features, min_obs=10)
    rows = []
    # Also need the pvalue distributions (normalized?)
    for feat in features:
        for val in rank[feat].unique():
            sub = rank[feat] == val
            case = np.array(rank[sub].index)
            control = np.array(rank[~sub].index)
            if len(case) < min_obs or len(control) < min_obs:
                continue
            pval, delta, one, fif, nin = permutation_test(case, control, tailed=tail)
            acc = m_data[m_data[feat] == val][state].mean()
            rows.append([feat, val, pval, acc, delta, one, fif, nin, counts.loc[feat, val]["count"]])

    plt_rank = pd.DataFrame(rows, columns=["feature", "value", "pval", "acc", "delta", "q1", "q50", "q99", "obs"])
    reject, adj_pval, x, y = multipletests(pvals=plt_rank['pval'], alpha=alpha, method="fdr_bh")
    adj = pd.concat([plt_rank,
                     pd.Series(reject, name="reject"),
                     pd.Series(adj_pval, name="adj_pval")], axis=1)
    return adj

def make_plotly(m_dat):
    """
    make interactive plotly charts from a surb_score data frame
    """
    custom_text = m_dat.apply(lambda x: "{}: {:.1f}% {:,d}".format(x["value"], x["acc"] * 100, x['obs']), axis=1)
    plts = [go.Scatter(x=m_dat.feature, 
                       y=m_dat.delta,
                       mode="markers", 
                       hovertemplate="%{text}",
                       text=custom_text,
                       name="",
                       marker=dict(size=10,
                                   line=dict(
                                       color='MediumPurple',
                                       width=2
                                   )
                        )
                    )
            ]
    dat = []
    for feat in m_dat.feature.unique():
        a = m_dat[m_dat["feature"] == feat]
        for pos, j in a.iterrows():
            dat.append([feat, j.q1])
            dat.append([feat, j.q50])
            dat.append([feat, j.q99])
            
        #plts.append(go.Scatter(x=(feat, feat), y=(a.q1.mean(), a.q99.mean()), 
        #                    hoverinfo='skip', name="", 
        #                    mode="lines",
        #                    marker=dict(size=1,color="black")))
    dat = pd.DataFrame(dat, columns=['feature', 'value'])
    plts.append(go.Box(x=dat['feature'], y=dat['value'], 
                       hoverinfo='skip', 
                       line_color='MediumPurple',
                       fillcolor='rgba(26,150,65,0)'))
    fig = go.FigureWidget(data=plts)
    fig.update_layout(showlegend=False) 
    return fig


def test():
    data = joblib.load('data.jl')
    features = ["svtype", "szbin", "cat", "remap", "cls"]
    class_var = "state"
    # top-down

    print("overall accuracy: %.0f%%" % (data["state"].mean() * 100).round(0))
    summary = (data.groupby("svtype")['state'].mean() * 100).round(0)
    print("accuracy by svtype: DEL(%.0f%%) INS(%.0f%%)" % (summary.loc['DEL'], summary.loc['INS']))

    tot = 0
    for i in features:
        tot += len(data[i].unique())
    print("%d features - %d values - ~%d groups" % (len(features), tot, len(data.groupby(features))))

    view = (data.groupby(['svtype', 'szbin'])['state'].mean() * 100).round(0)
    p = sb.barplot(data=view.reset_index(), x="szbin", y="state", hue="svtype")
    plt.xticks(rotation=45)
    hide = p.set(title="Accuracy by szbin/svtype")
    plt.show()


    # bottom-down

    rows = []
    for i in features:
        view = data.groupby([i])
        view = pd.DataFrame(view["state"].mean()).join(pd.DataFrame(view["state"].size()), lsuffix="_mu", rsuffix="_cnt")
        view['feat'] = i
        rows.append(view)

    rows = pd.concat(rows)
    rows.sort_values(['state_mu'])

    p = sb.boxplot(data=rows.reset_index(), x="feat", y="state_mu")
    hide = p.set(title="Per-feature accuracy")
    plt.show()

    mu = data.groupby(features)["state"].mean()
    cnt = data.groupby(features)["state"].size()
    view = pd.DataFrame(mu).join(pd.DataFrame(cnt), lsuffix="_mu", rsuffix="_cnt").dropna().sort_values('state_cnt')
    view['state_mu'] = (view['state_mu'] * 100).astype(int)
    hold = view.xs("[400,600)", level="szbin")
    hold.columns = ["mu", "cnt"]
    hold[hold["cnt"] >= 50]

    # SurbScore
    plt_rank = surb_score(data, "state", features)
    plt_rank[plt_rank['reject']][['feature', 'value', 'acc', 'obs', 'adj_pval']]
    fig = make_plotly(plt_rank)
    fig.show()
    #fig.write_html("biograph_figure.html")
