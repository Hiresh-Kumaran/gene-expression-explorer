"""Plotly charts for the app. Each function takes data and returns a figure."""

import math

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

GROUP_COLORS = {"Normal": "#4C9F70", "Tumour": "#D1495B"}
SIG_COLORS = {"Up in tumour": "#D1495B", "Down in tumour": "#3A6EA5", "Not significant": "#C8C8C8"}


def sample_distributions(expression: pd.DataFrame, samples: pd.DataFrame) -> go.Figure:
    """Box plot of every sample's expression values, a basic quality check."""
    long = expression.melt(var_name="sample", value_name="log2 expression")
    long = long.merge(samples[["sample", "group"]], on="sample")
    fig = px.box(long, x="sample", y="log2 expression", color="group",
                 color_discrete_map=GROUP_COLORS, points=False,
                 category_orders={"sample": list(samples["sample"])},
                 title="Expression distribution per sample")
    fig.update_layout(xaxis_title="", legend_title="")
    return fig


def pca_scatter(coords: pd.DataFrame, variance, samples: pd.DataFrame) -> go.Figure:
    data = coords.join(samples.set_index("sample")[["group", "patient"]])
    data = data.rename_axis("sample").reset_index()
    fig = px.scatter(data, x="PC1", y="PC2", color="group", text="sample",
                     color_discrete_map=GROUP_COLORS, hover_data=["patient"],
                     title="PCA of the 1,000 most variable genes")
    fig.update_traces(marker=dict(size=13), textposition="top center")
    fig.update_layout(xaxis_title=f"PC1 ({variance[0]:.0%} of variance)",
                      yaxis_title=f"PC2 ({variance[1]:.0%} of variance)", legend_title="")
    return fig


def volcano(results: pd.DataFrame, labels: pd.Series, fc_cutoff: float,
            fdr_cutoff: float, annotate: int = 5) -> go.Figure:
    data = results.assign(significance=labels).reset_index()
    fig = px.scatter(data, x="log2_fold_change", y="neg_log10_fdr", color="significance",
                     color_discrete_map=SIG_COLORS, hover_name="gene",
                     hover_data={"name": True, "fdr": ":.2e", "significance": False},
                     category_orders={"significance": list(SIG_COLORS)},
                     title="Volcano plot: tumour vs normal")
    fig.update_traces(marker=dict(size=5, opacity=0.7))

    # Dashed lines show the current cut-offs.
    for x in (-fc_cutoff, fc_cutoff):
        fig.add_vline(x=x, line_dash="dash", line_color="grey")
    fig.add_hline(y=-math.log10(fdr_cutoff), line_dash="dash", line_color="grey")

    # Label the strongest genes on each side (large change and strong evidence).
    data["score"] = data["log2_fold_change"].abs() * data["neg_log10_fdr"]
    to_label = pd.concat([data[data["significance"] == side].nlargest(annotate, "score")
                          for side in ("Up in tumour", "Down in tumour")])
    for _, row in to_label.iterrows():
        fig.add_annotation(x=row["log2_fold_change"], y=row["neg_log10_fdr"], text=row["gene"],
                           showarrow=False, yshift=10, font=dict(size=11))

    fig.update_layout(xaxis_title="log2 fold change (tumour vs normal)",
                      yaxis_title="-log10 FDR", legend_title="", height=600)
    return fig


def heatmap(expression: pd.DataFrame, samples: pd.DataFrame, genes: list) -> go.Figure:
    """Z-scored expression of selected genes, samples ordered normal then tumour."""
    order = list(samples.sort_values("group")["sample"])
    data = expression.loc[genes, order]
    z = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)
    fig = px.imshow(z, x=order, y=genes, color_continuous_scale="RdBu_r",
                    zmin=-2.5, zmax=2.5, aspect="auto",
                    title="Top differentially expressed genes")
    fig.update_layout(coloraxis_colorbar_title="z-score", height=max(450, 22 * len(genes)),
                      xaxis_title="Normal samples  |  Tumour samples", yaxis_title="")
    return fig


def gene_boxplot(expression: pd.DataFrame, samples: pd.DataFrame, gene: str) -> go.Figure:
    """Expression of one gene in each sample, with lines joining samples from the same patient."""
    data = samples.assign(value=expression.loc[gene, samples["sample"]].values)
    fig = px.box(data, x="group", y="value", color="group", points="all",
                 color_discrete_map=GROUP_COLORS, hover_data=["sample", "patient"],
                 category_orders={"group": ["Normal", "Tumour"]},
                 title=f"{gene} expression")
    fig.update_traces(pointpos=0, jitter=0, selector=dict(type="box"))  # dots sit on the boxes
    paired = data.groupby("patient").filter(lambda d: d["group"].nunique() == 2)
    for _, pair in paired.groupby("patient"):
        pair = pair.sort_values("group")
        fig.add_trace(go.Scatter(x=pair["group"], y=pair["value"], mode="lines",
                                 line=dict(color="lightgrey", width=1),
                                 showlegend=False, hoverinfo="skip"))
    fig.update_layout(xaxis_title="", yaxis_title="log2 expression", showlegend=False)
    return fig
