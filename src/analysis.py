"""Load the processed data and run the core analyses."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (expression, samples, genes).

    expression: log2 values, genes as rows and samples as columns
    samples:    one row per sample with its group (Normal or Tumour)
    genes:      full gene name for each gene symbol
    """
    expression = pd.read_csv(DATA_DIR / "expression.csv.gz", index_col="gene")
    # Read patient as text so "001" and "1" stay different patients.
    samples = pd.read_csv(DATA_DIR / "samples.csv", dtype={"patient": str})
    samples["fuhrman_grade"] = samples["fuhrman_grade"].astype("Int64")
    genes = pd.read_csv(DATA_DIR / "genes.csv", index_col="gene")
    return expression, samples, genes


def differential_expression(expression: pd.DataFrame, samples: pd.DataFrame,
                            genes: pd.DataFrame) -> pd.DataFrame:
    """Compare tumour with normal for every gene.

    - log2 fold change: mean(tumour) minus mean(normal) on the log2 scale.
      +1 means twice as high in tumour, -1 means half as high.
    - p-value: Welch's t-test, which does not assume equal variance in the two groups.
    - FDR: Benjamini-Hochberg adjusted p-value, which controls false discoveries
      when testing thousands of genes at once.
    """
    tumour = expression[samples.loc[samples["group"] == "Tumour", "sample"]]
    normal = expression[samples.loc[samples["group"] == "Normal", "sample"]]

    _, p_values = stats.ttest_ind(tumour, normal, axis=1, equal_var=False)
    results = pd.DataFrame({
        "log2_fold_change": tumour.mean(axis=1) - normal.mean(axis=1),
        "mean_normal": normal.mean(axis=1),
        "mean_tumour": tumour.mean(axis=1),
        "p_value": p_values,
    }, index=expression.index)
    results["fdr"] = stats.false_discovery_control(results["p_value"], method="bh")
    results["neg_log10_fdr"] = -np.log10(results["fdr"])
    results = results.join(genes["name"])
    return results.sort_values("p_value")


def label_significance(results: pd.DataFrame, fc_cutoff: float, fdr_cutoff: float) -> pd.Series:
    """Tag each gene as Up in tumour, Down in tumour, or Not significant."""
    up = (results["log2_fold_change"] >= fc_cutoff) & (results["fdr"] < fdr_cutoff)
    down = (results["log2_fold_change"] <= -fc_cutoff) & (results["fdr"] < fdr_cutoff)
    return pd.Series(np.select([up, down], ["Up in tumour", "Down in tumour"],
                               "Not significant"), index=results.index)


def run_pca(expression: pd.DataFrame, n_genes: int = 1000):
    """PCA on the most variable genes. Returns (coordinates, explained variance ratio)."""
    top = expression.loc[expression.var(axis=1).nlargest(n_genes).index]
    centred = top.T - top.T.mean()
    pca = PCA(n_components=2)
    coords = pd.DataFrame(pca.fit_transform(centred), index=top.columns, columns=["PC1", "PC2"])
    return coords, pca.explained_variance_ratio_


def top_genes(results: pd.DataFrame, n: int = 25) -> list:
    """The n most significant genes, split evenly between up and down in tumour."""
    sig = results[results["fdr"] < 0.05]
    up = sig[sig["log2_fold_change"] > 0].nlargest(n // 2 + n % 2, "log2_fold_change")
    down = sig[sig["log2_fold_change"] < 0].nsmallest(n // 2, "log2_fold_change")
    return list(up.index) + list(down.index)
