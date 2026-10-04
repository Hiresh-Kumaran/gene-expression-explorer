"""Gene Expression Explorer: tumour vs normal kidney tissue (GEO GSE781)."""

import streamlit as st

from src import plots
from src.analysis import (differential_expression, label_significance, load_data,
                          run_pca, top_genes)

st.set_page_config(page_title="Gene Expression Explorer", page_icon="🧬", layout="wide")

# Short notes on well-known genes in clear cell renal cell carcinoma (ccRCC).
GENE_NOTES = {
    "CA9": "Carbonic anhydrase IX. Switched on by HIF when cells sense low oxygen, and helps "
           "tumour cells cope with acidity. Pathologists use CA9 staining to help identify ccRCC.",
    "NDUFA4L2": "A HIF target that dials down mitochondrial oxygen use. One of the most strongly "
                "increased genes in ccRCC.",
    "ANGPTL4": "Angiopoietin-like 4, a HIF target involved in blood vessel growth and fat metabolism.",
    "VEGFA": "The main signal for new blood vessel growth. ccRCC tumours are rich in blood vessels, "
             "and drugs that block VEGF signalling (such as sunitinib) are standard treatments.",
    "EGLN3": "An oxygen-sensing enzyme that is itself switched on by HIF, part of a feedback loop.",
    "HK2": "Hexokinase 2, the first step of glucose breakdown. Higher levels reflect the shift to "
           "glycolysis that many tumours make (the Warburg effect).",
    "UMOD": "Uromodulin, made only in the kidney's loop of Henle and the most abundant protein in "
            "normal urine. Its loss shows the tumour no longer behaves like kidney tissue.",
    "ALDOB": "Aldolase B, a sugar metabolism enzyme of the proximal tubule. Low levels in ccRCC "
             "have been linked to worse outcomes.",
    "KNG1": "Kininogen 1, a blood pressure and clotting protein normally made by the kidney and liver.",
    "CALB1": "Calbindin, which handles calcium in the distal tubule. A marker of normal kidney tissue.",
    "SLC12A1": "The sodium-potassium-chloride transporter of the loop of Henle, and the target of "
               "loop diuretics such as furosemide.",
}


@st.cache_data
def get_data():
    expression, samples, genes = load_data()
    results = differential_expression(expression, samples, genes)
    coords, variance = run_pca(expression)
    return expression, samples, genes, results, coords, variance


expression, samples, genes, results, coords, variance = get_data()

st.title("🧬 Gene Expression Explorer")
st.markdown("Which genes are switched on or off in **kidney cancer**? Explore microarray data "
            "comparing clear cell renal cell carcinoma with normal kidney tissue.")

tab_overview, tab_qc, tab_de, tab_heatmap, tab_gene = st.tabs(
    ["Overview", "Quality and PCA", "Differential expression", "Heatmap", "Gene lookup"])

# ---------- Overview ----------
with tab_overview:
    n_tumour = int((samples["group"] == "Tumour").sum())
    paired = samples.groupby("patient")["group"].nunique().eq(2).sum()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", len(samples))
    c2.metric("Tumour / normal", f"{n_tumour} / {len(samples) - n_tumour}")
    c3.metric("Matched patient pairs", int(paired))
    c4.metric("Genes analysed", f"{len(expression):,}")

    st.markdown(
        "**The study.** Lenburg et al. (2003) measured gene expression in clear cell renal cell "
        "carcinoma (ccRCC), the most common type of kidney cancer, and in normal kidney tissue "
        "taken next to the tumour. Most patients gave both a tumour and a normal sample.\n\n"
        "**The biology.** Most ccRCC tumours lose a gene called *VHL*. Without it, a protein "
        "called HIF builds up and switches on the genes a cell normally uses when oxygen is "
        "low, even when oxygen is plentiful. This drives blood vessel growth and changes how "
        "the tumour uses energy. You can see this signature in the results.\n\n"
        "**The method.** Each gene is compared between tumour and normal samples with a "
        "t-test, and p-values are adjusted for testing thousands of genes at once (false "
        "discovery rate). Genes are called significant when they change at least two-fold "
        "with an FDR below 5%, and both cut-offs can be changed in the app.")

    st.subheader("Samples")
    st.dataframe(samples.rename(columns={
        "sample": "Sample", "gsm": "GEO ID", "patient": "Patient", "group": "Group",
        "age": "Age", "sex": "Sex", "fuhrman_grade": "Tumour grade (Fuhrman)"}),
        hide_index=True, width="stretch")

# ---------- Quality and PCA ----------
with tab_qc:
    st.plotly_chart(plots.sample_distributions(expression, samples), width="stretch")
    st.caption("All samples have similar distributions, so no sample stands out as a "
               "technical outlier. The data was already scaled to a common target by the "
               "original authors.")
    st.plotly_chart(plots.pca_scatter(coords, variance, samples), width="stretch")
    st.caption(f"PCA squeezes 1,000 genes into two dimensions. The first component "
               f"({variance[0]:.0%} of the variation) separates every tumour from every normal "
               f"sample, so the tumour vs normal difference is the strongest signal in the data.")

# ---------- Differential expression ----------
with tab_de:
    c1, c2 = st.columns(2)
    fc_cutoff = c1.slider("Minimum log2 fold change", 0.5, 4.0, 1.0, 0.5,
                          help="1 means at least two-fold, 2 means at least four-fold.")
    fdr_cutoff = c2.select_slider("Maximum FDR", [0.001, 0.01, 0.05, 0.1], value=0.05)
    labels = label_significance(results, fc_cutoff, fdr_cutoff)

    counts = labels.value_counts()
    m1, m2, m3 = st.columns(3)
    m1.metric("Up in tumour", int(counts.get("Up in tumour", 0)))
    m2.metric("Down in tumour", int(counts.get("Down in tumour", 0)))
    m3.metric("Not significant", int(counts.get("Not significant", 0)))

    st.plotly_chart(plots.volcano(results, labels, fc_cutoff, fdr_cutoff), width="stretch")
    st.caption("Each dot is a gene. Further right means higher in tumour, further left means "
               "lower, and higher up means stronger statistical evidence.")

    table = results.assign(significance=labels)
    table = table[table["significance"] != "Not significant"]
    st.dataframe(
        table[["name", "log2_fold_change", "mean_normal", "mean_tumour", "fdr", "significance"]]
        .rename(columns={"name": "Gene name", "log2_fold_change": "log2 FC",
                         "mean_normal": "Mean normal", "mean_tumour": "Mean tumour",
                         "fdr": "FDR", "significance": "Direction"})
        .style.format({"log2 FC": "{:.2f}", "Mean normal": "{:.2f}",
                       "Mean tumour": "{:.2f}", "FDR": "{:.2e}"}),
        width="stretch", height=400)
    st.download_button("Download significant genes (CSV)", table.to_csv().encode(),
                       "significant_genes.csv", "text/csv")

# ---------- Heatmap ----------
with tab_heatmap:
    n = st.slider("Number of genes", 10, 60, 30, 10)
    st.plotly_chart(plots.heatmap(expression, samples, top_genes(results, n)), width="stretch")
    st.caption("The genes with the largest significant changes, half up and half down. Each "
               "row is scaled so red means above that gene's average and blue means below. "
               "The clean block pattern shows these genes separate tumour from normal tissue.")

# ---------- Gene lookup ----------
with tab_gene:
    gene_list = list(expression.index)
    gene = st.selectbox("Choose a gene (type to search)", gene_list,
                        index=gene_list.index("CA9"))
    row = results.loc[gene]
    left, right = st.columns([3, 2])
    left.plotly_chart(plots.gene_boxplot(expression, samples, gene), width="stretch")
    right.subheader(gene)
    right.markdown(f"*{row['name']}*")
    right.metric("log2 fold change (tumour vs normal)", f"{row['log2_fold_change']:+.2f}",
                 help="+1 means twice as high in tumour, -1 means half as high.")
    right.metric("FDR", f"{row['fdr']:.1e}")
    if gene in GENE_NOTES:
        right.info(GENE_NOTES[gene])
    right.caption("Grey lines join tumour and normal samples from the same patient. "
                  "Try: " + ", ".join(GENE_NOTES))

st.divider()
st.caption("Data: Lenburg ME et al. (2003) BMC Cancer 3:31, GEO accession GSE781 "
           "(Affymetrix HG-U133A). For research and education only.")
