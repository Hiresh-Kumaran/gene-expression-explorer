"""Build the processed dataset used by the app from GEO series GSE781.

Run once from the project root:
    pip install GEOparse
    python scripts/prepare_data.py

Outputs (already included in data/, so the app does not need to run this):
    data/expression.csv.gz   log2 expression, one row per gene, one column per sample
    data/samples.csv         sample id, patient, group, age, sex, tumour grade
    data/genes.csv           gene symbol and full gene name
"""

import argparse
import re

import GEOparse
import numpy as np
import pandas as pd

PLATFORM = "GPL96"  # Affymetrix HG-U133A, the better annotated of the two arrays in GSE781
MIN_PRESENT = 3     # keep probes detected ("Present") in at least this many samples


def parse_description(text: str) -> dict:
    """Pull age, sex and Fuhrman grade out of the free-text sample description."""
    age = re.search(r"Age = (\d+)", text)
    sex = re.search(r"Gender = (\w+)", text)
    grade = re.search(r"Fuhrman Grade = (\d)", text)
    return {
        "age": int(age.group(1)) if age else None,
        "sex": sex.group(1) if sex else None,
        "fuhrman_grade": int(grade.group(1)) if grade else None,
    }


def build(gse) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    samples, values, calls = [], {}, {}
    for gsm_id, gsm in gse.gsms.items():
        if gsm.metadata["platform_id"][0] != PLATFORM:
            continue
        source = gsm.metadata["source_name_ch1"][0]
        title = gsm.metadata["title"][0]
        # Group comes from the tissue source, not the title: three normal samples
        # (N2, N3, N4) are mislabelled as carcinoma in their GEO titles.
        group = "Normal" if "normal tissue" in source.lower() else "Tumour"
        sample_name = title.split()[0]  # e.g. "C035" or "N035"
        samples.append({
            "sample": sample_name,
            "gsm": gsm_id,
            "patient": sample_name[1:],  # C035 and N035 come from the same patient
            "group": group,
            **parse_description(" ".join(gsm.metadata.get("description", []))),
        })
        if group == "Normal":  # for normal tissue the grade describes the nearby tumour
            samples[-1]["fuhrman_grade"] = None
        table = gsm.table.set_index("ID_REF")
        values[sample_name] = table["VALUE"]
        calls[sample_name] = table["ABS_CALL"]

    samples = pd.DataFrame(samples).sort_values(["group", "sample"]).reset_index(drop=True)
    values = pd.DataFrame(values)[samples["sample"]]
    calls = pd.DataFrame(calls)[samples["sample"]]

    # 1. Remove probes that are not reliably detected.
    keep = (calls == "P").sum(axis=1) >= MIN_PRESENT
    log2 = np.log2(values[keep].clip(lower=1))

    # 2. Map probes to gene symbols, dropping control probes and unannotated ones.
    annot = gse.gpls[PLATFORM].table.set_index("ID")[["Gene Symbol", "Gene Title"]].dropna()
    annot = annot[~annot.index.str.startswith("AFFX")]
    annot["Gene Symbol"] = annot["Gene Symbol"].str.split(" /// ").str[0]
    log2 = log2.join(annot["Gene Symbol"], how="inner")

    # 3. Several probes can measure one gene: keep the most highly expressed probe.
    log2["mean"] = log2[samples["sample"]].mean(axis=1)
    log2 = log2.sort_values("mean", ascending=False).drop_duplicates("Gene Symbol")
    expression = log2.set_index("Gene Symbol")[samples["sample"]].sort_index().round(3)
    expression.index.name = "gene"

    genes = (annot.drop_duplicates("Gene Symbol").set_index("Gene Symbol")
             .loc[expression.index, ["Gene Title"]].rename(columns={"Gene Title": "name"}))
    genes.index.name = "gene"
    return expression, samples, genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft", help="path to a downloaded GSE781_family.soft.gz (optional)")
    args = parser.parse_args()

    if args.soft:
        gse = GEOparse.get_GEO(filepath=args.soft, silent=True)
    else:
        gse = GEOparse.get_GEO(geo="GSE781", destdir=".", silent=True)

    expression, samples, genes = build(gse)
    expression.to_csv("data/expression.csv.gz")
    samples.to_csv("data/samples.csv", index=False)
    genes.to_csv("data/genes.csv")
    print(f"Saved {expression.shape[0]} genes x {expression.shape[1]} samples")


if __name__ == "__main__":
    main()
