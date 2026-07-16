from __future__ import annotations

import csv
import gzip
import tarfile
from io import TextIOWrapper
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import fisher_exact, mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
SOURCE_TAR = ROOT / "sources_dataset4" / "GSE183904_RAW.tar"
SAMPLE_MANIFEST = ROOT / "config" / "sample_manifest.tsv"
EPITHELIAL_MANIFEST = ROOT / "outputs" / "tables" / "epithelial_cell_manifest.tsv"
FINAL_METADATA = ROOT / "outputs" / "final" / "tables" / "final_cell_metadata.tsv.gz"
OUT_DIR = ROOT / "outputs" / "final" / "kcnq1_primary_normal_vs_tumor"
GENE = "KCNQ1"

COLORS = {
    "Primary normal epithelial": "#8ecae6",
    "Tumor intestinal-like": "#2b9348",
    "Tumor diffuse/EMT-like": "#ff5714",
    "Tumor intestinal patient": "#2b9348",
    "Tumor diffuse patient": "#ff5714",
}


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out.tolist()


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_manifest = pd.read_csv(SAMPLE_MANIFEST, sep="\t")
    epithelial = pd.read_csv(EPITHELIAL_MANIFEST, sep="\t")
    final = pd.read_csv(FINAL_METADATA, sep="\t")
    return sample_manifest, epithelial, final


def build_cell_sets(epithelial: pd.DataFrame, final: pd.DataFrame) -> dict[str, pd.DataFrame]:
    base_cols = [
        "cell_id",
        "barcode",
        "sample_id",
        "gsm",
        "patient",
        "tissue",
        "lauren",
        "total_counts",
        "n_features",
    ]

    normal = epithelial.loc[epithelial["tissue"] == "Primary_Normal", base_cols].copy()
    normal["comparison_group"] = "Primary normal epithelial"

    tumor_lauren = epithelial.loc[
        (epithelial["tissue"] == "Primary_Tumor") & epithelial["lauren"].isin(["Intestinal", "Diffuse"]),
        base_cols,
    ].copy()
    tumor_lauren["comparison_group"] = tumor_lauren["lauren"].map(
        {
            "Intestinal": "Tumor intestinal patient",
            "Diffuse": "Tumor diffuse patient",
        }
    )
    lauren = pd.concat([normal, tumor_lauren], ignore_index=True)
    lauren["analysis"] = "Normal epithelial vs tumor patients by Lauren type"

    all_cells = lauren.drop_duplicates("cell_id").reset_index(drop=True)
    return {"lauren": lauren, "all": all_cells}


def extract_gene_counts(sample_manifest: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    sample_to_gsm = dict(zip(sample_manifest["sample_id"], sample_manifest["gsm"]))
    wanted_by_sample = {
        sample_id: set(sub["barcode"].astype(str))
        for sample_id, sub in cells.groupby("sample_id")
    }
    rows = []

    with tarfile.open(SOURCE_TAR, "r") as tar:
        for sample_id, wanted_barcodes in wanted_by_sample.items():
            gsm = sample_to_gsm.get(sample_id)
            if gsm is None:
                raise ValueError(f"No GSM found for {sample_id}")
            member_name = f"{gsm}_{sample_id}.csv.gz"
            member = tar.getmember(member_name)
            raw = tar.extractfile(member)
            if raw is None:
                raise ValueError(f"Could not read {member_name}")

            with gzip.GzipFile(fileobj=raw) as gz:
                reader = csv.reader(TextIOWrapper(gz))
                header = next(reader)
                barcodes = header[1:]
                barcode_positions = {
                    idx: barcode
                    for idx, barcode in enumerate(barcodes)
                    if barcode in wanted_barcodes
                }
                gene_row = None
                for line in reader:
                    if line and line[0].strip('"') == GENE:
                        gene_row = line[1:]
                        break

            if gene_row is None:
                for barcode in wanted_barcodes:
                    rows.append({"sample_id": sample_id, "barcode": barcode, "raw_count": 0})
            else:
                for idx, barcode in barcode_positions.items():
                    rows.append(
                        {
                            "sample_id": sample_id,
                            "barcode": barcode,
                            "raw_count": int(float(gene_row[idx])),
                        }
                    )

    return pd.DataFrame(rows)


def attach_expression(cells: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    out = cells.merge(counts, on=["sample_id", "barcode"], how="left")
    out["raw_count"] = out["raw_count"].fillna(0).astype(int)
    out["detected"] = out["raw_count"] > 0
    out["log1p_cpm"] = np.log1p((out["raw_count"] / out["total_counts"].clip(lower=1)) * 10000)
    return out


def summarize_group(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("comparison_group", sort=False):
        positive = int(sub["detected"].sum())
        total = int(len(sub))
        pos_sub = sub.loc[sub["detected"]]
        rows.append(
            {
                "group": group,
                "n_cells": total,
                "n_patients": sub["patient"].nunique(),
                "detected_cells": positive,
                "detection_pct": 100 * positive / total if total else np.nan,
                "raw_transcripts_total": int(sub["raw_count"].sum()),
                "mean_raw_per_cell": sub["raw_count"].mean(),
                "mean_log1p_cpm": sub["log1p_cpm"].mean(),
                "mean_raw_among_detected": pos_sub["raw_count"].mean() if len(pos_sub) else 0,
                "mean_log1p_cpm_among_detected": pos_sub["log1p_cpm"].mean() if len(pos_sub) else 0,
            }
        )
    return pd.DataFrame(rows)


def summarize_patient(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (group, patient), sub in df.groupby(["comparison_group", "patient"], sort=False):
        rows.append(
            {
                "group": group,
                "patient": patient,
                "n_cells": len(sub),
                "detected_cells": int(sub["detected"].sum()),
                "detection_pct": 100 * sub["detected"].mean(),
                "raw_transcripts_total": int(sub["raw_count"].sum()),
                "mean_raw_per_cell": sub["raw_count"].mean(),
                "mean_log1p_cpm": sub["log1p_cpm"].mean(),
            }
        )
    return pd.DataFrame(rows)


def pairwise_tests(df: pd.DataFrame, group_order: list[str]) -> pd.DataFrame:
    rows = []
    for i, g1 in enumerate(group_order):
        for g2 in group_order[i + 1 :]:
            a = df.loc[df["comparison_group"] == g1]
            b = df.loc[df["comparison_group"] == g2]
            table = [
                [int(a["detected"].sum()), int((~a["detected"]).sum())],
                [int(b["detected"].sum()), int((~b["detected"]).sum())],
            ]
            _, fisher_p = fisher_exact(table)
            wilcox_p = mannwhitneyu(a["log1p_cpm"], b["log1p_cpm"], alternative="two-sided").pvalue
            rows.append(
                {
                    "comparison": f"{g1} vs {g2}",
                    "group_1": g1,
                    "group_2": g2,
                    "fisher_p_detection": fisher_p,
                    "wilcoxon_p_log1p_cpm": wilcox_p,
                    "group_1_detection_pct": 100 * a["detected"].mean(),
                    "group_2_detection_pct": 100 * b["detected"].mean(),
                    "group_1_mean_log1p_cpm": a["log1p_cpm"].mean(),
                    "group_2_mean_log1p_cpm": b["log1p_cpm"].mean(),
                }
            )
    out = pd.DataFrame(rows)
    out["fisher_fdr_detection"] = bh_adjust(out["fisher_p_detection"].tolist())
    out["wilcoxon_fdr_log1p_cpm"] = bh_adjust(out["wilcoxon_p_log1p_cpm"].tolist())
    return out


def patient_level_tests(patient_df: pd.DataFrame, group_order: list[str]) -> pd.DataFrame:
    rows = []
    for i, g1 in enumerate(group_order):
        for g2 in group_order[i + 1 :]:
            a = patient_df.loc[patient_df["group"] == g1]
            b = patient_df.loc[patient_df["group"] == g2]
            rows.append(
                {
                    "comparison": f"{g1} vs {g2}",
                    "group_1": g1,
                    "group_2": g2,
                    "n_patients_group_1": len(a),
                    "n_patients_group_2": len(b),
                    "mannwhitney_p_detection_pct": mannwhitneyu(
                        a["detection_pct"], b["detection_pct"], alternative="two-sided"
                    ).pvalue,
                    "mannwhitney_p_mean_log1p_cpm": mannwhitneyu(
                        a["mean_log1p_cpm"], b["mean_log1p_cpm"], alternative="two-sided"
                    ).pvalue,
                    "median_detection_pct_group_1": a["detection_pct"].median(),
                    "median_detection_pct_group_2": b["detection_pct"].median(),
                    "median_mean_log1p_cpm_group_1": a["mean_log1p_cpm"].median(),
                    "median_mean_log1p_cpm_group_2": b["mean_log1p_cpm"].median(),
                }
            )
    out = pd.DataFrame(rows)
    out["fdr_detection_pct"] = bh_adjust(out["mannwhitney_p_detection_pct"].tolist())
    out["fdr_mean_log1p_cpm"] = bh_adjust(out["mannwhitney_p_mean_log1p_cpm"].tolist())
    return out


def plot_analysis(df: pd.DataFrame, group_order: list[str], title: str, pdf: PdfPages) -> None:
    labels = {
        "Primary normal epithelial": "Primary normal\nepithelial",
        "Tumor intestinal-like": "Tumor\nintestinal-like",
        "Tumor diffuse/EMT-like": "Tumor\ndiffuse/EMT-like",
        "Tumor intestinal patient": "Tumor\nintestinal patient",
        "Tumor diffuse patient": "Tumor\ndiffuse patient",
    }
    summary = summarize_group(df).set_index("group").loc[group_order].reset_index()
    palette = [COLORS[g] for g in group_order]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    x_labels = [labels[g] for g in group_order]

    axes[0].bar(x_labels, summary["detection_pct"], color=palette, edgecolor="black", linewidth=0.4)
    for idx, value in enumerate(summary["detection_pct"]):
        axes[0].text(idx, value + 1, f"{value:.1f}%", ha="center", fontsize=9, fontweight="bold")
    axes[0].set_ylabel("KCNQ1+ cells (%)")
    axes[0].set_title("Detection rate")
    axes[0].set_ylim(0, max(5, summary["detection_pct"].max() * 1.25))

    axes[1].bar(x_labels, summary["mean_log1p_cpm"], color=palette, edgecolor="black", linewidth=0.4)
    axes[1].set_ylabel("Mean log1p(CPM)")
    axes[1].set_title("Mean normalized expression")

    sns.violinplot(
        data=df,
        x="comparison_group",
        y="log1p_cpm",
        order=group_order,
        palette=palette,
        cut=0,
        inner=None,
        linewidth=0.3,
        ax=axes[2],
    )
    sns.boxplot(
        data=df,
        x="comparison_group",
        y="log1p_cpm",
        order=group_order,
        width=0.18,
        showcaps=False,
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.5},
        whiskerprops={"linewidth": 0.5},
        medianprops={"color": "black", "linewidth": 0.7},
        showfliers=False,
        ax=axes[2],
    )
    axes[2].set_xticklabels(x_labels)
    axes[2].set_xlabel("")
    axes[2].set_ylabel("log1p(CPM)")
    axes[2].set_title("Cell-level distribution")

    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    patient_df = summarize_patient(df)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    sns.stripplot(
        data=patient_df,
        x="group",
        y="detection_pct",
        order=group_order,
        palette=palette,
        size=6,
        jitter=0.15,
        edgecolor="black",
        linewidth=0.3,
        ax=axes[0],
    )
    sns.boxplot(
        data=patient_df,
        x="group",
        y="detection_pct",
        order=group_order,
        color="white",
        width=0.35,
        showfliers=False,
        ax=axes[0],
    )
    axes[0].set_xticklabels(x_labels)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Patient-level KCNQ1+ cells (%)")
    axes[0].set_title("Detection by patient")

    sns.stripplot(
        data=patient_df,
        x="group",
        y="mean_log1p_cpm",
        order=group_order,
        palette=palette,
        size=6,
        jitter=0.15,
        edgecolor="black",
        linewidth=0.3,
        ax=axes[1],
    )
    sns.boxplot(
        data=patient_df,
        x="group",
        y="mean_log1p_cpm",
        order=group_order,
        color="white",
        width=0.35,
        showfliers=False,
        ax=axes[1],
    )
    axes[1].set_xticklabels(x_labels)
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Patient-level mean log1p(CPM)")
    axes[1].set_title("Mean expression by patient")
    fig.suptitle(title + " - patient summaries", fontsize=14, fontweight="bold")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sample_manifest, epithelial, final = read_inputs()
    cell_sets = build_cell_sets(epithelial, final)
    counts = extract_gene_counts(sample_manifest, cell_sets["all"])

    lauren_df = attach_expression(cell_sets["lauren"], counts)

    lauren_order = [
        "Primary normal epithelial",
        "Tumor intestinal patient",
        "Tumor diffuse patient",
    ]

    lauren_summary = summarize_group(lauren_df)
    lauren_patient = summarize_patient(lauren_df)
    lauren_tests = pairwise_tests(lauren_df, lauren_order)
    lauren_patient_tests = patient_level_tests(lauren_patient, lauren_order)

    lauren_df.to_csv(OUT_DIR / "kcnq1_primary_normal_vs_lauren_tumor_cell_level.tsv", sep="\t", index=False)
    lauren_summary.to_csv(OUT_DIR / "kcnq1_primary_normal_vs_lauren_tumor_summary.tsv", sep="\t", index=False)
    lauren_patient.to_csv(OUT_DIR / "kcnq1_primary_normal_vs_lauren_tumor_by_patient.tsv", sep="\t", index=False)
    lauren_tests.to_csv(OUT_DIR / "kcnq1_primary_normal_vs_lauren_tumor_tests.tsv", sep="\t", index=False)
    lauren_patient_tests.to_csv(OUT_DIR / "kcnq1_primary_normal_vs_lauren_tumor_patient_tests.tsv", sep="\t", index=False)

    with pd.ExcelWriter(OUT_DIR / "KCNQ1_primary_normal_vs_tumor_comparison.xlsx", engine="xlsxwriter") as writer:
        lauren_summary.to_excel(writer, sheet_name="Lauren_summary", index=False)
        lauren_tests.to_excel(writer, sheet_name="Lauren_cell_tests", index=False)
        lauren_patient.to_excel(writer, sheet_name="Lauren_patient", index=False)
        lauren_patient_tests.to_excel(writer, sheet_name="Lauren_patient_tests", index=False)

    pdf_path = OUT_DIR / "KCNQ1_primary_normal_vs_tumor_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        plot_analysis(
            lauren_df,
            lauren_order,
            "KCNQ1 expression: primary normal epithelial vs tumor samples by Lauren type",
            pdf,
        )

    print("Wrote", pdf_path)
    print("Wrote", OUT_DIR / "KCNQ1_primary_normal_vs_tumor_comparison.xlsx")
    print("\nLauren sample-level summary")
    print(lauren_summary.to_string(index=False))


if __name__ == "__main__":
    main()
