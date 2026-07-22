from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import fisher_exact, mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "outputs" / "final" / "kcnq1_primary_normal_vs_tumor" / "kcnq1_primary_normal_vs_lauren_tumor_cell_level.tsv"
OUT_DIR = ROOT / "outputs" / "final" / "kcnq1_lauren_downsampled"
N_ITER = 500
SEED = 12345

GROUP_ORDER = [
    "Primary normal epithelial",
    "Tumor intestinal patient",
    "Tumor diffuse patient",
]

GROUP_LABELS = {
    "Primary normal epithelial": "Primary normal\nepithelial",
    "Tumor intestinal patient": "Tumor\nintestinal",
    "Tumor diffuse patient": "Tumor\ndiffuse",
}

COLORS = {
    "Primary normal epithelial": "#8ecae6",
    "Tumor intestinal patient": "#2b9348",
    "Tumor diffuse patient": "#ff5714",
}


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
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


def summarize_cells(df):
    rows = []
    for group, sub in df.groupby("comparison_group", sort=False):
        positive = sub.loc[sub["detected"]]
        rows.append(
            {
                "group": group,
                "n_cells": len(sub),
                "n_patients": sub["patient"].nunique(),
                "detected_cells": int(sub["detected"].sum()),
                "detection_pct": 100 * sub["detected"].mean(),
                "raw_transcripts_total": int(sub["raw_count"].sum()),
                "mean_raw_per_cell": sub["raw_count"].mean(),
                "mean_log1p_cp10k_all_cells": sub["log1p_cp10k"].mean(),
                "n_kcnq1_positive": len(positive),
                "mean_raw_among_kcnq1_positive": positive["raw_count"].mean() if len(positive) else 0,
                "median_raw_among_kcnq1_positive": positive["raw_count"].median() if len(positive) else 0,
                "mean_log1p_cp10k_among_kcnq1_positive": positive["log1p_cp10k"].mean() if len(positive) else 0,
                "median_log1p_cp10k_among_kcnq1_positive": positive["log1p_cp10k"].median() if len(positive) else 0,
            }
        )
    return pd.DataFrame(rows)


def pairwise_tests(df):
    rows = []
    for i, g1 in enumerate(GROUP_ORDER):
        for g2 in GROUP_ORDER[i + 1 :]:
            a = df.loc[df["comparison_group"] == g1]
            b = df.loc[df["comparison_group"] == g2]
            a_pos = a.loc[a["detected"]]
            b_pos = b.loc[b["detected"]]

            fisher_table = [
                [int(a["detected"].sum()), int((~a["detected"]).sum())],
                [int(b["detected"].sum()), int((~b["detected"]).sum())],
            ]
            _, fisher_p = fisher_exact(fisher_table)
            wilcox_positive_p = np.nan
            if len(a_pos) > 0 and len(b_pos) > 0:
                wilcox_positive_p = mannwhitneyu(
                    a_pos["log1p_cp10k"], b_pos["log1p_cp10k"], alternative="two-sided"
                ).pvalue
            rows.append(
                {
                    "comparison": f"{g1} vs {g2}",
                    "group_1": g1,
                    "group_2": g2,
                    "fisher_p_detection": fisher_p,
                    "wilcoxon_p_expression_among_kcnq1_positive": wilcox_positive_p,
                    "group_1_detection_pct": 100 * a["detected"].mean(),
                    "group_2_detection_pct": 100 * b["detected"].mean(),
                    "group_1_mean_log1p_cp10k_among_positive": a_pos["log1p_cp10k"].mean(),
                    "group_2_mean_log1p_cp10k_among_positive": b_pos["log1p_cp10k"].mean(),
                }
            )
    out = pd.DataFrame(rows)
    out["fisher_fdr_detection"] = bh_adjust(out["fisher_p_detection"].tolist())
    valid = out["wilcoxon_p_expression_among_kcnq1_positive"].fillna(1).tolist()
    out["wilcoxon_fdr_expression_among_kcnq1_positive"] = bh_adjust(valid)
    return out


def make_downsample(df, n_per_group, rng):
    sampled = []
    for group in GROUP_ORDER:
        sub = df.loc[df["comparison_group"] == group]
        sampled.append(sub.sample(n=n_per_group, replace=False, random_state=int(rng.integers(0, 2**31 - 1))))
    return pd.concat(sampled, ignore_index=True)


def run_iterations(df, n_per_group):
    rng = np.random.default_rng(SEED)
    rows = []
    for iteration in range(1, N_ITER + 1):
        sampled = make_downsample(df, n_per_group, rng)
        for group, sub in sampled.groupby("comparison_group", sort=False):
            positive = sub.loc[sub["detected"]]
            rows.append(
                {
                    "iteration": iteration,
                    "group": group,
                    "n_cells": len(sub),
                    "detected_cells": int(sub["detected"].sum()),
                    "detection_pct": 100 * sub["detected"].mean(),
                    "mean_log1p_cp10k_all_cells": sub["log1p_cp10k"].mean(),
                    "n_kcnq1_positive": len(positive),
                    "mean_log1p_cp10k_among_kcnq1_positive": positive["log1p_cp10k"].mean() if len(positive) else 0,
                    "median_log1p_cp10k_among_kcnq1_positive": positive["log1p_cp10k"].median() if len(positive) else 0,
                }
            )
    return pd.DataFrame(rows)


def summarize_iterations(iter_df):
    metrics = [
        "detection_pct",
        "mean_log1p_cp10k_all_cells",
        "n_kcnq1_positive",
        "mean_log1p_cp10k_among_kcnq1_positive",
        "median_log1p_cp10k_among_kcnq1_positive",
    ]
    rows = []
    for group, sub in iter_df.groupby("group", sort=False):
        for metric in metrics:
            values = sub[metric].dropna()
            rows.append(
                {
                    "group": group,
                    "metric": metric,
                    "median": values.median(),
                    "mean": values.mean(),
                    "ci_2_5": values.quantile(0.025),
                    "ci_97_5": values.quantile(0.975),
                }
            )
    return pd.DataFrame(rows)


def plot_pdf(full_df, representative_df, iter_df, iter_summary, pdf_path):
    palette = [COLORS[g] for g in GROUP_ORDER]
    labels = [GROUP_LABELS[g] for g in GROUP_ORDER]

    with PdfPages(pdf_path) as pdf:
        rep_summary = summarize_cells(representative_df).set_index("group").loc[GROUP_ORDER].reset_index()

        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.7))
        axes[0].bar(labels, rep_summary["detection_pct"], color=palette, edgecolor="black", linewidth=0.4, width=0.48)
        for i, value in enumerate(rep_summary["detection_pct"]):
            axes[0].text(i, value + 0.45, f"{value:.1f}%", ha="center", fontsize=9, fontweight="bold")
        axes[0].set_ylabel("KCNQ1+ cells (%)")
        axes[0].set_title("Detection rate after downsampling")
        axes[0].set_ylim(0, max(rep_summary["detection_pct"]) * 1.22)

        sns.violinplot(
            data=representative_df.loc[representative_df["detected"]],
            x="comparison_group",
            y="log1p_cp10k",
            order=GROUP_ORDER,
            palette=palette,
            cut=0,
            inner=None,
            linewidth=0.4,
            ax=axes[1],
        )
        sns.boxplot(
            data=representative_df.loc[representative_df["detected"]],
            x="comparison_group",
            y="log1p_cp10k",
            order=GROUP_ORDER,
            width=0.18,
            showfliers=False,
            boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.5},
            medianprops={"color": "black", "linewidth": 0.8},
            whiskerprops={"linewidth": 0.5},
            ax=axes[1],
        )
        axes[1].set_xticklabels(labels)
        axes[1].set_xlabel("")
        axes[1].set_ylabel("log1p(CP10K), KCNQ1+ cells only")
        axes[1].set_title("Expression among KCNQ1+ cells")

        sns.boxplot(
            data=iter_df,
            x="group",
            y="detection_pct",
            order=GROUP_ORDER,
            palette=palette,
            showfliers=False,
            ax=axes[2],
        )
        axes[2].set_xticklabels(labels)
        axes[2].set_xlabel("")
        axes[2].set_ylabel("KCNQ1+ cells (%)")
        axes[2].set_title(f"Downsampling stability ({N_ITER} runs)")
        fig.suptitle("KCNQ1: primary normal vs intestinal and diffuse Lauren tumors", fontsize=14, fontweight="bold")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.6, 3.9))
        plot_data = iter_df.copy()
        sns.boxplot(
            data=plot_data,
            x="group",
            y="mean_log1p_cp10k_among_kcnq1_positive",
            order=GROUP_ORDER,
            palette=palette,
            showfliers=False,
            ax=ax,
        )
        ax.set_xticklabels(labels)
        ax.set_xlabel("")
        ax.set_ylabel("Mean log1p(CP10K), KCNQ1+ cells only")
        ax.set_title("Expression among KCNQ1+ cells across downsampling runs")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t")
    df = df.loc[df["comparison_group"].isin(GROUP_ORDER)].copy()
    df["comparison_group"] = pd.Categorical(df["comparison_group"], categories=GROUP_ORDER, ordered=True)

    n_per_group = int(df["comparison_group"].value_counts().min())
    rng = np.random.default_rng(SEED)
    representative = make_downsample(df, n_per_group, rng)
    iter_df = run_iterations(df, n_per_group)
    iter_summary = summarize_iterations(iter_df)

    full_summary = summarize_cells(df)
    representative_summary = summarize_cells(representative)
    representative_tests = pairwise_tests(representative)

    full_summary.to_csv(OUT_DIR / "KCNQ1_lauren_full_group_summary.tsv", sep="\t", index=False)
    representative.to_csv(OUT_DIR / "KCNQ1_lauren_downsampled_representative_cells.tsv", sep="\t", index=False)
    representative_summary.to_csv(OUT_DIR / "KCNQ1_lauren_downsampled_representative_summary.tsv", sep="\t", index=False)
    representative_tests.to_csv(OUT_DIR / "KCNQ1_lauren_downsampled_representative_tests.tsv", sep="\t", index=False)
    iter_df.to_csv(OUT_DIR / "KCNQ1_lauren_downsampling_iterations.tsv", sep="\t", index=False)
    iter_summary.to_csv(OUT_DIR / "KCNQ1_lauren_downsampling_iteration_summary.tsv", sep="\t", index=False)

    with pd.ExcelWriter(OUT_DIR / "KCNQ1_lauren_downsampled_analysis.xlsx", engine="xlsxwriter") as writer:
        pd.DataFrame(
            {
                "field": [
                    "gene",
                    "groups",
                    "downsampling",
                    "iterations",
                    "seed",
                    "expression_test",
                    "detection_test",
                ],
                "value": [
                    "KCNQ1",
                    "Primary normal epithelial; Tumor intestinal patient; Tumor diffuse patient",
                    f"Random downsampling to {n_per_group} cells per group, matching the smallest group",
                    N_ITER,
                    SEED,
                    "Exploratory Wilcoxon/Mann-Whitney on log1p(CP10K), KCNQ1-positive cells only",
                    "Exploratory Fisher exact test on pooled KCNQ1-positive versus KCNQ1-negative cells",
                ],
            }
        ).to_excel(writer, sheet_name="README", index=False)
        full_summary.to_excel(writer, sheet_name="Full_summary", index=False)
        representative_summary.to_excel(writer, sheet_name="Downsample_summary", index=False)
        representative_tests.to_excel(writer, sheet_name="Downsample_tests", index=False)
        iter_summary.to_excel(writer, sheet_name="Iteration_summary", index=False)
        iter_df.to_excel(writer, sheet_name="Iterations", index=False)

    plot_pdf(
        full_df=df,
        representative_df=representative,
        iter_df=iter_df,
        iter_summary=iter_summary,
        pdf_path=OUT_DIR / "KCNQ1_lauren_downsampled_analysis.pdf",
    )

    print(f"Downsampled to {n_per_group} cells per group for {N_ITER} iterations")
    print("\nRepresentative downsample summary")
    print(representative_summary.to_string(index=False))
    print("\nRepresentative tests")
    print(representative_tests.to_string(index=False))
    print("\nIteration summary")
    print(iter_summary.to_string(index=False))


if __name__ == "__main__":
    main()
