"""Stream GSE183904 CSV matrices and retain broad epithelial candidates."""

from __future__ import annotations

import csv
import gzip
import io
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "sources_dataset4" / "GSE183904_RAW.tar"
MANIFEST_PATH = ROOT / "config" / "sample_manifest.tsv"
TABLE_DIR = ROOT / "outputs" / "tables"
INTERMEDIATE_DIR = ROOT / "outputs" / "intermediate"
MATRIX_PATH = INTERMEDIATE_DIR / "epithelial_counts.mtx"
BODY_PATH = ROOT / "tmp" / "epithelial_counts.entries.tmp"

MARKER_SETS = {
    "epithelial": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT20", "MUC1", "CDH1", "TACSTD2", "KRT17"],
    "gastric_lineage": ["MUC5AC", "MUC6", "TFF1", "TFF2", "PGC", "PGA3", "LIPF"],
    "intestinal_lineage": ["CDH17", "REG4", "MUC13", "TFF3", "CDX1", "CDX2", "KRT20"],
    "immune": ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD79A", "NKG7", "LST1", "TYROBP", "FCER1G", "CD68"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "PDGFRA", "PDGFRB", "COL5A1", "SPARC", "FN1", "S100A4", "VIM"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ENG", "RAMP2", "PLVAP"],
    "pericyte": ["RGS5", "CSPG4", "MCAM", "ACTA2", "TAGLN", "MYL9"],
    "emt": ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "COL1A1", "COL1A2", "S100A4"],
    "malignant_support": ["CEACAM5", "CEACAM6", "MSLN", "MKI67", "TOP2A", "CLDN18", "MUC16", "SOX9", "MYC", "ERBB2"],
}
CORE_EPITHELIAL = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT20", "MUC1", "CDH1", "TACSTD2"]
ALL_MARKERS = sorted({gene for genes in MARKER_SETS.values() for gene in genes})


def parse_values(line: str) -> tuple[str, np.ndarray]:
    gene_field, values_field = line.rstrip("\r\n").split(",", 1)
    gene = gene_field.strip('"').upper()
    return gene, np.fromstring(values_field, dtype=np.int32, sep=",")


def open_member(tar: tarfile.TarFile, member_name: str):
    raw = tar.extractfile(member_name)
    if raw is None:
        raise FileNotFoundError(member_name)
    compressed = gzip.GzipFile(fileobj=raw)
    return io.TextIOWrapper(compressed, newline="")


def score_sample(tar: tarfile.TarFile, member_name: str, sample_row: pd.Series):
    with open_member(tar, member_name) as stream:
        header = next(csv.reader([stream.readline()]))
        barcodes = [cell.strip('"') for cell in header[1:]]
        n_cells = len(barcodes)
        total_counts = np.zeros(n_cells, dtype=np.int64)
        n_features = np.zeros(n_cells, dtype=np.int32)
        mt_counts = np.zeros(n_cells, dtype=np.int64)
        marker_counts = {gene: np.zeros(n_cells, dtype=np.int32) for gene in ALL_MARKERS}
        features = []
        for line in stream:
            gene, values = parse_values(line)
            if values.size != n_cells:
                raise ValueError(f"Unexpected cell count in {member_name}, gene {gene}")
            features.append(gene)
            total_counts += values
            n_features += values > 0
            if gene.startswith("MT-"):
                mt_counts += values
            if gene in marker_counts:
                marker_counts[gene] = values

    normalized = {
        gene: np.log1p(values / np.maximum(total_counts, 1) * 10000)
        for gene, values in marker_counts.items()
    }
    scores = {}
    for group, genes in MARKER_SETS.items():
        scores[group] = np.vstack([normalized[gene] for gene in genes]).mean(axis=0)

    epithelial_detected = np.vstack([marker_counts[gene] > 0 for gene in CORE_EPITHELIAL]).sum(axis=0)
    competitor = np.maximum.reduce([scores["immune"], scores["fibroblast"], scores["endothelial"], scores["pericyte"]])
    epithelial_candidate = (epithelial_detected >= 3) | ((epithelial_detected >= 2) & (scores["epithelial"] > 0.65 * competitor))
    pct_mt = mt_counts / np.maximum(total_counts, 1) * 100

    frame = pd.DataFrame({
        "cell_id": [f"{sample_row.sample_id}_{barcode}" for barcode in barcodes],
        "barcode": barcodes,
        "sample_id": sample_row.sample_id,
        "gsm": sample_row.gsm,
        "patient": sample_row.patient,
        "tissue": sample_row.tissue,
        "lauren": sample_row.lauren,
        "total_counts": total_counts,
        "n_features": n_features,
        "percent_mt": pct_mt,
        "epithelial_detected": epithelial_detected,
        "epithelial_candidate": epithelial_candidate,
    })
    for group, values in scores.items():
        frame[f"{group}_score"] = values
    return frame, features


def write_selected_entries(tar, member_name, selected_indices, column_offset, body_handle, expected_features):
    kept_nnz = 0
    with open_member(tar, member_name) as stream:
        stream.readline()
        for row_index, line in enumerate(stream, start=1):
            gene, values = parse_values(line)
            if gene != expected_features[row_index - 1]:
                raise ValueError(f"Feature order mismatch in {member_name} at row {row_index}")
            selected_values = values[selected_indices]
            nonzero = np.flatnonzero(selected_values)
            for local_column in nonzero:
                body_handle.write(f"{row_index} {column_offset + local_column + 1} {selected_values[local_column]}\n")
            kept_nnz += nonzero.size
    return kept_nnz


def main():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
    BODY_PATH.parent.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
    manifest["analysis_scope"] = manifest["analysis_scope"].astype(bool)
    scope = manifest.loc[manifest.analysis_scope].copy()
    scores = []
    first_features = None

    with tarfile.open(ARCHIVE) as tar:
        members = {Path(name).name.split("_")[1].replace(".csv.gz", ""): name for name in tar.getnames()}
        for row in scope.itertuples(index=False):
            member_name = members[row.sample_id]
            frame, features = score_sample(tar, member_name, row)
            if first_features is None:
                first_features = features
            elif features != first_features:
                raise ValueError(f"Feature order differs in {row.sample_id}")
            scores.append(frame)
            print(f"Scored {row.sample_id}: {len(frame)} cells; {frame.epithelial_candidate.sum()} epithelial candidates", flush=True)

        all_scores = pd.concat(scores, ignore_index=True)
        all_scores.to_csv(TABLE_DIR / "analysis_scope_cell_scores.tsv.gz", sep="\t", index=False, compression="gzip")
        selected = all_scores.loc[all_scores.epithelial_candidate].copy()
        selected["matrix_column"] = np.arange(1, len(selected) + 1)
        selected.to_csv(TABLE_DIR / "epithelial_cell_manifest.tsv", sep="\t", index=False)

        sample_summary = all_scores.groupby(["sample_id", "patient", "tissue", "lauren"], as_index=False).agg(
            input_cells=("cell_id", "size"),
            epithelial_candidates=("epithelial_candidate", "sum"),
            median_features=("n_features", "median"),
            median_counts=("total_counts", "median"),
            median_percent_mt=("percent_mt", "median"),
        )
        sample_summary["epithelial_fraction"] = sample_summary.epithelial_candidates / sample_summary.input_cells
        sample_summary.to_csv(TABLE_DIR / "input_and_epithelial_counts_by_sample.tsv", sep="\t", index=False)

        registry = pd.DataFrame([
            {"signature": group, "gene": gene, "available": gene in set(first_features), "source": "Broad lineage gate adapted from dataset2"}
            for group, genes in MARKER_SETS.items() for gene in genes
        ])
        registry.to_csv(TABLE_DIR / "broad_gate_signature_registry.tsv", sep="\t", index=False)
        pd.DataFrame({"row_index": np.arange(1, len(first_features) + 1), "gene": first_features}).to_csv(
            TABLE_DIR / "features.tsv", sep="\t", index=False
        )

        total_nnz = 0
        column_offset = 0
        with BODY_PATH.open("w", encoding="ascii", newline="\n") as body:
            for row in scope.itertuples(index=False):
                member_name = members[row.sample_id]
                sample_scores = all_scores.loc[all_scores.sample_id == row.sample_id].reset_index(drop=True)
                selected_indices = np.flatnonzero(sample_scores.epithelial_candidate.to_numpy())
                total_nnz += write_selected_entries(
                    tar, member_name, selected_indices, column_offset, body, first_features
                )
                column_offset += selected_indices.size
                print(f"Extracted {row.sample_id}: {selected_indices.size} cells", flush=True)

    with MATRIX_PATH.open("wb") as output:
        output.write(b"%%MatrixMarket matrix coordinate integer general\n")
        output.write(f"{len(first_features)} {len(selected)} {total_nnz}\n".encode("ascii"))
        with BODY_PATH.open("rb") as body:
            shutil.copyfileobj(body, output, length=1024 * 1024)
    BODY_PATH.unlink(missing_ok=True)

    flow = pd.DataFrame([
        {"step": "GEO processed matrices", "cells": int(manifest.merge(sample_summary, on=["sample_id", "patient", "tissue", "lauren"], how="inner").input_cells.sum())},
        {"step": "Broad epithelial candidates", "cells": int(len(selected))},
    ])
    flow.to_csv(TABLE_DIR / "cell_flow_initial.tsv", sep="\t", index=False)
    print(f"Final epithelial matrix: {len(first_features)} genes x {len(selected)} cells; {total_nnz} nonzero entries")


if __name__ == "__main__":
    main()
