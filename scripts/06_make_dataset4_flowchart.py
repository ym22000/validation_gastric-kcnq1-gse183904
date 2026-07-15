"""Create an A4 portrait flowchart summarizing the GSE183904 dataset4 pipeline."""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "final" / "figures" / "GSE183904_dataset4_pipeline_flowchart_A4.pdf"

PAGE_W, PAGE_H = A4
MARGIN = 28
GREY_DARK = colors.HexColor("#4a4a4a")
GREY_MID = colors.HexColor("#bdbdbd")
GREY_LIGHT = colors.HexColor("#c5edac")
GREY_HEADER = colors.HexColor("#aee38c")
GREY_TABLE = colors.HexColor("#f4f4f4")
BLACK = colors.HexColor("#111111")


def draw_wrapped(c, text, x, y, max_width, font="Helvetica", size=7.0, leading=8.2, bold=False):
    c.setFont("Helvetica-Bold" if bold else font, size)
    approx = max(12, int(max_width / (size * 0.47)))
    lines = []
    for part in str(text).split("\n"):
        lines.extend(wrap(part, width=approx) or [""])
    for line in lines:
        c.drawString(x, y, line)
        y -= leading
    return y


def box(c, x, y, w, h, title, body, fill=GREY_LIGHT, title_fill=GREY_HEADER):
    c.setStrokeColor(GREY_DARK)
    c.setLineWidth(0.8)
    c.setFillColor(fill)
    c.roundRect(x, y, w, h, 7, fill=1, stroke=1)
    c.setFillColor(title_fill)
    c.roundRect(x, y + h - 17, w, 17, 7, fill=1, stroke=0)
    c.setFillColor(BLACK)
    c.setFont("Helvetica-Bold", 7.4)
    c.drawString(x + 7, y + h - 12, title)
    draw_wrapped(c, body, x + 7, y + h - 24, w - 14, size=6.6, leading=7.45)


def arrow(c, x1, y1, x2, y2):
    c.setStrokeColor(GREY_DARK)
    c.setLineWidth(1.0)
    c.line(x1, y1, x2, y2)
    if y2 < y1:
        c.line(x2, y2, x2 - 4, y2 + 6)
        c.line(x2, y2, x2 + 4, y2 + 6)
    elif x2 > x1:
        c.line(x2, y2, x2 - 6, y2 - 4)
        c.line(x2, y2, x2 - 6, y2 + 4)


def result_table(c, x, y, w, h):
    rows = [
        ["Gene", "Detection Int.", "Detection Diff.", "Raw counts Int./Diff.", "Fisher FDR", "Direction"],
        ["KCNQ1", "25.7%", "1.7%", "1714 / 46", "7.1e-178", "intestinal"],
        ["KCNE2", "4.1%", "0.8%", "316 / 54", "3.5e-17", "weak intestinal"],
        ["KCNE3", "28.2%", "5.7%", "2282 / 190", "1.4e-122", "intestinal"],
    ]
    widths = [0.13, 0.16, 0.16, 0.22, 0.16, 0.17]
    c.setStrokeColor(GREY_DARK)
    c.setLineWidth(0.8)
    c.setFillColor(GREY_LIGHT)
    c.roundRect(x, y, w, h, 7, fill=1, stroke=1)
    c.setFillColor(GREY_HEADER)
    c.roundRect(x, y + h - 18, w, 18, 7, fill=1, stroke=0)
    c.setFillColor(BLACK)
    c.setFont("Helvetica-Bold", 7.4)
    c.drawString(x + 7, y + h - 12, "Final KCN readout in high-confidence malignant states")

    top = y + h - 25
    row_h = 13
    c.setFont("Helvetica-Bold", 5.9)
    xx = x + 6
    for label, frac in zip(rows[0], widths):
        c.drawString(xx, top, label)
        xx += w * frac
    c.setStrokeColor(GREY_MID)
    c.line(x + 5, top - 3, x + w - 5, top - 3)
    c.setFont("Helvetica", 5.9)
    for idx, row in enumerate(rows[1:], start=1):
        yy = top - idx * row_h
        if idx % 2:
            c.setFillColor(GREY_TABLE)
            c.rect(x + 5, yy - 3, w - 10, row_h, fill=1, stroke=0)
        c.setFillColor(BLACK)
        xx = x + 6
        for label, frac in zip(row, widths):
            c.drawString(xx, yy, label)
            xx += w * frac

    c.setFont("Helvetica", 5.8)
    c.setFillColor(BLACK)
    c.drawString(x + 7, y + 8, "Patient-aware Spearman supports KCNQ1 and KCNE3 intestinal orientation; KCNE2 remains low and less robust.")


def cell_flow_table(c, x, y, w, h):
    rows = [
        ["Selected GEO matrices", "113,470 cells", "primary tumors + matched normals"],
        ["Broad epithelial candidates", "28,845 cells", "relaxed epithelial gate"],
        ["Malignant epithelial candidates", "11,328 cells", "tumor program + inferCNV rescue"],
        ["After immune-contamination filtering", "8,901 cells", "immune marker threshold >=2 removed"],
        ["Final high-confidence states", "6,573 cells", "4,170 intestinal-like; 2,403 diffuse/EMT-like"],
    ]
    c.setStrokeColor(GREY_DARK)
    c.setLineWidth(0.8)
    c.setFillColor(GREY_LIGHT)
    c.roundRect(x, y, w, h, 7, fill=1, stroke=1)
    c.setFillColor(GREY_HEADER)
    c.roundRect(x, y + h - 17, w, 17, 7, fill=1, stroke=0)
    c.setFillColor(BLACK)
    c.setFont("Helvetica-Bold", 7.4)
    c.drawString(x + 7, y + h - 12, "Cell-flow metrics retained in the final analysis")

    row_h = 10.2
    top = y + h - 28
    col1 = x + 8
    col2 = x + 205
    col3 = x + 315
    c.setFont("Helvetica", 6.1)
    for idx, row in enumerate(rows):
        yy = top - idx * row_h
        if idx % 2 == 0:
            c.setFillColor(GREY_TABLE)
            c.rect(x + 5, yy - 3, w - 10, row_h, fill=1, stroke=0)
        c.setFillColor(BLACK)
        c.drawString(col1, yy, row[0])
        c.drawString(col2, yy, row[1])
        c.drawString(col3, yy, row[2])


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(OUT), pagesize=A4)
    c.setTitle("GSE183904 dataset4 pipeline flowchart")

    c.setFillColor(GREY_HEADER)
    c.roundRect(MARGIN, PAGE_H - 72, PAGE_W - 2 * MARGIN, 45, 8, fill=1, stroke=0)
    c.setFillColor(BLACK)
    c.setFont("Helvetica-Bold", 13)
    c.drawString(MARGIN + 10, PAGE_H - 45, "GSE183904 gastric cancer dataset4 - KCN intestinal vs diffuse/EMT pipeline")
    c.setFont("Helvetica", 7.2)
    c.drawString(MARGIN + 10, PAGE_H - 59, "Biological question: are KCNQ1, KCNE2 and KCNE3 preferentially associated with malignant intestinal-like or diffuse/EMT-like gastric cancer cells?")

    x = MARGIN
    w = PAGE_W - 2 * MARGIN
    y = PAGE_H - 126
    h = 42
    box(c, x, y, w, h, "1. Public dataset and biological scope",
        "GSE183904, Kumar et al., Cancer Discovery 2022, DOI 10.1158/2159-8290.CD-21-0683. Data type: author-processed raw scRNA-seq gene-count matrices from Cell Ranger. Selected scope: primary tumors with Lauren label plus matched normals for CNA reference.",
        fill=GREY_LIGHT)

    y2 = y - 58
    col_gap = 10
    col_w = (w - col_gap) / 2
    box(c, x, y2, col_w, 49, "2A. Samples retained",
        "20 primary tumors: 14 intestinal and 6 diffuse. 5 matched primary normal samples used only as inferCNV reference. Excluded: mixed tumors, metastases, peritoneal samples and missing Lauren labels.")
    box(c, x + col_w + col_gap, y2, col_w, 49, "2B. Input matrix",
        "Compressed CSV count matrices inside GSE183904_RAW.tar. Gene x cell raw UMI counts streamed sample-by-sample and exported as sparse Matrix Market after epithelial candidate selection.")
    arrow(c, x + w / 2, y, x + w / 2, y2 + 49)

    flow = [
        ("3. Broad epithelial candidate gate",
         "Rule: >=3 core epithelial genes, or >=2 core epithelial genes with epithelial score > 65% of strongest non-epithelial score. Core genes include EPCAM, KRT7/8/18/19/20, MUC1, CDH1, TACSTD2. Output: 28,845 / 113,470 cells."),
        ("4. Seurat representation",
         "LogNormalize scale factor 10,000; 2,000 variable genes; ScaleData regressing total UMI; PCA 30 components; neighbors/t-SNE on PCs 1-20; k=30; resolution=0.8; seed=12345."),
        ("5. RNA-inferred CNA support",
         "inferCNV using 750 matched-normal epithelial reference cells, balanced as 150 cells x 5 normal patients. GRCh38 gene order; cutoff=0.1; denoise=TRUE; HMM=FALSE; 8 threads; high CNA threshold = normal 95th percentile, 0.0552533."),
        ("6. Malignant epithelial selection",
         "Tumor and non-malignant epithelial programs from Zhou et al. Supp. Table S4. Primary rule: tumor sample and tumor-program difference > normal 95th percentile, 1.10733. 8,543 selected by tumor program + 2,785 rescued by inferCNV. After removing 2,427 immune-contaminated cells: 8,901 malignant epithelial candidates."),
        ("7. Intestinal vs diffuse/EMT scoring",
         "UCell max rank=1,500. Intestinal signature: 22 genes. Diffuse/EMT score = z(36 EMT-up genes) - z(10 epithelial-junction genes). KCNQ1, KCNE2 and KCNE3 are absent from classifiers. High-confidence rule: abs(delta) >= 0.5 SD and dominant score >= 0."),
        ("8. Final high-confidence cells",
         "6,573 final malignant epithelial cells: 4,170 intestinal-like, 2,403 diffuse/EMT-like, 2,328 indeterminate excluded. Lauren sanity check: diffuse patients have higher median diffuse/EMT-like fraction, 0.670 vs 0.289; Wilcoxon p=0.0433."),
        ("9. KCN tests",
         "Detection = raw count >=1. Two-sided Fisher exact test compares expressing vs non-expressing cells; BH FDR across 3 KCN genes. Continuous association: Spearman on normalized expression, plus within-patient rank permutation with 1,000 permutations."),
    ]

    start_y = y2 - 54
    box_h = 48
    box_gap = 8
    current_y = start_y
    for title, body in flow:
        box(c, x, current_y, w, box_h, title, body)
        next_y = current_y - box_gap - box_h
        if title != flow[-1][0]:
            arrow(c, x + w / 2, current_y, x + w / 2, next_y + box_h)
        current_y = next_y

    cell_flow_table(c, x, 145, w, 74)
    arrow(c, x + w / 2, current_y + box_gap + box_h, x + w / 2, 219)
    result_table(c, x, 54, w, 82)

    c.setFont("Helvetica", 5.6)
    c.setFillColor(GREY_DARK)
    c.drawString(MARGIN, 34, "Tools: Python 3.12.13, numpy 2.3.5, pandas 3.0.1; R 4.5.3, Seurat 5.4.0, inferCNV 1.23.0, UCell 2.14.0, Matrix 1.7.4, dplyr 1.2.0, ggplot2 4.0.2, patchwork 1.3.2.")
    c.drawRightString(PAGE_W - MARGIN, 22, "A4 portrait workflow summary")

    c.showPage()
    c.save()
    print(OUT)


if __name__ == "__main__":
    main()
