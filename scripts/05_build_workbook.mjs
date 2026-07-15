import fs from "node:fs/promises";
import path from "node:path";
import zlib from "node:zlib";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(scriptDir, "..");
const finalDir = path.join(root, "outputs", "final");
const previewDir = path.join(root, "tmp", "excel_previews");

function parseValue(value) {
  if (value === "") return null;
  if (value === "TRUE") return true;
  if (value === "FALSE") return false;
  if (/^-?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(value)) return Number(value);
  return value;
}

async function readTsv(relativePath) {
  const fullPath = path.join(root, relativePath);
  let buffer = await fs.readFile(fullPath);
  if (relativePath.endsWith(".gz")) buffer = zlib.gunzipSync(buffer);
  const lines = buffer.toString("utf8").replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  return lines.map((line) => line.split("\t").map(parseValue));
}

const sheetSpecs = [
  ["Cell Flow", "outputs/tables/cell_flow_complete.tsv"],
  ["Samples", "config/sample_manifest.tsv"],
  ["Input by Sample", "outputs/tables/input_and_epithelial_counts_by_sample.tsv"],
  ["State Summary", "outputs/tables/state_summary.tsv"],
  ["Patient States", "outputs/final/tables/final_state_composition_by_patient.tsv"],
  ["Lauren Validation", "outputs/final/tables/state_classifier_lauren_validation.tsv"],
  ["KCN Summary", "outputs/final/tables/KCN_final_summary.tsv"],
  ["KCN Detection", "outputs/final/tables/KCN_detection_by_state.tsv"],
  ["Detection Tests", "outputs/final/tables/KCN_detection_tests.tsv"],
  ["KCN Correlations", "outputs/final/tables/KCN_signature_correlations.tsv"],
  ["Patient Detection", "outputs/final/tables/KCN_detection_by_patient_state.tsv"],
  ["Signatures", "outputs/tables/final_signature_gene_registry.tsv"],
  ["Parameters", "config/analysis_constants.tsv"],
  ["Malignancy", "outputs/tables/malignant_selection_summary.tsv"],
  ["inferCNV", "outputs/tables/infercnv_run_summary.tsv"],
];

const methods = [
  ["Stage", "Method", "Purpose", "Primary parameter or rule"],
  ["Data scope", "GEO author-processed raw counts", "Use primary gastric tumors with usable Lauren diagnosis", "20 tumors: 14 intestinal, 6 diffuse; 5 matched normals for inferCNV"],
  ["Epithelial gate", "Compartment signature gate", "Retain epithelial and plastic epithelial candidates", ">=3 core genes, or >=2 and epithelial score >0.65 x strongest competitor"],
  ["Representation", "Seurat", "Normalize, reduce dimension and visualize", "LogNormalize 10,000; 2,000 HVGs; 30 PCs; neighbors/t-SNE PCs 1:20; k=30; resolution=0.8"],
  ["CNA support", "inferCNV", "Rescue malignant cells with weak epithelial programs", "750 balanced normal references; cutoff=0.1; q95 normal threshold; denoise; HMM off"],
  ["Malignancy", "Zhou malignant/nonmalignant scores plus inferCNV", "Select malignant epithelial cells without using CNA alone", "Tumor-program delta > normal q95, then inferCNV rescue"],
  ["State scores", "UCell", "Score intestinal and diffuse/EMT programs per cell", "maxRank=1,500; confidence margin=0.5 SD"],
  ["Detection and expression", "Fisher exact test plus raw-count summaries", "Compare expressing versus non-expressing cells and report expression intensity", "Two-sided Fisher test; raw transcripts and mean raw counts reported by state"],
  ["Correlation", "Within-patient rank permutation", "Test continuous KCN-signature association beyond patient composition", "1,000 permutations; seed=12345"],
  ["Multiplicity", "Benjamini-Hochberg", "Control false discovery rate", "Applied within each result family"],
];

const references = [
  ["Use", "Reference", "URL"],
  ["Source atlas and counts", "Kumar et al. 2022, Cancer Discovery 12:670-691", "https://doi.org/10.1158/2159-8290.CD-21-0683"],
  ["GEO archive", "GSE183904", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183904"],
  ["Malignant/nonmalignant programs", "Zhou et al. 2023, Cellular and Molecular Life Sciences 80:52", "https://doi.org/10.1007/s00018-023-04702-1"],
  ["Intestinal tumor-state program", "Kim et al. 2022, npj Precision Oncology 6:9", "https://doi.org/10.1038/s41698-022-00251-1"],
  ["Diffuse/EMT tumor-state program", "Tanabe et al. 2014, International Journal of Oncology 44:1953-1960", "https://doi.org/10.3892/ijo.2014.2387"],
  ["RNA-inferred CNA tool", "inferCNV", "https://github.com/broadinstitute/infercnv"],
  ["Single-cell signature scoring", "UCell", "https://doi.org/10.1016/j.csbj.2021.06.043"],
];

const signatureSources = [
  ["Signature", "Role in pipeline", "Source", "DOI or URL", "Note"],
  ["Zhou_malignant", "Malignant epithelial program used to select tumor-like epithelial cells", "Zhou et al. 2023, Cellular and Molecular Life Sciences, Supplementary Table S4", "https://doi.org/10.1007/s00018-023-04702-1", "Genes filtered from the source table before UCell scoring."],
  ["Zhou_nonmalignant", "Non-malignant epithelial counter-program used to avoid keeping normal epithelial cells", "Zhou et al. 2023, Cellular and Molecular Life Sciences, Supplementary Table S4", "https://doi.org/10.1007/s00018-023-04702-1", "Genes filtered from the source table before UCell scoring."],
  ["Intestinal", "Intestinal malignant epithelial state score", "Kim et al. 2022, npj Precision Oncology, GSE150290 intestinal tumor markers and gastric intestinal-lineage genes", "https://doi.org/10.1038/s41698-022-00251-1", "KCNQ1, KCNE2 and KCNE3 are deliberately absent from this classifier."],
  ["EMT_up", "Diffuse/EMT malignant epithelial state score", "Tanabe et al. 2014, International Journal of Oncology diffuse-type gastric cancer EMT signature; supported by Kim et al. 2022 GSE150290 diffuse/EMT controls", "https://doi.org/10.3892/ijo.2014.2387", "KCNQ1, KCNE2 and KCNE3 are deliberately absent from this classifier."],
  ["Epithelial_junction", "Epithelial-junction counter-score subtracted from EMT-up signal", "Canonical epithelial and junction markers used as EMT counter-score; motivated by Tanabe et al. 2014 and Kim et al. 2022", "https://doi.org/10.3892/ijo.2014.2387; https://doi.org/10.1038/s41698-022-00251-1", "Used to reduce the chance that diffuse/EMT calls reflect epithelial preservation rather than EMT activation."],
];

const colors = {
  navy: "#123047",
  teal: "#0A9396",
  red: "#AE2012",
  gold: "#E9C46A",
  pale: "#EAF4F4",
  paleRed: "#FBE9E7",
  light: "#F5F8FA",
  border: "#CBD5E1",
  white: "#FFFFFF",
  text: "#1F2937",
};

function columnName(index) {
  let name = "";
  for (let n = index + 1; n > 0; n = Math.floor((n - 1) / 26)) {
    name = String.fromCharCode(65 + ((n - 1) % 26)) + name;
  }
  return name;
}

function styleDataSheet(sheet, data) {
  const rows = data.length;
  const cols = data[0].length;
  const end = columnName(cols - 1);
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  const used = sheet.getRange(`A1:${end}${rows}`);
  used.format.font = { name: "Aptos", size: 10, color: colors.text };
  used.format.verticalAlignment = "center";
  const header = sheet.getRange(`A1:${end}1`);
  header.format.fill = colors.navy;
  header.format.font = { name: "Aptos", size: 10, bold: true, color: colors.white };
  header.format.wrapText = true;
  header.format.rowHeight = 30;
  header.format.borders = { preset: "doubleBottom", style: "medium", color: colors.teal };
  if (rows > 1) {
    sheet.getRange(`A2:${end}${rows}`).format.borders = {
      insideHorizontal: { style: "thin", color: colors.border },
    };
  }
  for (let c = 0; c < cols; c += 1) {
    const values = data.slice(0, Math.min(rows, 80)).map((row) => String(row[c] ?? ""));
    const longest = Math.max(...values.map((value) => value.length), 8);
    const width = Math.min(Math.max(longest + 2, 11), c === 0 ? 28 : 34);
    sheet.getRange(`${columnName(c)}1:${columnName(c)}${rows}`).format.columnWidth = width;
  }
  const headerValues = data[0].map((value) => String(value).toLowerCase());
  headerValues.forEach((headerValue, c) => {
    const range = sheet.getRange(`${columnName(c)}2:${columnName(c)}${rows}`);
    if (/fraction|rate|percent/.test(headerValue)) range.format.numberFormat = "0.00%";
    else if (/p_value|fdr/.test(headerValue)) range.format.numberFormat = "0.00E+00";
    else if (/rho|odds_ratio|median_|mean_|score|delta/.test(headerValue)) range.format.numberFormat = "0.000";
    else if (/cells|transcripts|features|patients|_n$|^n$/.test(headerValue)) range.format.numberFormat = "#,##0";
  });
  const stateColumn = headerValues.findIndex((value) => value.includes("state"));
  if (stateColumn >= 0 && rows > 1) {
    const stateRange = sheet.getRange(`${columnName(stateColumn)}2:${columnName(stateColumn)}${rows}`);
    stateRange.conditionalFormats.add("containsText", { text: "Intestinal", format: { fill: colors.pale, font: { color: colors.teal, bold: true } } });
    stateRange.conditionalFormats.add("containsText", { text: "Diffuse", format: { fill: colors.paleRed, font: { color: colors.red, bold: true } } });
  }
}

await fs.mkdir(finalDir, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

const workbook = Workbook.create();
const summary = workbook.worksheets.add("Summary");
for (const [name] of sheetSpecs) workbook.worksheets.add(name);
workbook.worksheets.add("Methods");
workbook.worksheets.add("References");
workbook.worksheets.add("Signature Sources");

for (const [name, relativePath] of sheetSpecs) {
  const data = await readTsv(relativePath);
  const sheet = workbook.worksheets.getItem(name);
  sheet.getRangeByIndexes(0, 0, data.length, data[0].length).values = data;
  styleDataSheet(sheet, data);
  sheet.tables.add(`A1:${columnName(data[0].length - 1)}${data.length}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
}

for (const [name, data] of [["Methods", methods], ["References", references], ["Signature Sources", signatureSources]]) {
  const sheet = workbook.worksheets.getItem(name);
  sheet.getRangeByIndexes(0, 0, data.length, data[0].length).values = data;
  styleDataSheet(sheet, data);
  sheet.getRange(`B2:${columnName(data[0].length - 1)}${data.length}`).format.wrapText = true;
}

summary.showGridLines = false;
summary.getRange("A1:H2").merge();
summary.getRange("A1").values = [["GSE183904 | Gastric malignant epithelial states and KCN genes"]];
summary.getRange("A1:H2").format = {
  fill: colors.navy,
  font: { name: "Aptos Display", size: 18, bold: true, color: colors.white },
  verticalAlignment: "center",
};
summary.getRange("A4:B4").merge();
summary.getRange("C4:D4").merge();
summary.getRange("E4:F4").merge();
summary.getRange("G4:H4").merge();
summary.getRange("A4").values = [["Input cells"]];
summary.getRange("C4").values = [["Malignant candidates"]];
summary.getRange("E4").values = [["Intestinal-like"]];
summary.getRange("G4").values = [["Diffuse/EMT-like"]];
summary.getRange("A5:B6").merge();
summary.getRange("C5:D6").merge();
summary.getRange("E5:F6").merge();
summary.getRange("G5:H6").merge();
summary.getRange("A5").formulas = [["='Cell Flow'!B2"]];
summary.getRange("C5").formulas = [["='Cell Flow'!B4"]];
summary.getRange("E5").formulas = [["='State Summary'!B4"]];
summary.getRange("G5").formulas = [["='State Summary'!B2"]];
summary.getRange("A4:H4").format = { fill: colors.teal, font: { bold: true, color: colors.white }, horizontalAlignment: "center" };
summary.getRange("A5:H6").format = { fill: colors.pale, font: { size: 18, bold: true, color: colors.navy }, horizontalAlignment: "center", verticalAlignment: "center", numberFormat: "#,##0", borders: { preset: "outside", style: "thin", color: colors.border } };
summary.getRange("A8:H8").merge();
summary.getRange("A8").values = [["Final KCN result"]];
summary.getRange("A8:H8").format = { fill: colors.navy, font: { bold: true, color: colors.white, size: 12 } };
summary.getRange("A9:H13").values = [
  ["Gene", "Intestinal detection", "Diffuse/EMT detection", "Fisher odds ratio", "Fisher FDR", "Within-patient intestinal rho", "Within-patient diffuse/EMT rho", "Interpretation"],
  ["KCNQ1", 0.256594724220624, 0.0166458593424886, 20.3878609753901, 7.07040189053279e-178, 0.10280944865733, -0.12663798662369, "Strong intestinal-oriented association"],
  ["KCNE2", 0.0412470023980815, 0.00790678318768206, 5.39670833140489, 3.46314242274162e-17, -0.00811183927042548, -0.0515291855762772, "Weaker; intestinal correlation not supported within patients"],
  ["KCNE3", 0.281534772182254, 0.0574282147315855, 6.42963364011674, 1.35139275722745e-122, 0.144056574127286, -0.148437711571662, "Strong intestinal-oriented association"],
  ["Clinical check", null, null, null, 0.043308142810792, null, null, "Diffuse patients have a higher diffuse/EMT-like fraction (median 0.670 vs 0.289)"],
];
summary.getRange("A9:H9").format = { fill: colors.teal, font: { bold: true, color: colors.white }, wrapText: true };
summary.getRange("A10:H13").format.borders = { insideHorizontal: { style: "thin", color: colors.border } };
summary.getRange("B10:C12").format.numberFormat = "0.00%";
summary.getRange("D10:D12").format.numberFormat = "0.00";
summary.getRange("E10:E13").format.numberFormat = "0.00E+00";
summary.getRange("F10:G12").format.numberFormat = "0.000";
summary.getRange("H10:H13").format.wrapText = true;
summary.getRange("A15:H16").merge();
summary.getRange("A15").values = [["Interpretation: KCNQ1 and KCNE3 are preferentially associated with an intestinal malignant epithelial program. KCNE2 follows the detection direction but has weaker within-patient evidence. These are state associations, not causal or functional channel evidence."]];
summary.getRange("A15:H16").format = { fill: colors.gold, font: { color: colors.navy, italic: true }, wrapText: true, verticalAlignment: "center" };
summary.getRange("A18:H18").merge();
summary.getRange("A18").values = [["Expression context from raw counts"]];
summary.getRange("A18:H18").format = { fill: colors.navy, font: { bold: true, color: colors.white, size: 12 } };
summary.getRange("A19:H22").values = [
  ["Gene", "Int raw transcripts", "Diffuse raw transcripts", "Int mean raw/cell", "Diffuse mean raw/cell", "Int mean among expressing", "Diffuse mean among expressing", "Expression reading"],
  ["KCNQ1", 1714, 46, 0.411031175059952, 0.0191427382438618, 1.6018691588785, 1.15, "More often detected and higher total expression in intestinal-like cells"],
  ["KCNE2", 316, 54, 0.075779376498801, 0.0224719101123596, 1.83720930232558, 2.84210526315789, "Low expression overall; diffuse expressing cells can have slightly higher per-positive-cell counts"],
  ["KCNE3", 2282, 190, 0.547242206235012, 0.0790678318768206, 1.94378194207836, 1.3768115942029, "More often detected and higher total expression in intestinal-like cells"],
];
summary.getRange("A19:H19").format = { fill: colors.teal, font: { bold: true, color: colors.white }, wrapText: true };
summary.getRange("A20:H22").format.borders = { insideHorizontal: { style: "thin", color: colors.border } };
summary.getRange("B20:C22").format.numberFormat = "#,##0";
summary.getRange("D20:G22").format.numberFormat = "0.000";
summary.getRange("H20:H22").format.wrapText = true;
summary.getRange("A1:H22").format.font.name = "Aptos";
summary.getRange("A1:A22").format.columnWidth = 18;
summary.getRange("B1:G22").format.columnWidth = 17;
summary.getRange("H1:H22").format.columnWidth = 42;
summary.freezePanes.freezeRows(2);

let outputPath = path.join(finalDir, "GSE183904_KCN_intestinal_diffuse_pipeline_registry.xlsx");
const output = await SpreadsheetFile.exportXlsx(workbook);
try {
  await output.save(outputPath);
} catch (error) {
  if (error && error.code === "EBUSY") {
    outputPath = path.join(finalDir, "GSE183904_KCN_intestinal_diffuse_pipeline_registry_with_expression.xlsx");
    await output.save(outputPath);
  } else {
    throw error;
  }
}

const summaryCheck = await workbook.inspect({
  kind: "table",
  range: "Summary!A1:H22",
  include: "values,formulas",
  tableMaxRows: 20,
  tableMaxCols: 10,
});
console.log(summaryCheck.ndjson);

const errorCheck = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errorCheck.ndjson);

for (const sheetName of ["Summary", ...sheetSpecs.map(([name]) => name), "Methods", "References", "Signature Sources"]) {
  const sheet = workbook.worksheets.getItem(sheetName);
  const used = sheet.getUsedRange(true);
  const preview = await workbook.render({
    sheetName,
    range: used ? undefined : "A1:H20",
    autoCrop: "all",
    scale: sheetName === "Summary" ? 1.4 : 1,
    format: "png",
  });
  await fs.writeFile(path.join(previewDir, `${sheetName.replace(/[^A-Za-z0-9]+/g, "_")}.png`), new Uint8Array(await preview.arrayBuffer()));
}

console.log(`Workbook saved to ${outputPath}`);
