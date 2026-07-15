# Critical assessment

## Strengths

- The input is a large public gastric cancer atlas with primary intestinal and diffuse tumors and matched normal samples.
- Malignancy is not assigned from one marker or from inferCNV alone. Tumor transcriptional programs provide the primary evidence and RNA-inferred CNA is used only as a rescue signal.
- The intestinal and diffuse/EMT classifiers do not contain KCNQ1, KCNE2 or KCNE3, avoiding circularity.
- A confidence margin prevents every malignant cell from being forced into one of two states.
- Detection is tested with Fisher exact tests, and continuous associations are checked with within-patient permutations.
- Lauren diagnosis is used only afterward as a clinical sanity check. Diffuse tumors contain a higher median diffuse/EMT fraction, which supports the biological direction.

## Limitations and possible biases

- GEO provides author-processed count matrices rather than FASTQ files. The original alignment and initial QC cannot be rerun here.
- The exact author malignant-cell labels and CONICSmat results are not available at cell level. This is an adapted reconstruction using inferCNV, not a strict reproduction of the publication.
- inferCNV estimates CNA-like RNA patterns. It is not direct DNA copy-number measurement and can be influenced by expression programs or technical variation.
- Broad epithelial preselection can still lose very advanced EMT cells with almost no epithelial signal, or include some activated stromal cells. The inferCNV rescue and immune filter reduce but do not eliminate this risk.
- The final state labels are supervised by predefined signatures. Their separation demonstrates consistency with the chosen programs, not the discovery of completely independent cell classes.
- Single-cell observations from one patient are not independent biological replicates. Fisher exact tests are useful for detection-rate comparison, but they remain cell-level tests; the within-patient rank permutations are therefore important supporting analyses.
- KCN transcripts are sparse and affected by dropout. Absence of a count does not prove absence of expression or channel activity.
- t-SNE preserves local neighborhoods imperfectly and is used only for visualization.

## Strength of the gene conclusions

**KCNQ1: strong exploratory evidence.** The detection difference is large, Fisher exact test strongly supports the difference, and both patient-aware signature directions agree. The data support preferential association with the intestinal-like state.

**KCNE3: strong exploratory evidence.** Detection and both patient-aware correlations agree with an intestinal orientation. The evidence is similar in direction to KCNQ1, although expression of an auxiliary subunit does not establish the composition of a functional channel complex.

**KCNE2: limited evidence.** Detection is higher in intestinal-like cells and lower with the diffuse/EMT score, but the within-patient intestinal correlation is null. It should not be presented as robustly intestinal-specific from this dataset alone.

## Overall conclusion

The pipeline provides solid single-cell evidence for an association of KCNQ1 and KCNE3 with an intestinal malignant epithelial program in GSE183904. It is sufficient for a well-supported exploratory result, not for a causal or biomarker claim. Replication in an independent cohort, spatial or protein validation, and functional experiments remain necessary before publication-level biological causality can be argued.
