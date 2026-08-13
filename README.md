# Single-Cell Analysis of the Breast Cancer Tumor Microenvironment

This repository extends an earlier bulk RNA-seq project — [Breast Cancer Subtype Transcriptomics](https://github.com/sana-unes/breastcancer-subtype-expression) — to single-cell resolution, using the Wu et al. 2021 breast cancer scRNA-seq atlas.

## Research Question

The bulk RNA-seq subtype analysis found a mitotic/cell-cycle gene signature enriched in the more aggressive breast cancer subtype relative to the less aggressive one. Bulk RNA-seq averages expression across every cell in a tumor sample, so that finding alone can't distinguish two different explanations: malignant cells themselves being more proliferative, or the aggressive subtype simply containing a different mix of cell types that makes the tumor look more proliferative on average without any single cell type actually behaving differently. This project asks which explanation holds up once cell type and malignancy status are accounted for directly, at single-cell resolution — and whether a benchmarked set of annotation methods can even make that distinction reliably in the first place.

## Background

Bulk transcriptomics is a measurement of an entire tissue sample at once — every cell type present gets folded into a single averaged signal. A shift in that signal between two conditions can come from individual cells changing their expression, from the underlying mixture of cell types shifting, or both at once, and bulk data has no way to tell these apart. Single-cell RNA-seq resolves expression per cell, which makes it possible to ask the composition question and the cell-intrinsic question separately rather than as one confounded signal. That's the motivation for this extension: not a new question, but the tool needed to actually answer the one the bulk work raised.

## Dataset

- **Source:** Wu et al. 2021, *Nature Genetics* — "A single-cell and spatially resolved atlas of human breast cancers" (GSE176078), accessed via CZ CELLxGENE
- **Scope:** the full scRNA-seq atlas — all patients, all annotated cell types, not a pre-filtered subset
- **Size:** 100,064 cells as downloaded, 98,964 after quality control and doublet removal
- **Patients:** 26 donors across three molecular subtypes — TNBC, ER+, HER2+
- **Cell types:** nine broad categories from the original authors' annotation (T-cells, Cancer Epithelial, Myeloid, Endothelial, CAFs, PVL, Normal Epithelial, Plasmablasts, B-cells)

## Methods

The analysis is organized into three notebooks, each building on a checkpoint saved by the one before it.

**1. Preprocessing and clustering** (`01_preprocessing_and_clustering.ipynb`)
Raw counts and gene symbols are recovered from the CELLxGENE object's `.raw` attribute, since `.X` ships pre-normalized. Quality control thresholds are chosen from the actual observed distributions rather than copied from a general-purpose default. Doublet detection is run per-donor with Scrublet, since doublets form within a single sequencing run and simulating them across pooled patients would produce biologically meaningless pairs. Highly variable gene selection is batch-aware (`batch_key="donor_id"`), and PCA is computed on a restricted HVG subset to avoid densifying the full sparse matrix. Harmony integrates across the 26 donors, and Leiden clustering runs on the batch-corrected representation.

**2. Annotation benchmarking** (`02_annotation_benchmarking.ipynb`)
Three cell-type annotation approaches are compared against the original authors' expert labels: marker-gene scoring, a logistic regression classifier trained directly on the data with a patient-level (not cell-level) train/test split, and CellTypist, a pretrained immune-cell classifier. CellTypist is evaluated only on the immune subset of the ground truth, since its models have no concept of epithelial or stromal identity — scoring it against the full label set would penalize it for a distinction it was never built to make.

**3. Malignant signature analysis** (`03_malignant_signature_analysis.ipynb`)
The benchmarking step showed that none of the three annotation methods reliably separates malignant from normal epithelial cells — that distinction depends on copy-number variation, not marker expression. This notebook instead uses the dataset's own CNA-based `normal_cell_call` field to define the malignant population directly, scores the original bulk mitotic signature within that population only, and tests TNBC vs. ER+ malignant cells at the patient level rather than the cell level, since treating thousands of correlated cells from the same tumor as independent samples would inflate significance artificially.

## Key Results

**Clustering:** 30 Leiden clusters recovered from the harmony-corrected representation. Donor mixing within clusters confirms batch correction worked rather than leaving patient identity as the dominant signal. A handful of very small, single-donor clusters were investigated individually rather than dismissed — each resolved to a single, biologically plausible cell type (PVL, Myeloid, Cancer/Normal Epithelial), consistent with genuine patient-specific biology (particularly for malignant epithelial cells, which carry patient-specific genomic alterations) rather than technical artifacts.

**Annotation benchmarking:**

| Method | Scope | Accuracy | ARI |
|---|---|---|---|
| Marker scoring | all 9 categories | 84.5% | 0.738 |
| Trained classifier (held-out patients) | all 9 categories | 95.4% | 0.950 |
| Marker scoring | immune only | 87.9% | 0.702 |
| Trained classifier (held-out patients) | immune only | 99.2% | 0.986 |
| CellTypist | immune only | 98.8% | 0.978 |

Marker scoring's weakest point is Cancer vs. Normal Epithelial (40.9% accuracy) — a real, informative limitation rather than a coding error, since there's no canonical marker gene for malignancy.

**Malignant signature analysis:** restricted to CNA-confirmed malignant cells and tested at the patient level (n=8 TNBC, n=9 ER+ patients), the original bulk mitotic signature holds up — Mann-Whitney U, p=0.0002, with 9/9 ER+ patients scoring negative and 7/8 TNBC patients scoring positive. Tumor composition differs independently: TNBC tumors show a smaller malignant-cell fraction and roughly double the myeloid infiltration of ER+ tumors.

## Interpretation

The most defensible reading of these results is that the original bulk finding reflects two mechanisms compounding each other rather than one explaining away the other: TNBC malignant cells appear genuinely more proliferative at the single-cell level, and TNBC tumors independently have a different, more myeloid-shifted microenvironment. Bulk RNA-seq could not have distinguished which of these was driving the signal; single-cell resolution, tested correctly at the patient level, distinguishes them directly.

## Limitations

- This dataset's `subtype` field (TNBC / ER+ / HER2+) reflects clinical receptor status, not the PAM50 molecular subtypes (Basal / Luminal A) used in the original bulk comparison. TNBC and Basal overlap substantially in the literature but aren't interchangeable, and ER+ is a broader category that would include both Luminal A and Luminal B. TNBC vs. ER+ is used here as the closest available proxy, not an exact match.
- HER2+ patients (n=3) were reported for context only and excluded from the primary statistical comparison given the small sample size.
- CNA-based malignancy calls are only available for the epithelial compartment, consistent with how the original authors applied copy-number inference.
- CellTypist's evaluation is intentionally scoped to immune cells; it was not tested against epithelial or stromal ground truth, since its pretrained models don't cover those categories.

## Repository Contents

- `01_preprocessing_and_clustering.ipynb` — QC, doublet removal, normalization, batch correction, clustering
- `02_annotation_benchmarking.ipynb` — marker scoring, trained classifier, CellTypist, and their comparison
- `03_malignant_signature_analysis.ipynb` — the core analysis: malignant-cell signature scoring and patient-level statistics
- `environment.yml` — package versions for local reproduction
- `figures/` — exported plots (QC distributions, UMAP overview, confusion matrices, signature score comparisons)

## Tools and Libraries

- Python, scanpy, anndata
- harmonypy (batch correction), leidenalg/igraph (clustering)
- celltypist (automated annotation)
- scikit-learn (classification, evaluation metrics)
- scipy (statistical testing)
- pandas, numpy, matplotlib, seaborn

## Author

**Sana Younes**
