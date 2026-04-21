# 🧠 neuro-deg-scanner

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![GEO](https://img.shields.io/badge/Data-NCBI%20GEO-orange)

> A reusable differential gene expression pipeline for neuroimmune and autoimmune disease datasets from NCBI GEO. Point it at any GEO accession number — it handles everything else.

**Author:** Hafsah Shamsi | Microbiology + Data Science | Mithibai College, Mumbai

---

## What it does

Most differential expression pipelines are written once for one dataset and never reused. **neuro-deg-scanner** is built differently — it's a tool you can run on *any* GEO microarray dataset with a single function call.

```python
from neuro_deg_scanner import run

results = run(
    geo_id           = 'GSE26927',       # any NCBI GEO accession
    disease_label    = 'MS',
    disease_keywords = ['ms', 'multiple sclerosis', 'patient'],
    healthy_keywords = ['healthy', 'control'],
)
```

That's it. It fetches the data, normalises it, runs statistics, and produces three publication-quality figures.

---

## Pipeline

```
GEO Accession ID
      │
      ▼
 Fetch via GEOparse (auto-cached)
      │
      ▼
 Build expression matrix (probes × samples)
      │
      ▼
 Auto-label samples from metadata keywords
      │
      ▼
 log2 normalisation + variance filtering (top 25%)
      │
      ├─────────────────────────┐
      ▼                         ▼
 PCA + Scree Plot        Welch t-test across all probes
 (unsupervised check)    + Benjamini-Hochberg FDR
      │                         │
      │                    Volcano Plot
      │                    Heatmap (Z-scored)
      │                         │
      └────────────┬────────────┘
                   ▼
         Summary report + CSV export
```

---

## Outputs

For every run, neuro-deg-scanner saves:

| File | Description |
|---|---|
| `pca.png` | PCA scatter (disease vs healthy) + scree plot |
| `volcano.png` | Volcano plot with top significant probes labelled |
| `heatmap.png` | Clustered heatmap of top 40 DEGs, Z-score normalised |
| `{GEO_ID}_deg_results.csv` | Full results table sorted by adjusted p-value |

---

## Installation

```bash
git clone https://github.com/HafsahShamsi/neuro-deg-scanner.git
cd neuro-deg-scanner
pip install -r requirements.txt
```

---

## Usage

### Minimal

```python
from neuro_deg_scanner import run

results = run(geo_id='GSE26927', disease_label='MS')
```

### Full options

```python
results = run(
    geo_id              = 'GSE26927',
    disease_label       = 'MS',
    disease_keywords    = ['ms', 'multiple sclerosis', 'rrms', 'patient'],
    healthy_keywords    = ['healthy', 'control', 'hc'],
    fdr_threshold       = 0.05,     # adjusted p-value cutoff
    lfc_threshold       = 0.5,      # |log2 fold change| cutoff
    variance_percentile = 75.0,     # keep top 25% most variable probes
    output_dir          = './outputs/GSE26927',
    geo_cache           = './geo_cache',
)
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `geo_id` | str | required | NCBI GEO accession e.g. `'GSE26927'` |
| `disease_label` | str | `'Disease'` | Display name for disease group |
| `disease_keywords` | list | auto | Keywords to identify disease samples in metadata |
| `healthy_keywords` | list | auto | Keywords to identify healthy samples in metadata |
| `fdr_threshold` | float | `0.05` | Adjusted p-value significance cutoff |
| `lfc_threshold` | float | `0.5` | log₂ fold change magnitude cutoff |
| `variance_percentile` | float | `75.0` | Discard probes below this variance percentile |
| `output_dir` | str | `'./outputs'` | Where to save figures and CSV |
| `geo_cache` | str | `'./geo_cache'` | Where to cache downloaded GEO files |

---

## Example datasets

| Disease | GEO ID | Tissue | Notes |
|---|---|---|---|
| Multiple Sclerosis | GSE26927 | Whole blood | RRMS patients vs healthy controls |
| Systemic Lupus | GSE50772 | Whole blood | SLE patients vs healthy controls |
| Parkinson's Disease | GSE7621 | Substantia nigra | PD vs control brain tissue |
| Alzheimer's Disease | GSE5281 | Entorhinal cortex | AD vs control brain tissue |

See the `examples/` folder for ready-to-run notebooks for each.

---

## Project structure

```
neuro-deg-scanner/
│
├── neuro_deg_scanner/
│   ├── __init__.py          # exposes run()
│   └── pipeline.py          # full pipeline logic
│
├── examples/
│   ├── ms_GSE26927.ipynb    # Multiple Sclerosis
│   └── sle_GSE50772.ipynb   # Systemic Lupus
│
├── outputs/                 # generated figures and CSVs (gitignored)
├── geo_cache/               # cached GEO downloads (gitignored)
├── requirements.txt
├── setup.py
└── README.md
```

---

## Why this exists

Differential expression pipelines get rewritten from scratch constantly — one script per paper, per lab, per dataset. This project packages the core logic into a tool anyone can run on any autoimmune/neurological GEO dataset in minutes, without writing boilerplate.

Built from a real MS analysis ([GSE26927 project](https://github.com/HafsahShamsi/ms_gene_expression)) and generalised into a reusable tool. The diseases this targets — MS, lupus, ALS, Parkinson's, Alzheimer's — all involve immune-gene expression signatures in blood or tissue that are meaningfully different from healthy controls, making them well-suited to this approach.

---

## Biological context

In neuroimmune diseases, blood or tissue gene expression reflects systemic immune dysregulation. This tool is designed around that signal:

- **Upregulated in disease** → typically interferon-stimulated genes, inflammatory cytokines, innate immune activation
- **Downregulated in disease** → regulatory T cell markers, immune tolerance genes, neuroprotective signals

The pipeline uses **Welch's t-test** (unequal variance) + **Benjamini-Hochberg FDR correction** — appropriate for small, unbalanced clinical cohorts typical in GEO datasets.

---

## Roadmap

- [ ] GPL annotation mapping — probe IDs → gene symbols
- [ ] Pathway enrichment (KEGG/GO) on significant DEGs
- [ ] Support for RNA-seq count data (DESeq2-style normalisation)
- [ ] Multi-dataset comparison — same disease, different cohorts
- [ ] GWAS integration — flag probes overlapping known risk loci
- [ ] CLI interface — run from terminal without writing Python

---

## Requirements

```
pandas>=2.0
numpy>=1.26
matplotlib>=3.7
seaborn>=0.13
scipy>=1.11
scikit-learn>=1.3
statsmodels>=0.14
GEOparse>=2.0
adjustText>=1.0   # optional, for cleaner volcano labels
```

---

## License

MIT License — free to use, modify, and build on.

---

*Built as part of a computational biology portfolio focused on neuroimmune disease genomics.*  
*Hafsah Shamsi | Mithibai College, Mumbai*
