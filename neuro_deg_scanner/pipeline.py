"""
neuro-deg-scanner · pipeline.py
Core analysis pipeline: fetch → normalise → DEG → plots → report
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')


# ── Colour palette ────────────────────────────────────────────────────────────
PALETTE = {
    'up'      : '#E05C5C',
    'down'    : '#5C9BE0',
    'ns'      : '#CCCCCC',
    'healthy' : '#5C9BE0',
    'disease' : '#E05C5C',
}


def fetch_dataset(geo_id: str, destdir: str = './geo_cache') -> object:
    """Download a GEO dataset by accession number."""
    try:
        import GEOparse
    except ImportError:
        raise ImportError(
            "GEOparse is required. Install it with: pip install GEOparse"
        )
    os.makedirs(destdir, exist_ok=True)
    print(f"[1/5] Fetching {geo_id} from NCBI GEO...")
    gse = GEOparse.get_GEO(geo=geo_id, destdir=destdir, silent=True)
    print(f"      Title   : {gse.metadata['title'][0]}")
    print(f"      Samples : {len(gse.gsms)}")
    return gse


def build_expression_matrix(gse) -> pd.DataFrame:
    """Build a probes × samples expression matrix from a GSE object."""
    print("[2/5] Building expression matrix...")
    expr_dict = {}
    for gsm_id, gsm in gse.gsms.items():
        if gsm.table is not None and not gsm.table.empty:
            expr_dict[gsm_id] = gsm.table.set_index('ID_REF')['VALUE']

    expr_df = pd.DataFrame(expr_dict)
    expr_df = expr_df.apply(pd.to_numeric, errors='coerce').dropna()
    print(f"      Matrix  : {expr_df.shape[0]:,} probes × {expr_df.shape[1]} samples")
    return expr_df


def extract_labels(
    gse,
    expr_df: pd.DataFrame,
    disease_keywords: list = None,
    healthy_keywords: list = None,
) -> pd.Series:
    """
    Extract sample labels from GEO metadata.

    Parameters
    ----------
    disease_keywords : list, optional
        Extra keywords that identify disease samples.
    healthy_keywords : list, optional
        Extra keywords that identify healthy/control samples.
    """
    DEFAULT_HEALTHY  = ['healthy', 'control', 'normal', 'hc', 'ctrl']
    DEFAULT_DISEASE  = ['patient', 'disease', 'case', 'affected']

    healthy_kw = DEFAULT_HEALTHY  + (healthy_keywords or [])
    disease_kw = DEFAULT_DISEASE  + (disease_keywords or [])

    labels = {}
    for gsm_id, gsm in gse.gsms.items():
        if gsm_id not in expr_df.columns:
            continue
        combined = ' '.join([
            gsm.metadata.get('title',             [''])[0],
            gsm.metadata.get('source_name_ch1',   [''])[0],
            *gsm.metadata.get('characteristics_ch1', [])
        ]).lower()

        if any(w in combined for w in healthy_kw):
            labels[gsm_id] = 'Healthy'
        elif any(w in combined for w in disease_kw):
            labels[gsm_id] = 'Disease'
        else:
            labels[gsm_id] = 'Unknown'

    label_series = pd.Series(labels)
    counts = label_series.value_counts()
    print(f"      Labels  : {counts.get('Disease', 0)} disease / "
          f"{counts.get('Healthy', 0)} healthy / "
          f"{counts.get('Unknown', 0)} unknown")

    keep         = [s for s in expr_df.columns if labels.get(s) in ['Disease', 'Healthy']]
    label_series = label_series[keep]
    return label_series


def normalise_and_filter(
    expr_df: pd.DataFrame,
    label_series: pd.Series,
    variance_percentile: float = 75.0,
) -> pd.DataFrame:
    """
    Log2-transform if needed, then keep top variance probes.

    Parameters
    ----------
    variance_percentile : float
        Keep probes above this variance percentile (default 75 = top 25%).
    """
    expr = expr_df[label_series.index]

    # Log2 transform
    median_val = np.nanmedian(expr.values)
    if median_val > 100:
        expr = np.log2(expr + 1)
        print(f"      log2(x+1) applied  (median was {median_val:.1f})")
    else:
        print(f"      Already log-scale  (median = {median_val:.2f})")

    # Variance filter
    gene_var      = expr.var(axis=1)
    threshold     = gene_var.quantile(variance_percentile / 100)
    expr_filtered = expr[gene_var >= threshold]

    print(f"      Probes  : {expr.shape[0]:,} → {expr_filtered.shape[0]:,} "
          f"(top {100 - variance_percentile:.0f}% by variance)")
    return expr_filtered


def run_differential_expression(
    expr_filtered: pd.DataFrame,
    label_series: pd.Series,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
) -> pd.DataFrame:
    """Welch t-test + Benjamini-Hochberg FDR correction across all probes."""
    print("[3/5] Running differential expression...")

    disease_samples = label_series[label_series == 'Disease'].index
    healthy_samples = label_series[label_series == 'Healthy'].index

    d_expr = expr_filtered[disease_samples]
    h_expr = expr_filtered[healthy_samples]

    # Vectorised Welch t-test
    _, pvals = stats.ttest_ind(d_expr.values, h_expr.values, axis=1, equal_var=False)
    log2fc   = d_expr.mean(axis=1).values - h_expr.mean(axis=1).values

    _, padj, _, _ = multipletests(pvals, method='fdr_bh')

    results = pd.DataFrame({
        'probe'      : expr_filtered.index,
        'log2FC'     : log2fc,
        'pval'       : pvals,
        'padj'       : padj,
        '-log10padj' : -np.log10(np.clip(padj, 1e-300, None)),
    })

    results['significance'] = 'Not significant'
    results.loc[
        (results['padj'] < fdr_threshold) & (results['log2FC'] >  lfc_threshold),
        'significance'
    ] = 'Upregulated'
    results.loc[
        (results['padj'] < fdr_threshold) & (results['log2FC'] < -lfc_threshold),
        'significance'
    ] = 'Downregulated'

    sig = results['significance'].value_counts()
    print(f"      Upregulated   : {sig.get('Upregulated',   0):>5,}")
    print(f"      Downregulated : {sig.get('Downregulated', 0):>5,}")
    print(f"      Not sig.      : {sig.get('Not significant', 0):>5,}")
    return results


def plot_pca(
    expr_filtered: pd.DataFrame,
    label_series: pd.Series,
    output_dir: str = '.',
    disease_label: str = 'Disease',
) -> None:
    """PCA scatter + scree plot."""
    X        = expr_filtered.T
    X_scaled = StandardScaler().fit_transform(X)
    pca      = PCA(n_components=10, random_state=42)
    coords   = pca.fit_transform(X_scaled)

    pca_df          = pd.DataFrame(coords[:, :2], columns=['PC1', 'PC2'], index=X.index)
    pca_df['Group'] = label_series[pca_df.index].map(
        {'Disease': disease_label, 'Healthy': 'Healthy'}
    )

    colors  = {disease_label: PALETTE['disease'], 'Healthy': PALETTE['healthy']}
    markers = {disease_label: 'o',                'Healthy': 's'}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    for grp, col in colors.items():
        sub = pca_df[pca_df['Group'] == grp]
        ax.scatter(sub['PC1'], sub['PC2'], c=col, marker=markers[grp],
                   alpha=0.8, s=60, label=grp, edgecolors='white', linewidth=0.5)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)", fontsize=11)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)", fontsize=11)
    ax.set_title(f"PCA — {disease_label} vs Healthy Controls", fontsize=12)
    ax.legend(fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)

    ax2 = axes[1]
    pcs = range(1, 11)
    ax2.bar(pcs, pca.explained_variance_ratio_ * 100, color='#7CB9E8', edgecolor='white')
    ax2.plot(pcs, np.cumsum(pca.explained_variance_ratio_ * 100),
             color=PALETTE['disease'], marker='o', linewidth=1.5, markersize=4, label='Cumulative')
    ax2.set_xlabel('Principal Component', fontsize=11)
    ax2.set_ylabel('Variance Explained (%)', fontsize=11)
    ax2.set_title('Scree Plot', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.set_xticks(list(pcs))
    ax2.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    path = os.path.join(output_dir, 'pca.png')
    plt.savefig(path, bbox_inches='tight', dpi=150)
    plt.show()
    print(f"      Saved  → {path}")


def plot_volcano(
    results: pd.DataFrame,
    output_dir: str = '.',
    disease_label: str = 'Disease',
    top_n: int = 10,
) -> None:
    """Volcano plot with top-N labelled probes."""
    color_map = {
        'Upregulated'    : PALETTE['up'],
        'Downregulated'  : PALETTE['down'],
        'Not significant': PALETTE['ns'],
    }

    fig, ax = plt.subplots(figsize=(10, 8))
    for sig, color in color_map.items():
        sub   = results[results['significance'] == sig]
        alpha = 0.35 if sig == 'Not significant' else 0.8
        size  = 6    if sig == 'Not significant' else 20
        ax.scatter(sub['log2FC'], sub['-log10padj'],
                   c=color, alpha=alpha, s=size,
                   label=f"{sig} (n={len(sub):,})",
                   edgecolors='none', rasterized=True)

    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline( 0.5, color='grey', linestyle=':', linewidth=0.8, alpha=0.7)
    ax.axvline(-0.5, color='grey', linestyle=':', linewidth=0.8, alpha=0.7)

    top = results[results['significance'] != 'Not significant'].nlargest(top_n, '-log10padj')
    texts = []
    for _, row in top.iterrows():
        txt = ax.text(row['log2FC'], row['-log10padj'], row['probe'],
                      fontsize=7.5, fontweight='bold', color='#333333')
        texts.append(txt)

    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
    except ImportError:
        for txt in texts:
            x, y = txt.get_position()
            txt.set_position((x + 0.05, y + 0.3))

    ax.set_xlabel(f'log₂ Fold Change ({disease_label} vs Healthy)', fontsize=12)
    ax.set_ylabel('−log₁₀ adjusted p-value', fontsize=12)
    ax.set_title(f'Volcano Plot — {disease_label} vs Healthy', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left', framealpha=0.8)
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    path = os.path.join(output_dir, 'volcano.png')
    plt.savefig(path, bbox_inches='tight', dpi=150)
    plt.show()
    print(f"      Saved  → {path}")


def plot_heatmap(
    expr_filtered: pd.DataFrame,
    label_series: pd.Series,
    results: pd.DataFrame,
    output_dir: str = '.',
    top_n: int = 40,
) -> None:
    """Clustered heatmap of top DEGs, Z-score normalised."""
    sig_degs = results[results['significance'] != 'Not significant']

    if len(sig_degs) >= 5:
        top_probes = sig_degs.nlargest(top_n, '-log10padj')['probe'].tolist()
    else:
        top_probes = (results
                      .assign(abs_fc=results['log2FC'].abs())
                      .nlargest(top_n, 'abs_fc')['probe'].tolist())
        print("      Note: using top |FC| genes (few passed FDR threshold)")

    healthy_s    = label_series[label_series == 'Healthy'].index
    disease_s    = label_series[label_series == 'Disease'].index
    sample_order = [s for s in list(healthy_s) + list(disease_s)
                    if s in expr_filtered.columns]

    heatmap_data = expr_filtered.loc[top_probes, sample_order]
    heatmap_z    = heatmap_data.apply(zscore, axis=1)

    col_colors = pd.Series(
        [PALETTE['healthy'] if label_series[s] == 'Healthy' else PALETTE['disease']
         for s in sample_order],
        index=sample_order
    )

    g = sns.clustermap(
        heatmap_z,
        col_colors=col_colors,
        col_cluster=False,
        row_cluster=True,
        cmap='RdBu_r',
        center=0, vmin=-2, vmax=2,
        xticklabels=False,
        yticklabels=True,
        figsize=(14, 10),
        cbar_kws={'label': 'Z-score (log₂ expression)'}
    )
    g.ax_heatmap.set_ylabel('Probe', fontsize=10)
    g.fig.suptitle('Top Differentially Expressed Genes', y=1.01, fontsize=13, fontweight='bold')

    healthy_patch = mpatches.Patch(color=PALETTE['healthy'], label='Healthy')
    disease_patch = mpatches.Patch(color=PALETTE['disease'], label='Disease')
    g.ax_col_colors.legend(
        handles=[healthy_patch, disease_patch],
        loc='upper right', bbox_to_anchor=(1.12, 2.5),
        frameon=False, fontsize=10
    )

    path = os.path.join(output_dir, 'heatmap.png')
    plt.savefig(path, bbox_inches='tight', dpi=150)
    plt.show()
    print(f"      Saved  → {path}")


def generate_report(
    results: pd.DataFrame,
    label_series: pd.Series,
    geo_id: str,
    disease_label: str,
    output_dir: str = '.',
) -> None:
    """Print summary + save CSV results."""
    up   = results[results['significance'] == 'Upregulated']
    down = results[results['significance'] == 'Downregulated']

    print("\n" + "=" * 60)
    print("  NEURO-DEG-SCANNER  ·  RESULTS SUMMARY")
    print("=" * 60)
    print(f"  Dataset         : {geo_id}")
    print(f"  Condition       : {disease_label} vs Healthy")
    print(f"  Disease samples : {(label_series == 'Disease').sum()}")
    print(f"  Healthy samples : {(label_series == 'Healthy').sum()}")
    print(f"  Probes tested   : {len(results):,}")
    print(f"  Upregulated     : {len(up):,}   (FDR<0.05, log₂FC>0.5)")
    print(f"  Downregulated   : {len(down):,}   (FDR<0.05, log₂FC<-0.5)")
    print("=" * 60)

    if len(up):
        print("\n  Top 10 upregulated:")
        print(up.nlargest(10, 'log2FC')[['probe', 'log2FC', 'padj']]
                .to_string(index=False, float_format='{:.4f}'.format))

    if len(down):
        print("\n  Top 10 downregulated:")
        print(down.nsmallest(10, 'log2FC')[['probe', 'log2FC', 'padj']]
                .to_string(index=False, float_format='{:.4f}'.format))

    csv_path = os.path.join(output_dir, f'{geo_id}_deg_results.csv')
    results.sort_values('padj').to_csv(csv_path, index=False)
    print(f"\n  Full results → {csv_path}")
    print("=" * 60)


def run(
    geo_id: str,
    disease_label: str = 'Disease',
    disease_keywords: list = None,
    healthy_keywords: list = None,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
    variance_percentile: float = 75.0,
    output_dir: str = './outputs',
    geo_cache: str = './geo_cache',
) -> pd.DataFrame:
    """
    Run the full neuro-deg-scanner pipeline.

    Parameters
    ----------
    geo_id              : GEO accession number e.g. 'GSE26927'
    disease_label       : Display name for the disease group e.g. 'MS'
    disease_keywords    : Keywords identifying disease samples in GEO metadata
    healthy_keywords    : Keywords identifying healthy/control samples
    fdr_threshold       : Adjusted p-value cutoff (default 0.05)
    lfc_threshold       : |log2FC| cutoff (default 0.5)
    variance_percentile : Filter low-variance probes below this percentile (default 75)
    output_dir          : Folder where plots and CSV are saved
    geo_cache           : Folder where GEO files are cached

    Returns
    -------
    pd.DataFrame of differential expression results
    """
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  neuro-deg-scanner  ·  {geo_id}  ·  {disease_label}")
    print(f"{'='*60}\n")

    gse          = fetch_dataset(geo_id, destdir=geo_cache)
    expr_df      = build_expression_matrix(gse)
    label_series = extract_labels(gse, expr_df, disease_keywords, healthy_keywords)

    print("[2/5] Normalising and filtering...")
    expr_filtered = normalise_and_filter(expr_df, label_series, variance_percentile)

    results = run_differential_expression(expr_filtered, label_series, fdr_threshold, lfc_threshold)

    print("[4/5] Generating plots...")
    sns.set_theme(style='whitegrid')
    plt.rcParams['figure.dpi'] = 130
    plot_pca(expr_filtered, label_series, output_dir, disease_label)
    plot_volcano(results, output_dir, disease_label)
    plot_heatmap(expr_filtered, label_series, results, output_dir)

    print("[5/5] Generating report...")
    generate_report(results, label_series, geo_id, disease_label, output_dir)

    return results
