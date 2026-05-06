"""
velocyto_integration.py
=======================
Runs velocyto on per-sample BAM files, then transfers spliced/unspliced layers
into an existing annotated h5ad object (which contains scTransform data).

Samples
-------
  anagen   → sample_alignments_IGF133568.bam
  catagen  → sample_alignments_IGF133569.bam
  telogen  → sample_alignments_IGF133570.bam

Output
------
  annotated_with_velocity.h5ad  — original object + spliced/unspliced layers

Notes on transformed vs raw data
---------------------------------
scVelo MUST run on raw (un-normalised, un-log-transformed) counts.
The RNA velocity model estimates transcriptional kinetics by fitting
spliced/unspliced count ratios — log-normalised or variance-stabilised
values would break the steady-state and dynamical model assumptions.

Your h5ad contains counts_SCT (scTransform corrected counts).
We therefore:
  1. Keep counts_SCT intact for all clustering/DE purposes.
  2. Add raw spliced / unspliced layers for velocity only.
  3. Run scVelo preprocessing ONLY on those raw layers.

Usage
-----
  python velocyto_integration.py \\
      --h5ad       path/to/your_object.h5ad \\
      --bam_dir    path/to/bam_files/ \\
      --gtf        path/to/genome.gtf \\
      --output_dir ./velocity_output \\
      [--repeat_gtf path/to/repeat_masker.gtf] \\
      [--n_jobs 8]
"""

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
from scipy.sparse import csr_matrix
import anndata as ad
import loompy
import scvelo as scv
import matplotlib.pyplot as plt


#################################################################
# Logging
#################################################################

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


#################################################################
# Sample Map
#################################################################

SAMPLE_MAP = {
    "Anagen":  "IGF133568",
    "Catagen": "IGF133569",
    "Telogen": "IGF133570",
}


#################################################################
# 1. Export per-sample barcode whitelists from the h5ad object
#################################################################

def inspect_and_export_barcodes(adata: ad.AnnData, output_dir: Path) -> dict:
    """
    Examine barcodes in the h5ad object and export per-sample whitelists.

    The function tries to infer the barcode prefix format used in the object
    so that barcodes can later be matched against velocyto loom output.

    Returns a dict mapping sample_id -> (barcode_csv_path, prefix_used)
    """
    log.info("=== Barcode Inspection ===")
    log.info(f"Total cells in object : {adata.n_obs}")
    log.info(f"Sample column values  : {adata.obs['orig.ident'].value_counts().to_dict()}")

    example_barcodes = adata.obs_names[:10].tolist()
    log.info(f"Example barcodes from object:\n  " + "\n  ".join(example_barcodes))

    barcode_paths = {}

    for ident, igf_id in SAMPLE_MAP.items():
        mask = adata.obs["orig.ident"] == ident
        sample_barcodes = adata.obs_names[mask]

        log.info(f"\n--- {ident} ({igf_id}) ---")
        log.info(f"  Cell count        : {mask.sum()}")
        log.info(f"  Example barcodes  : {sample_barcodes[:5].tolist()}")

        # Detect prefix format
        # Common patterns:
        #   (a) 'ACGTACGT-1'               -> no prefix, bare 10x barcode
        #   (b) 'IGF133568_ACGTACGT-1'     -> IGF prefix
        #   (c) 'Anagen_ACGTACGT-1'        -> capitalised ident prefix
        #   (d) '_Anagen_ACGTACGT-1'       -> leading underscore + ident prefix
        #   (e) 'anagen_ACGTACGT-1'        -> lowercase ident prefix
        first_bc = sample_barcodes[0]

        # Build candidate prefixes in order of specificity (most specific first)
        candidates = [
            f"_{ident}_",           # '_Anagen_'  (leading underscore)
            f"_{ident.lower()}_",   # '_anagen_'
            f"{igf_id}_",           # 'IGF133568_'
            f"{ident}_",            # 'Anagen_'
            f"{ident.lower()}_",    # 'anagen_'
        ]

        prefix = None
        for candidate in candidates:
            if first_bc.startswith(candidate):
                prefix = candidate
                break

        log.info(f"  Detected prefix   : {repr(prefix)}")

        # Strip prefix to get the raw 10x barcode for velocyto matching
        if prefix:
            raw_barcodes = pd.Index([
                bc[len(prefix):] if bc.startswith(prefix) else bc
                for bc in sample_barcodes
            ])
        else:
            raw_barcodes = sample_barcodes.copy()

        log.info(f"  Stripped barcodes : {raw_barcodes[:3].tolist()}")

        # Write whitelist CSV (no header — velocyto expects bare barcodes)
        csv_path = output_dir / f"barcodes_{ident}.csv"
        raw_barcodes.to_frame().to_csv(csv_path, index=False, header=False)
        log.info(f"  Whitelist written : {csv_path}")

        barcode_paths[ident] = {
            "csv_path"         : csv_path,
            "prefix"           : prefix,
            "raw_barcodes"     : raw_barcodes,
            "original_barcodes": sample_barcodes,
            "igf_id"           : igf_id,
        }

    return barcode_paths


#################################################################
# 2. Run velocyto on each BAM
#################################################################

def run_velocyto(
    bam_path: Path,
    gtf_path: Path,
    barcode_csv: Path,
    output_dir: Path,
    repeat_gtf: Path = None,
    n_jobs: int = 8,
) -> Path:
    """
    Run velocyto run on a single BAM file and return the path to the loom output.
    """
    cmd = [
        "velocyto", "run",
        "-b", str(barcode_csv),
        "-o", str(output_dir),
        f"-@{n_jobs}",
    ]

    if repeat_gtf is not None:
        cmd += ["-m", str(repeat_gtf)]

    cmd += [str(bam_path), str(gtf_path)]

    log.info(f"Running velocyto:\n  {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"velocyto failed for {bam_path.name}. "
            f"Check logs above for details."
        )

    # velocyto names the loom after the BAM stem
    loom_path = output_dir / f"{bam_path.stem}.loom"
    if not loom_path.exists():
        # Sometimes velocyto writes into a subdirectory
        candidates = list(output_dir.rglob("*.loom"))
        if not candidates:
            raise FileNotFoundError(
                f"Could not find loom output in {output_dir}. "
                f"velocyto may have failed silently."
            )
        loom_path = candidates[0]

    log.info(f"Loom written : {loom_path}")
    return loom_path


#################################################################
# 3. Parse loom and align to h5ad barcodes
#################################################################

def parse_loom(loom_path: Path, barcode_info: dict) -> dict:
    """
    Read a velocyto loom file and return spliced/unspliced matrices indexed
    by the original (prefixed) barcodes used in the h5ad object.

    Returns
    -------
    dict with keys: spliced, unspliced, genes, matched_original_barcodes,
                    coverage_pct
    """
    log.info(f"\nParsing loom: {loom_path}")

    with loompy.connect(str(loom_path), "r") as ds:
        raw_loom_barcodes = ds.ca["CellID"]
        log.info(f"  Raw loom CellID example: {raw_loom_barcodes[:3].tolist()}")

        # velocyto CellID format depends on whether a whitelist was supplied:
        #   With whitelist (-b flag):  "filename:ACGTACGTACGTACGT-1"   (suffix = -1)
        #   Without whitelist:         "filename:ACGTACGTACGTACGTx"    (suffix = x)
        #
        # In both cases:
        #   Step 1 — split on ':' -> take right side
        #   Step 2 — normalise suffix to '-1' (strip x / keep -1)
        def clean_barcode(bc: str) -> str:
            # Strip filename prefix
            if ":" in bc:
                bc = bc.split(":")[-1]
            # Normalise suffix
            if bc.endswith("x"):
                bc = bc[:-1] + "-1"
            return bc

        cleaned_loom_barcodes = pd.Index([clean_barcode(bc) for bc in raw_loom_barcodes])

        # Log format so any remaining mismatch is immediately visible
        log.info(f"  Cleaned loom barcodes example: {cleaned_loom_barcodes[:3].tolist()}")
        log.info(f"  h5ad raw barcodes example    : {barcode_info['raw_barcodes'][:3].tolist()}")

        spliced   = csr_matrix(ds["spliced"][:, :].T)   # -> cells x genes
        unspliced = csr_matrix(ds["unspliced"][:, :].T)
        genes     = pd.Index(ds.ra["Gene"])

    # Match cleaned loom barcodes to raw_barcodes (stripped of prefix)
    raw_barcodes      = barcode_info["raw_barcodes"]       # stripped from h5ad
    original_barcodes = barcode_info["original_barcodes"]  # as in h5ad object

    common_raw = raw_barcodes.intersection(cleaned_loom_barcodes)
    coverage   = len(common_raw) / len(raw_barcodes) * 100

    log.info(f"  h5ad cells (this sample)  : {len(raw_barcodes)}")
    log.info(f"  Loom cells                : {len(cleaned_loom_barcodes)}")
    log.info(f"  Matched cells             : {len(common_raw)}")
    log.info(f"  Coverage                  : {coverage:.1f}%")

    if coverage < 50:
        log.warning(
            "Coverage is below 50% — likely a barcode format mismatch. "
            "Inspect the examples above and adjust the cleaning logic."
        )

    # Build index arrays for subsetting
    loom_idx  = [cleaned_loom_barcodes.get_loc(bc) for bc in common_raw]

    # Map back to original (prefixed) barcodes
    prefix            = barcode_info["prefix"]
    matched_originals = pd.Index([
        f"{prefix}{bc}" if prefix else bc for bc in common_raw
    ])

    return {
        "spliced"                  : spliced[loom_idx, :],
        "unspliced"                : unspliced[loom_idx, :],
        "genes"                    : genes,
        "matched_original_barcodes": matched_originals,
        "coverage_pct"             : coverage,
    }


#################################################################
# 4. Transfer layers into the h5ad object
#################################################################

def transfer_velocity_layers(
    adata: ad.AnnData,
    loom_data_per_sample: dict,
) -> ad.AnnData:
    """
    Add spliced and unspliced layers to adata from the parsed loom data.
    Cells with no velocity coverage receive zeros (scVelo handles sparsity).
    """
    log.info("\n=== Transferring velocity layers ===")

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # Initialise empty sparse matrices — zeros for unmatched cells
    from scipy.sparse import lil_matrix
    spliced_full   = lil_matrix((n_cells, n_genes), dtype=np.float32)
    unspliced_full = lil_matrix((n_cells, n_genes), dtype=np.float32)

    # Build obs and var lookup indices once
    obs_lookup = pd.Series(np.arange(n_cells), index=adata.obs_names)
    var_lookup = pd.Series(np.arange(n_genes), index=adata.var_names)

    total_matched = 0

    for ident, loom_data in loom_data_per_sample.items():
        log.info(f"\n  {ident}:")
        log.info(f"    Coverage : {loom_data['coverage_pct']:.1f}%")

        matched_bcs   = loom_data["matched_original_barcodes"]
        loom_genes    = loom_data["genes"]
        spliced_mat   = loom_data["spliced"]
        unspliced_mat = loom_data["unspliced"]

        # Intersect genes
        # velocyto loom gene index often contains duplicates (same gene name
        # on multiple strands). We keep only the first occurrence for indexing
        # and restrict common_genes to those that are unique in the loom.
        loom_gene_counts  = loom_genes.value_counts()
        unique_loom_genes = loom_gene_counts[loom_gene_counts == 1].index
        common_genes      = adata.var_names.intersection(unique_loom_genes)
        log.info(f"    Common genes (unique in loom) : {len(common_genes)} / {len(loom_genes)}")

        # Manual mapping: gene_name -> first positional index in loom
        # Sidesteps get_indexer's requirement for a unique Index
        loom_gene_to_idx = {}
        for i, g in enumerate(loom_genes):
            if g not in loom_gene_to_idx:
                loom_gene_to_idx[g] = i

        loom_gene_idx  = np.array([loom_gene_to_idx[g] for g in common_genes])
        adata_gene_idx = np.array(var_lookup[common_genes].values, dtype=int)

        # Drop any -1 entries (genes not found — safety check)
        valid_gene_mask = loom_gene_idx >= 0
        loom_gene_idx   = loom_gene_idx[valid_gene_mask]
        adata_gene_idx  = adata_gene_idx[valid_gene_mask]

        # Intersect cells
        valid_mask           = matched_bcs.isin(adata.obs_names)
        valid_loom_cell_idx  = np.where(valid_mask)[0]
        valid_adata_cell_idx = obs_lookup[matched_bcs[valid_mask]].values

        spliced_sub   = spliced_mat[valid_loom_cell_idx, :][:, loom_gene_idx]
        unspliced_sub = unspliced_mat[valid_loom_cell_idx, :][:, loom_gene_idx]

        # Insert into full matrices
        for i, (adata_i, loom_i) in enumerate(
            zip(valid_adata_cell_idx, range(spliced_sub.shape[0]))
        ):
            spliced_full[adata_i, adata_gene_idx]   = spliced_sub[loom_i, :].toarray()
            unspliced_full[adata_i, adata_gene_idx] = unspliced_sub[loom_i, :].toarray()

        valid_bcs = matched_bcs[valid_mask]
        total_matched += len(valid_bcs)
        log.info(f"    Cells inserted : {len(valid_bcs)}")

    log.info(f"\n  Total cells with velocity data : {total_matched} / {n_cells} "
             f"({total_matched/n_cells*100:.1f}%)")

    adata.layers["spliced"]   = csr_matrix(spliced_full)
    adata.layers["unspliced"] = csr_matrix(unspliced_full)

    # Flag cells with velocity coverage
    has_velocity = np.array(adata.layers["spliced"].sum(axis=1)).flatten() > 0
    adata.obs["has_velocity"] = has_velocity
    log.info(f"  obs['has_velocity'] added ({has_velocity.sum()} True)")

    return adata


#################################################################
# 5. scVelo preprocessing — on raw spliced/unspliced (NOT SCT data)
#################################################################

def run_scvelo(adata: ad.AnnData, output_dir: Path) -> ad.AnnData:
    """
    Run scVelo RNA velocity pipeline.

    Layer selection rationale:
      scVelo's filter_and_normalize and moments operate on adata.X by default.
      We do NOT want to overwrite counts_SCT. Instead we:
        1. Set adata.X temporarily to the raw spliced layer for scVelo normalisation.
        2. Run scVelo preprocessing entirely through layers= arguments.
        3. Keep counts_SCT layer untouched throughout.

    The dynamical model is most accurate; stochastic is faster and still robust.
    """
    log.info("\n=== scVelo preprocessing ===")

    scv.settings.verbosity = 2
    scv.settings.figdir    = str(output_dir / "figures")
    os.makedirs(scv.settings.figdir, exist_ok=True)

    # Monkey-patch scVelo's make_unique_list for pandas >=2.0
    # pandas.unique() rejects plain lists, and wrapping in np.array breaks
    # downstream iteration (0-d arrays are unhashable). This full replacement
    # uses pure Python dict.fromkeys for order-preserving dedup.
    import scvelo.tools.utils as _scv_utils
    import scvelo.tools._em_model_core as _scv_em

    def _patched_make_unique(key, allow_array=False):
        if isinstance(key, pd.Index):
            key = key.tolist()
        is_list = isinstance(key, (list, tuple, np.record)) or (
            allow_array and isinstance(key, np.ndarray)
        )
        is_list_of_str = is_list and len(key) > 0 and isinstance(key[0], str)
        if is_list_of_str:
            return list(dict.fromkeys(key))
        elif is_list and len(key) < 20:
            return key
        else:
            return [key]

    # Must patch both the source module AND _em_model_core which already
    # imported the function by name (from ... import creates a local ref)
    _scv_utils.make_unique_list = _patched_make_unique
    _scv_em.make_unique_list = _patched_make_unique
    log.info("  Patched scVelo make_unique_list for pandas >=2.0 compat")

    # Set X to raw spliced counts
    log.info("  Setting adata.X = spliced layer (raw counts for scVelo)")
    adata.X = adata.layers["spliced"].copy()

    # Preprocessing via scanpy (replaces deprecated scVelo wrappers)
    scv.pp.filter_genes(adata, min_shared_counts=20)
    log.info(f"  Genes after filtering : {adata.n_vars}")

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Neighbors via scanpy (scVelo >=0.4.0 requires this explicitly)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    log.info("  PCA + neighbors computed via scanpy")

    # Moments (now uses the precomputed neighbor graph)
    scv.pp.moments(adata)

    # Velocity estimation (dynamical model)
    log.info("  Running dynamical model (this may take several minutes)...")
    scv.tl.recover_dynamics(adata, n_jobs=4)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

    # Velocity confidence and latent time
    scv.tl.velocity_confidence(adata)
    scv.tl.latent_time(adata)

    log.info("  Velocity computation complete.")
    log.info(f"  New obs columns : {[c for c in adata.obs.columns if 'velocity' in c or 'latent' in c]}")
    log.info(f"  New layers      : {list(adata.layers.keys())}")

    # Restore SCT layer as adata.X
    if "counts_SCT" in adata.layers:
        log.info("  Restoring adata.X = counts_SCT")
        adata.X = adata.layers["counts_SCT"].copy()
    else:
        log.warning(
            "  counts_SCT layer not found — adata.X left as scVelo-normalised spliced."
        )

    return adata


#################################################################
# 6. Diagnostic summary
#################################################################

def print_diagnostic_summary(adata: ad.AnnData, output_dir: Path):
    log.info("\n=== Final Object Summary ===")
    log.info(f"  Shape        : {adata.shape}")
    log.info(f"  Layers       : {list(adata.layers.keys())}")
    log.info(f"  Obs columns  : {adata.obs.columns.tolist()}")

    for ident in SAMPLE_MAP:
        mask = adata.obs["orig.ident"] == ident
        has_v = adata.obs.loc[mask, "has_velocity"].sum()
        total = mask.sum()
        log.info(f"  {ident:10s} : {has_v}/{total} cells have velocity "
                 f"({has_v/total*100:.1f}%)")

    # Quick coverage plot
    fig, ax = plt.subplots(figsize=(6, 4))
    coverage_vals = []
    labels = []
    for ident in SAMPLE_MAP:
        mask = adata.obs["orig.ident"] == ident
        pct  = adata.obs.loc[mask, "has_velocity"].mean() * 100
        coverage_vals.append(pct)
        labels.append(ident)
    ax.bar(labels, coverage_vals, color=["#4C72B0", "#DD8452", "#55A868"])
    ax.set_ylabel("% cells with velocity data")
    ax.set_title("Velocity coverage per sample")
    ax.set_ylim(0, 105)
    for i, v in enumerate(coverage_vals):
        ax.text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=10)
    plt.tight_layout()
    fig_path = output_dir / "velocity_coverage.png"
    fig.savefig(fig_path, dpi=150)
    log.info(f"  Coverage plot saved : {fig_path}")
    plt.close()


#################################################################
# CLI
#################################################################

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad",        required=True,  help="Path to annotated .h5ad object")
    p.add_argument("--bam_dir",     required=True,  help="Directory containing the BAM files")
    p.add_argument("--gtf",         required=True,  help="Genome GTF (same as used for Cell Ranger)")
    p.add_argument("--output_dir",  default="./velocity_output", help="Output directory")
    p.add_argument("--repeat_gtf",  default=None,   help="Repeat masker GTF (optional but recommended)")
    p.add_argument("--n_jobs",      default=8, type=int, help="Threads for velocyto and scVelo")
    p.add_argument("--skip_velocyto", action="store_true",
                   help="Skip velocyto run (use if looms already exist in output_dir)")
    return p.parse_args()


#################################################################
# Main
#################################################################

def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_dir    = Path(args.bam_dir)
    gtf_path   = Path(args.gtf)
    repeat_gtf = Path(args.repeat_gtf) if args.repeat_gtf else None

    # Load h5ad
    log.info(f"Loading h5ad: {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    log.info(f"  Shape: {adata.shape}")
    log.info(f"  Layers: {list(adata.layers.keys())}")
    log.info(f"  obs['orig.ident'] values: {adata.obs['orig.ident'].unique().tolist()}")

    # Inspect barcodes and export per-sample whitelists
    barcode_info = inspect_and_export_barcodes(adata, output_dir)

    # Run velocyto (or skip if looms already exist)
    loom_paths = {}

    bam_name_map = {
        "Anagen":  "sample_alignments_IGF133568.bam",
        "Catagen": "sample_alignments_IGF133569.bam",
        "Telogen": "sample_alignments_IGF133570.bam",
    }

    for ident, bam_name in bam_name_map.items():
        bam_path  = bam_dir / bam_name
        loom_stem = bam_path.stem

        # velocyto appends a random 5-char ID to the loom filename on every run
        # so we search by prefix pattern rather than exact name
        existing_looms = sorted(output_dir.glob(f"{loom_stem}*.loom"))
        log.info(f"  Loom search pattern : {loom_stem}*.loom")
        log.info(f"  Looms found         : {[l.name for l in existing_looms]}")

        if args.skip_velocyto:
            if existing_looms:
                loom_path = existing_looms[-1]
                log.info(f"Skipping velocyto for {ident} — using existing loom: {loom_path}")
                loom_paths[ident] = loom_path
                continue
            else:
                log.warning(
                    f"--skip_velocyto set but no loom found matching "
                    f"'{loom_stem}*.loom' in {output_dir}. "
                    f"Falling through to run velocyto."
                )

        if not bam_path.exists():
            raise FileNotFoundError(f"BAM not found: {bam_path}")

        loom_paths[ident] = run_velocyto(
            bam_path    = bam_path,
            gtf_path    = gtf_path,
            barcode_csv = barcode_info[ident]["csv_path"],
            output_dir  = output_dir,
            repeat_gtf  = repeat_gtf,
            n_jobs      = args.n_jobs,
        )

    # Parse looms and align barcodes
    loom_data_per_sample = {}
    for ident, loom_path in loom_paths.items():
        loom_data_per_sample[ident] = parse_loom(loom_path, barcode_info[ident])

    # Transfer velocity layers into h5ad
    adata = transfer_velocity_layers(adata, loom_data_per_sample)

    # Run scVelo
    adata = run_scvelo(adata, output_dir)

    # Diagnostics
    print_diagnostic_summary(adata, output_dir)

    # Save
    out_path = output_dir / "annotated_with_velocity.h5ad"
    adata.write_h5ad(out_path)
    log.info(f"\n✓ Saved: {out_path}")
    log.info(f"  Final layers: {list(adata.layers.keys())}")


if __name__ == "__main__":
    main()