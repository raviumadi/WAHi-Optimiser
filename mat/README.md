# MATLAB Analysis Scripts (`mat/`)

This folder contains MATLAB scripts used for batch optimisation runs, post-hoc analyses, and figure generation associated with the **WAH-*i*** framework.

The scripts here are intended for:
- Reproducing analyses reported in the accompanying paper
- Running batch optimisation experiments
- Analysing and visualising results exported from the WAH-*i* GUI

---

## Quick Start (GUI users)

If you are primarily using the **WAH-*i* GUI**, the most relevant script is:

### `iwah_singleRun_analyses.m`

Use this script to analyse a **single optimisation run exported from the GUI** (`.mat` file).

Typical use case:
1. Run an optimisation in the GUI
2. Export the run (`.mat`) via the *Export Run* functionality
3. Set the `RUN_FILE` path in `iwah_singleRun_analyses.m`
4. Run the script to generate:
   - Before/after accuracy distributions
   - Distance-dependent error plots
   - Pass-rate trajectories
   - Geometry drift diagnostics

This script is ideal for inspecting and validating individual designs.

---

## Script Overview

### Core experiment control
- **`run_iwah_experiment_batch.m`**  
  
  Runs large-scale batch experiments across array geometries, microphone counts, and random initialisations.  
  
  Used to generate the datasets stored in `data/runs/`.

### Batch-level analyses
- **`batch_run_analyses.m`**  
  
  Aggregates batch results, computes summary statistics, and prepares data for downstream visualisation and tables.
  
- **`economy_gain_table.m`**  
  
  Computes economy-of-design metrics (performance gain vs. microphone count and geometry).

### Geometry and stability diagnostics
- **`geometry_drift_figure.m`**  
  
  Analyses and visualises microphone position drift across optimisation iterations.

### Single-run analysis (GUI-export friendly)
- **`iwah_singleRun_analyses.m`**  
  Detailed analysis of a single optimisation run (see *Quick Start* above).

---

## Outputs

These scripts generate figures and tables that are typically written to:
- `fig/` – publication-quality figures
- `results/` – CSV and LaTeX-ready tables
- `reports/` – run-specific diagnostic plots

Exact output paths are defined within each script and can be customised.

---

## License

These scripts are released under **GPLv3**, consistent with the main WAH-*i* software license.