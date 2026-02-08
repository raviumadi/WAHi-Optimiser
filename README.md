# WAH-*i* -  Microphone Array Geometry Optimisation Algorithm

**Companion repository for the research paper:**

**WAH-*i*: Optimising Microphone Array Geometry for Customised Localisation Accuracy**

[Read the Paper](https://doi.org/10.64898/2026.02.07.704547)

This repository contains the complete software, analyses, datasets, and application developed as part of the WAH-*i* research programme. The project introduces a principled optimisation framework for designing microphone array geometries that achieve a desired spatial field of localisation accuracy.

 

### [Quickstart and Usage Guide](https://raviumadi.github.io/WAHi-Optimiser/quickstart/)



---

## Overview

**WAH-*i*** (Widefield Acoustics Heuristic — inverse, iterative) is an optimisation framework that:

- Directly targets **volumetric localisation accuracy**, rather than heuristic geometric criteria
- Supports **task-specific accuracy thresholds** and pass-rate objectives
- Enables **systematic comparison of array geometries** under realistic constraints
- Scales from **small, portable arrays** to larger multichannel systems

The repository includes:

- A full-featured **graphical user interface**
- **Standalone installable application packages**
- Batch-run pipelines used in the paper
- Analysis scripts, figures, and exported results

See the folders `app` and `app_install`

## Requirements

- MATLAB (e.g. 2023b)
  - Optimization Toolbox
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Communications Toolbox

---

## **Graphical User Interface (GUI)**

The interactive WAH-*i* GUI allows full control over optimisation parameters, constraints, and diagnostics.

📁 **Relevant folders**

- app/ — GUI source and documentation
- app_install/ — standalone, installable application packages

![iwah_gui_done_dense_grid](img/iwah_gui_done_dense_grid.png)

The GUI supports:

- Iterative optimisation with live diagnostics
- Manual geometry refinement and restart strategies
- Export to CSV (CAD workflows) and MAT (further analysis)
- Automated report generation

## Batch Runs and Analyses

Batch experiments and analysis scripts used in the paper are provided for full reproducibility.

📁 **MATLAB scripts**

- mat/ — batch execution, analysis, and figure generation scripts

📁 **Data**

- data/runs/ — batch optimisation runs used in the manuscript
- data/single_runs/ — example GUI-driven optimisation runs

📁 **Outputs**

- fig/ — publication figures
- reports/ — GUI-generated summary reports
- results/ — processed tables and statistics

## Array Geometries

📁 **configs/**

Contains CSV files defining microphone array geometries used throughout the study, including:

- Canonical polyhedral configurations (tetrahedron, octahedron, icosahedron, etc.)
- Random initialisations
- Optimised geometries exported directly from the GUI

These files can be loaded directly into the GUI or used for batch experiments.

## Source Code

📁 **src/**

Core algorithmic and utility components, including:

- Optimisation engine
- Localisation model
- Plotting and export utilities

## License

- **Software:** GNU General Public License v3 (GPLv3)
- **Documentation:** CC BY-NC-ND 4.0

See LICENSE.md for details.

## Disclaimer

This software and associated analyses are provided **as-is**, without warranty of any kind.

No guarantee is made regarding suitability for any particular application or experimental setup.

Users are responsible for validating results in their own deployment contexts.

## Support & Acknowledgement

This work is self-funded and developed as part of open, independent research.

If you find WAH-*i* useful and would like to support further development:

☕ **Buy me a coffee:** https://buymeacoffee.com/raviumadi

You can find my other research and open tools at:

🌍 **https://biosonix.io**

---

Whenever possible, consider supporting or promoting **science education in developing countries**.