---
layout: default
title: Install
permalink: /install/
nav_order: 2
---

# Install

## Option A — Use the packaged installer (recommended)
1. Download the latest release from the repository’s **Releases** page.
2. Run the installer for your platform.
3. Install MATLAB Runtime when prompted (the installer will guide you).

## Option B — Run from MATLAB (development)
1. Clone the repo:
   ```bash
   git clone https://github.com/raviumadi/WAHi-Optimiser.git

2. In MATLAB, add the repo to your path:
```matlab
 addpath(genpath(pwd));
savepath;
````
3. Launch the GUI:
   ```matlab
   wahi_gui
   ```

## Dependencies
	•	Optimisation Toolbox (lsqnonlin) for optimisation runs
	•	The localisation engine class (e.g., BatCallLocaliser.m) must be on the MATLAB path (included in the repo/app bundle)
