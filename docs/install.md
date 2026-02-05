---
layout: default
title: Install
permalink: /install/
nav_order: 2
---

# Install

WAHi-Optimiser can be used in two ways:

- **Packaged installer (recommended)** — no full MATLAB licence required (uses MATLAB Runtime).
- **From MATLAB (development)** — best if you plan to edit the code or contribute.

> **Note**
>
> For the app to look its best, [install the Lato font](https://fonts.google.com/specimen/Lato) on your computer.  If you do not have the recommended font, it defaults to the system fonts after checking for a few fallback options.

---

## Option A — Use the packaged installer (recommended)

1. Download the latest release from the repository’s **Releases** page.
2. Run the installer for your platform.
3. Install **MATLAB Runtime** when prompted (the installer will guide you).

<p style="margin-top: 1rem;">
  <a class="btn" href="https://github.com/raviumadi/WAHi-Optimiser/releases">Go to Releases</a>
  <a class="btn" href="https://www.mathworks.com/products/compiler/matlab-runtime.html">MATLAB Runtime (MathWorks)</a>
</p>
<figure style="margin: 1.25rem 0;">
  <img src="{{ 'assets/img/install/mac_installer.png' | relative_url }}"
       alt="WAHi-Optimiser installed on macOS"
       style="width:100%; max-width:900px; height:auto; border-radius:10px; box-shadow: 0 6px 18px rgba(0,0,0,0.08);">
  <figcaption style="margin-top:10px; font-size:0.95em; color:#666;">
    <strong>WAH<em>i</em>:</strong> Installer.
  </figcaption>
</figure>

### Verify installation

After installation, you should be able to launch WAHi-Optimiser from your Applications/Start Menu (depending on platform).

<figure style="margin: 1.25rem 0;">
  <img src="{{ '/assets/img/install/mac_app_start.png' | relative_url }}"
       alt="WAHi-Optimiser version check placeholder"
       style="width:100%; max-width:900px; height:auto; border-radius:10px; box-shadow: 0 6px 18px rgba(0,0,0,0.08);">
  <figcaption style="margin-top:10px; font-size:0.95em; color:#666;">
    <strong>macOS Version:</strong> Installed app.
  </figcaption>
</figure>


---

## Option B — Run from MATLAB (development)

Use this route if you want to modify the code, run experiments, or contribute improvements.

### 1) Clone the repository

```bash
git clone https://github.com/raviumadi/WAHi-Optimiser.git
cd WAHi-Optimiser
```

### 2) Add WAHi-Optimiser to your MATLAB path
In MATLAB (from the repo root):

```matlab
addpath(genpath(pwd));
savepath;
```

### 3) Launch the GUI
```matlab
wahi_gui
```
<figure style="margin: 1.25rem 0;">
  <img src="{{ '/assets/img/install/matlab_gui_start.png' | relative_url }}"
       alt="MATLAB launch placeholder"
       style="width:100%; max-width:900px; height:auto; border-radius:10px; box-shadow: 0 6px 18px rgba(0,0,0,0.08);">
  <figcaption style="margin-top:10px; font-size:0.95em; color:#666;">
    <strong>GUI Launch:</strong> MATLAB `wahi_gui` launch and successful first Run WAH with 6 mics.
  </figcaption>
</figure>


## Dependencies

WAHi-Optimiser is designed to run on standard MATLAB installations, but some workflows require specific toolboxes.
	•	Optimisation Toolbox (for optimisation runs using lsqnonlin).
MathWorks link: https://www.mathworks.com/products/optimization.html
	•	Localisation engine class (e.g., BatCallLocaliser.m) must be on the MATLAB path.
This is included in the repo / app bundle, but if you are integrating your own engine, ensure it’s discoverable via addpath.

Quick dependency check (MATLAB)

```matlab
ver % lists installed toolboxes
which lsqnonlin -all
which BatCallLocaliser -all
```
<figure style="margin: 1.25rem 0;">
  <img src="{{ '/assets/img/install/requirement_check.png' | relative_url }}"
       alt="Toolbox check placeholder"
       style="width:100%; max-width:900px; height:auto; border-radius:10px; box-shadow: 0 6px 18px rgba(0,0,0,0.08);">
  <figcaption style="margin-top:10px; font-size:0.95em; color:#666;">
    <strong>Check Requirements:</strong> `ver` output and `which` checks showing dependencies found. You will not need to replicate the shown list, but ensure all necessary toolboxes and functions are available.
  </figcaption>
</figure>

## GUI-Limited Features

The `Export Workspace` option is only enabled when running the GUI via MATLAB. This allows easy export of the data to the base workspace and supports on-the-fly analyses. In the deployed version, this is disabled. 

## Troubleshooting
1. I get “undefined function” errors
	Make sure you ran addpath(genpath(pwd)); from the repo root, and restart MATLAB once.
2. Optimisation features are disabled / lsqnonlin not found
	Install the Optimisation Toolbox (or run in a MATLAB environment where it is available).
3. Runtime installed but the app won’t start
	Confirm your Runtime version matches the release requirements (2023b).