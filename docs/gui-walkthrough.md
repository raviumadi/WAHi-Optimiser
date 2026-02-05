---
title: GUI walkthrough
nav_order: 4
---

# GUI walkthrough

## Panels overview

### Array & Field Parameters
Defines the geometry constraints and the evaluation region:
- **Lmax**: maximum mic radius (arm length limit)
- **Min mic spacing**: prevents clustering
- **FoA (R, Rin, step)**: evaluation shell
- **Threshold (Target + tolerance)**: defines pass/fail criterion

### Optimiser Parameters
Controls how WAH-*i* explores geometry space:
- **Iterations**
- **Probes per mic**
- **Step size schedule**
- **Stall + shake/reset behaviour**

### Call Synthesis Parameters
Sets the synthetic signal used for TDOA estimation:
- sampling rate, duration, sweep band, SNR, tail padding, call type (FM/CF)

### Manual Management
Manually move microphones (sliders/edit boxes), re-centre, and apply changes.

### Display & Handling
- Metric selection: **Position error** vs **Angular error**
- Octant selection for optimisation region

## Suggested workflow
1. Random geometry → Apply
2. Coarse FoA evaluation
3. Coarse optimisation
4. Tighten grid / tighten threshold
5. Re-optimise with smaller steps