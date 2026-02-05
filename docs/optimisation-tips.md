---
layout: default
title: Optimisation tips
permalink: /optimisation-tips/
nav_order: 5
---

# Optimisation tips

## Use staged optimisation
Don’t chase 95–99% pass rate on the first run.

A good pattern:
1. **Coarse grid** (larger step) + moderate threshold → find a strong basin
2. **Finer grid** (smaller step) → refine
3. Tighten **threshold** gradually while reducing step size

## Watch what “improving” actually means
- **Pass rate** can rise while **tail failures** (catastrophic errors) remain.
- Track **P95** and **Max** alongside pass rate.

## When you get stuck
- Increase probes per mic a little
- Allow occasional **shake/reset**
- Try small manual mic adjustments to break symmetry
- Re-run WAH before re-optimising after major parameter changes

## Practical geometry advice
- Avoid near-coplanar layouts
- Ensure diversity in 3D baselines (not all mics on one ring)
- Keep minimum spacing realistic for your mount/manufacturing tolerances