# Fisher-KPP Tip/Stalk Reaction-Diffusion Model

MATLAB code for fitting and analyzing a Fisher-type reaction-diffusion
model of tip/stalk/collision growth dynamics against image-derived
density data.

## Folders

- `Optimization/` - fits the 8 model parameters to a target image
  (`s39.png`) via constrained optimization. Entry point: `test.m`.
- `Computer_model/` - a standalone variant with configurable initial
  conditions (head-on collision, cuts, synthetic T-shapes) and a
  mid-run domain "restart" (bisection) feature. Entry point:
  `fisher4_restart.m`.
- `Sensitivity_analysis/` - one-parameter-at-a-time sensitivity sweep
  of the 8 model parameters around a reference set. Entry point:
  `sensitivity1.m`.

Each folder is self-contained: it has its own copy of the input image
(`s3.png` / `s39.png`) and every `.m` file it depends on, so it can be
run independently of the others.

## Requirements

- MATLAB with the Optimization Toolbox (`Optimization/` uses `fmincon`)
  and the Statistics and Machine Learning Toolbox (`Computer_model/`
  uses `ksdensity`).
