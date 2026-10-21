# Manuscript reproducibility materials

This directory contains analysis scripts retained separately from the `winratioReg` package source. They document the simulation and C9633 workflows used during manuscript development; they are excluded from the built R package via `.Rbuildignore`.

## Simulation scripts

`simulations/` contains the final scripts supplied during package/manuscript development:

- `pointwise_dense_grid_three_class_final.R`: dense-grid three-method pointwise simulation used for figure construction. Its header explicitly notes that it does **not** replace the finalized five-point inferential simulation.
- `local_coefficients_final.R`: final Local IPCW coefficient/BLP simulation.
- `phcopula_coefficients_final.R`: final PH--copula coefficient/BLP simulation.
- `local_loglinearity_1d_final.R` and `local_loglinearity_2d_final.R`: Local IPCW testing simulations.
- `phcopula_loglinearity_1d_final.R`: PH--copula 1D log-linearity simulation.
- `phcopula_loglinearity_2d_alternatives.R`: PH--copula 2D script used for the alternative scenarios.
- `phcopula_loglinearity_2d_matched_H0_final.R`: final H0-only rerun using the 2D near-null DGP matched to the Local IPCW simulation. The manuscript's final 2D PH--copula H0 rejection rates come from this rerun.

These programs are computationally intensive and may use parallel processing, checkpoint files, and long runtimes. Review the settings at the top of each script before running.

## C9633 application

`C9633/` contains the Local IPCW overall/subgroup scripts and the PH--BB1 overall age/tumor-size scripts used during manuscript development. The participant-level CALGB 9633 data are **not included** because data access and redistribution are governed separately.

See `C9633/README.md` for setup instructions.

## Software validation

Development validation scripts are retained under `inst/validation/` in the repository. They are also excluded from the built package. These scripts were used to verify equivalence with the finalized Local IPCW algorithms, the original BB1/Clayton/Gumbel PH--copula implementation, independent Frank/Gaussian copula calculations, Bernstein-degree selection, and the primary user workflow.
