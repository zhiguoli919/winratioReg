# C9633 application scripts

The scripts in this folder document the C9633 analyses used during manuscript development. **Participant-level C9633 data are not distributed in this repository.** Authorized access to the trial data is required.

## Data path

The scripts read the data path from the `C9633_DATA` environment variable. In an R session, set it to your authorized local copy before sourcing a script, for example:

```r
Sys.setenv(C9633_DATA = "C:/secure/path/C9633.csv.xlsx")
```

The scripts use `readxl::read_excel()` and expect the variable names used in the original analysis files.

## PH--BB1 helper file

The two PH--BB1 application scripts also source `PH_BB1_recovered_functions.R`. If the working directory is this folder, the default works directly. Otherwise set:

```r
Sys.setenv(
  WINRATIOREG_BB1_HELPERS = "C:/path/to/winratioReg/reproducibility/C9633/PH_BB1_recovered_functions.R"
)
```

## Included analyses

Local IPCW:

- overall age
- age within PS 0 and PS 1
- overall tumor size
- tumor size within PS 0 and PS 1

PH--copula:

- overall age, BB1 working copula
- overall tumor size, BB1 working copula

No PH--copula subgroup analyses are included because they were not part of the finalized manuscript analysis.

The application scripts were retained for exact manuscript provenance. The package functions provide the cleaner user-facing interface for new analyses.
