# winratioReg 1.0.0 release checklist

From the repository root in RStudio:

```r
devtools::load_all()
devtools::test()
source("inst/validation/local_ipcw_validation.R")
source("inst/validation/phcopula_validation.R")
source("inst/validation/frank_gaussian_validation.R")
source("inst/validation/degree_selection_validation.R")
source("inst/validation/user_workflow_validation.R")
devtools::check()
citation("winratioReg")
```

The target check result is:

```text
0 errors | 0 warnings | 0 notes
```

Then build the source package:

```r
devtools::build()
```

This should create a source archive named `winratioReg_1.0.0.tar.gz`.

Before pushing, review `git status` or GitHub Desktop to ensure that no raw C9633 data, local output files, checkpoint files, or credentials are staged. Commit the release, push to GitHub, and create a GitHub release/tag `v1.0.0` when ready.

A permanent archive (for example, via Zenodo linked to the GitHub release) can be added later if desired.
