# Nunez-Mir et al. — Analysis Code and Data

Code and data underlying "Invasive plants display divergent trait syndromes across invasion outcomes and growth forms."

## Contents

- `Scripts/Analyses_MainText_Script.R` — full analysis script (PGLS
  models reported in the main text and Tables S2-S10)
- `Data/CombinedData_draft_4_6_26.csv` — species trait dataset
- `Data/output_tree5.tre` — phylogenetic tree
- `Data/rownames_8_28_25.csv` — species name order matching the tree tips

## Running the analysis

1. Open `Analyses_MainText_ForReviewers.Rproj` in RStudio.
2. Run `renv::restore()` to install the exact package versions used
   (see `renv.lock`).
3. Open and run `Scripts/Analyses_MainText_Script.R`.

All file paths in the script are relative to the project root, so no
changes are needed regardless of where the project is located.

## Requirements

- R >= 4.4.0 (this project was developed under R 4.5.1). If `renv::restore()`
  fails for `MASS` or `Matrix`, your R version is likely too old. Update R and try again.
