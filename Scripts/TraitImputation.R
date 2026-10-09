################################################################################
# Phylogenetic trait imputation
#
# Purpose:
# Transform functional trait values, impute missing values using phylogenetic
# information with Rphylopars, and back-transform the resulting values to
# their original measurement scales.
#
# Inputs:
#   Data/AnalysisPhylogeny.tre       - Phylogenetic tree
#   Data/TraitImputationInput.csv    - Untransformed trait dataset
#
# Outputs:
#   Data/ImputedTraitData.csv        - Trait values following phylogenetic
#                                     imputation and back-transformation
#
# Notes:
# Observed trait values are retained, while missing values are estimated
# using phylogenetic information where possible.
#
# Following imputation, imputed values falling outside the observed range of
# the corresponding trait were replaced with missing values (blanks). This
# screening step was performed manually after running this script and is not
# included in the code. The screened trait values were subsequently merged
# with the remaining traits and species-level variables from the compiled
# dataset to produce AnalysisData.csv.
################################################################################

library(ape)
library(Rphylopars)
library(MASS)

# Read input files
myTree <- read.tree("Data/AnalysisPhylogeny.tre")

traitdata <- read.csv("Data/TraitImputationInput.csv")

# Apply Box-Cox transformation (log transform for min elevation) to variables
traitdata$ElevationMin_m_t = log1p(traitdata$ElevationMin_m)

b <- boxcox(lm(traitdata$ElevationMax_m ~ 1))
elevmax_lambda <- b$x[which.max(b$y)]
traitdata$ElevationMax_m_t = (traitdata$ElevationMax_m ^ elevmax_lambda - 1) / elevmax_lambda

b <- boxcox(lm((traitdata$PrecipMin_in+0.5) ~ 1))
precipmin_lambda <- b$x[which.max(b$y)]
traitdata$PrecipMin_in_t = ((traitdata$PrecipMin_in+0.5) ^ precipmin_lambda - 1) / precipmin_lambda

b <- boxcox(lm(traitdata$PrecipMax_in ~ 1))
precipmax_lambda <- b$x[which.max(b$y)]
traitdata$PrecipMax_in_t = (traitdata$PrecipMax_in ^ precipmax_lambda - 1) / precipmax_lambda

b <- boxcox(lm(traitdata$HeightMax_m ~ 1))
heightmax_lambda <- b$x[which.max(b$y)]
traitdata$HeightMax_m_t = (traitdata$HeightMax_m ^ heightmax_lambda - 1) / heightmax_lambda

b <- boxcox(lm(traitdata$LeafWidthMax_cm ~ 1))
leafw_lambda <- b$x[which.max(b$y)]
traitdata$LeafWidthMax_cm_t = (traitdata$LeafWidthMax_cm ^ leafw_lambda - 1) / leafw_lambda

b <- boxcox(lm(traitdata$LeafLengthMax_cm ~ 1))
leafl_lambda <- b$x[which.max(b$y)]
traitdata$LeafLengthMax_cm_t = (traitdata$LeafLengthMax_cm ^ leafl_lambda - 1) / leafl_lambda

b <- boxcox(lm(traitdata$holoploid ~ 1))
holoploid_lambda <- b$x[which.max(b$y)]
traitdata$holoploid_t = (traitdata$holoploid ^ holoploid_lambda - 1) / holoploid_lambda

b <- boxcox(lm(traitdata$monoploid ~ 1))
monoploid_lambda <- b$x[which.max(b$y)]
traitdata$monoploid_t = (traitdata$monoploid ^ monoploid_lambda - 1) / monoploid_lambda

# Perform phylogenetic trait imputation

imputation_data <- traitdata[, c(
  "species",
  "pHMin", "pHMax",
  "HardinessZoneMin", "HardinessZoneMax",
  "ElevationMin_m_t", "ElevationMax_m_t",
  "PrecipMin_in_t", "PrecipMax_in_t",
  "HeightMax_m_t",
  "LeafWidthMax_cm_t", "LeafLengthMax_cm_t",
  "holoploid_t", "monoploid_t"
)]

p_BM <- phylopars(trait_data = imputation_data, tree = myTree)

# Extract imputed trait estimates for the study species
imp_data <- as.data.frame(p_BM$anc_recon[seq_len(nrow(imputation_data)),])

# Back-transform data for merging with analysis dataset

imp_data$ElevationMin_m <- expm1(imp_data$ElevationMin_m_t)

imp_data$ElevationMax_m <- (imp_data$ElevationMax_m_t * elevmax_lambda + 1)^(1/elevmax_lambda)

imp_data$PrecipMin_in <- ((imp_data$PrecipMin_in_t * precipmin_lambda + 1)^(1/precipmin_lambda)) - 0.5

imp_data$PrecipMax_in <- (imp_data$PrecipMax_in_t * precipmax_lambda + 1)^(1/precipmax_lambda)

imp_data$HeightMax_m <- (imp_data$HeightMax_m_t * heightmax_lambda + 1)^(1/heightmax_lambda)

imp_data$LeafWidthMax_cm <- (imp_data$LeafWidthMax_cm_t * leafw_lambda + 1)^(1/leafw_lambda)

imp_data$LeafLengthMax_cm <- (imp_data$LeafLengthMax_cm_t * leafl_lambda + 1)^(1/leafl_lambda)

imp_data$holoploid <- (imp_data$holoploid_t * holoploid_lambda + 1)^(1/holoploid_lambda)

imp_data$monoploid <- (imp_data$monoploid_t * monoploid_lambda + 1)^(1/monoploid_lambda)

# Export imputed and back-transformed trait estimates
write.csv(imp_data[, c(
  "species",
  "pHMin", "pHMax",
  "HardinessZoneMin", "HardinessZoneMax",
  "ElevationMin_m", "ElevationMax_m",
  "PrecipMin_in", "PrecipMax_in",
  "HeightMax_m",
  "LeafWidthMax_cm", "LeafLengthMax_cm",
  "holoploid", "monoploid"
)],"Data/ImputedTraitData.csv",row.names = F)

