################################################################################
# Create phylogenetic tree
#
# Purpose:
# Generate the phylogenetic tree used in subsequent phylogenetic generalized
# least squares (PGLS) analyses.
#
# Inputs:
#   PhylogenySpeciesList.csv   - Species list for taxa included in the study
#   PlantMegatree.tre          - Megatree used as the phylogenetic backbone
#   PlantGenusList.csv         - Genus-family reference file provided by U.PhyloMaker
#
# Outputs:
#   AnalysisPhylogeny.tre      - Phylogenetic tree used in subsequent analyses
#   
#
# Phylogeny construction:
# Species were placed onto the supplied megatree using U.PhyloMaker with
# nodes.type = 1 and scenario = 3.
################################################################################

library("U.PhyloMaker")

# Read input files
# Read input files
sp.list <- read.csv("Data/PhylogenySpeciesList.csv")
megatree <- read.tree("Data/PlantMegatree.tre")
gen.list <- read.csv("Data/PlantGenusList.csv")

# Generate phylogenetic tree
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)

# Export phylogenetic tree
write.tree(result$phylo,"Data/AnalysisPhylogeny.tre")