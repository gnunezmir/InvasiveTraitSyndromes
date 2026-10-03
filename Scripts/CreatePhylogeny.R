################################################################################
# Create phylogenetic tree
#
# Purpose:
# Generate the phylogenetic tree used in subsequent phylogenetic generalized
# least squares (PGLS) analyses.
#
# Inputs:
#   u.phylo_species.csv   - Species list for taxa included in the study
#   plant_megatree.tre    - Megatree used as the phylogenetic backbone
#   plant_genus_list.csv  - Genus-family reference file provided by U.PhyloMaker
#
# Outputs:
#   output_tree5.tre      - Phylogenetic tree used in subsequent analyses
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