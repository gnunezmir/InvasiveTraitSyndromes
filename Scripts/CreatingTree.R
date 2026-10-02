# load the package
library("U.PhyloMaker")

# input the sample species list, the megatree and genus-family relationship files
sp.list <- read.csv("C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\U.Phylo Inputs\\u.phylo_species.csv")
megatree <-read.tree("C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\U.Phylo Inputs\\plant_megatree.tre")
gen.list <- read.csv("C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\U.Phylo Inputs\\plant_genus_list.csv")

# generate a phylogeny for the sample species list
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
write.tree(result$phylo, "C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\U.Phylo Inputs\\output_tree3.tre")
write.csv(result$sp.list, "C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\U.Phylo Inputs\\output_splist.csv")
