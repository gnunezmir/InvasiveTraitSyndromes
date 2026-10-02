
library(ape)
library(Rphylopars)

#Inputs
myTree <- read.tree("C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\RPhylopars Inputs\\output_tree3.tre")
traitdata=read.csv("C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\Round2\\RPhylopars Inputs\\PreImputationData_transformed2.csv")

#Run imputation
p_BM <- phylopars(trait_data = traitdata,tree = myTree)

#Extract imputed data including 95% CI
imp_data <- as.data.frame(p_BM$anc_recon[1:930,])
imp_data_low95 <- as.data.frame(p_BM$anc_recon[1:930,] - sqrt(p_BM$anc_var[1:930,])*1.96)
imp_data_up95 <- as.data.frame(p_BM$anc_recon[1:930,] + sqrt(p_BM$anc_var[1:930,])*1.96)

#Export
write.csv(imp_data,"C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\Round2\\imputedtraitsdata.csv")
write.csv(imp_data_low95,"C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\Round2\\imputedtraitsdata_lower95.csv")
write.csv(imp_data_up95,"C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Data imputation\\Round2\\imputedtraitsdata_upper95.csv")


