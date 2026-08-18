library(MASS)
library(ggplot2)
library(broom.mixed)
library(dplyr)
library(nlme)
library(ape)
library(purrr)
library(performance)
library(ggpubr)



data=read.csv("Data\\CombinedData_draft_4.6.26.csv",header=T,stringsAsFactors = T)


###############################################################################
####################### GROUPING TRAITS BY MECHANISM ##########################
###############################################################################

#################### Transforming traits for normality ########################

b <- boxcox(lm(data$InvasiveRangeSize ~ 1))
lambda <- b$x[which.max(b$y)]
data$InvasiveRangeSize_t <- (data$InvasiveRangeSize ^ lambda - 1) / lambda

b <- boxcox(lm(data$HabitatBreadth ~ 1))
lambda <- b$x[which.max(b$y)]
data$HabitatBreadth_t <- (data$HabitatBreadth ^ lambda - 1) / lambda

b <- boxcox(lm((data$LocalAbundance+0.5) ~ 1))
lambda <- b$x[which.max(b$y)]
data$LocalAbundance_t <- ((data$LocalAbundance+0.5) ^ lambda - 1) / lambda

b <- boxcox(lm((data$IntroYear) ~ 1))
lambda <- b$x[which.max(b$y)]
data$IntroYear_t <- ((data$IntroYear) ^ lambda - 1) / lambda

b <- boxcox(lm((data$ResidenceTime) ~ 1))
lambda <- b$x[which.max(b$y)]
data$ResidenceTime_t <- ((data$ResidenceTime) ^ lambda - 1) / lambda


b <- boxcox(lm((data$SLA_m2kg1) ~ 1))
lambda <- b$x[which.max(b$y)]
data$SLA_t <- ((data$SLA_m2kg1) ^ lambda - 1) / lambda

b <- boxcox(lm((data$SeedWeightMean_g) ~ 1))
lambda <- b$x[which.max(b$y)]
data$SeedWeightMean_g_t <- ((data$SeedWeightMean_g) ^ lambda - 1) / lambda

b <- boxcox(lm(data$ElevationMax_m ~ 1))
lambda <- b$x[which.max(b$y)]
data$ElevationMax_m_t = (data$ElevationMax_m ^ lambda - 1) / lambda
hist(data$ElevationMax_m_t)

b <- boxcox(lm((data$PrecipMin_in+0.5) ~ 1))
lambda <- b$x[which.max(b$y)]
data$PrecipMin_in_t = ((data$PrecipMin_in+0.5) ^ lambda - 1) / lambda
hist(data$PrecipMin_in_t)

b <- boxcox(lm(data$PrecipMax_in ~ 1))
lambda <- b$x[which.max(b$y)]
data$PrecipMax_in_t = (data$PrecipMax_in ^ lambda - 1) / lambda
hist(data$PrecipMax_in_t)

b <- boxcox(lm(data$HeightMax_m ~ 1))
lambda <- b$x[which.max(b$y)]
data$HeightMax_m_t = (data$HeightMax_m ^ lambda - 1) / lambda
hist(data$HeightMax_m_t)

b <- boxcox(lm(data$LeafWidthMax_cm ~ 1))
lambda <- b$x[which.max(b$y)]
data$LeafWidthMax_cm_t = (data$LeafWidthMax_cm ^ lambda - 1) / lambda
hist(data$LeafWidthMax_cm_t)

b <- boxcox(lm(data$LeafLengthMax_cm ~ 1))
lambda <- b$x[which.max(b$y)]
data$LeafLengthMax_cm_t = (data$LeafLengthMax_cm ^ lambda - 1) / lambda
hist(data$LeafLengthMax_cm_t)

b <- boxcox(lm(data$holoploid ~ 1))
lambda <- b$x[which.max(b$y)]
data$holoploid_t = (data$holoploid ^ lambda - 1) / lambda
hist(data$holoploid_t)

b <- boxcox(lm(data$monoploid ~ 1))
lambda <- b$x[which.max(b$y)]
data$monoploid_t = (data$monoploid ^ lambda - 1) / lambda
hist(data$monoploid_t)

#################### Principal Component Analyses #############################


##------------------ Specific environmental tolerance (SET) --------------------------##
data_SET=na.omit(cbind.data.frame(ScientificName=data$ScientificName,data$ElevationMax_m_t,data$PrecipMin_in_t,
                                  data$PrecipMax_in_t,data$pHMin,data$pHMax,data$HardinessZoneMin,data$HardinessZoneMax))

pca_SET=prcomp(data_SET[-1],center=T,scale.=T)
summary(pca_SET)

###figure out what the influential components represent through barplots
loadings=as.data.frame(pca_SET$rotation)
p <-ggplot(loadings, aes(row.names(loadings),PC3))
p +geom_bar(stat = "identity")+
  coord_flip()

###obtain scores for all principal components
scores_SET=cbind.data.frame(data_SET$ScientificName, pca_SET$x)
colnames(scores_SET)=c("ScientificName","SETPC1","SETPC2","SETPC3")

data = merge(data,scores_SET[1:4],by="ScientificName",all.x=T)


##---------------------------- Superior competitive ability (AGR) -------------------------------##
data_AGR=na.omit(cbind.data.frame(ScientificName=data$ScientificName,data$LeafWidthMax_cm_t,
                                    data$HeightMax_m_t,data$LeafLengthMax_cm_t,data$SLA_t,data$SeedWeightMean_g_t))
                 
pca_AGR=prcomp(data_AGR[-1],center=T,scale.=T)
summary(pca_AGR)

###figure out what the influential components represent through barplots
loadings=as.data.frame(pca_AGR$rotation)
p <-ggplot(loadings, aes(row.names(loadings),PC1))
p +geom_bar(stat = "identity")+
  coord_flip()

###obtain scores for all principal components
scores_AGR=cbind.data.frame(data_AGR$ScientificName, pca_AGR$x)
colnames(scores_AGR)=c("ScientificName","AGRPC1","AGRPC2")


data = merge(data,scores_AGR[1:3],by="ScientificName",all.x=T)



############################## PGLS Models ##################################

myTree <- read.tree("Data\\output_tree5.tre")
rnames <- read.csv("Data\\rownames.8.28.25.csv",header=T)
data$species <-rnames$species

corStruct <- corPagel(1, phy = myTree, fixed = FALSE,form = ~species)

## re-orienting PCs for interpretability
data$SETPC3 <- data$SETPC3*-1
data$AGRPC1 <- data$AGRPC1*-1


## modifying local abundance to capture probability of high abundance per established hexagon
data$AbundanceRate <- log1p(data$LocalAbundance) -  log1p(data$InvasiveRangeSize)

##------------------------- Invasive Range Size -----------------------------##

model1=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=data,correlation = corStruct,na.action = na.omit)

summary(model1)
performance::r2(model1)


##--------------------------- Habitat Breadth -------------------------------##

model2=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
               +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=data,correlation = corStruct,na.action = na.omit)
summary(model2)
performance::r2(model2)


##--------------------------- Local Abundance -------------------------------##

model3=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
               +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=data,correlation = corStruct,na.action = na.omit)
summary(model3)
performance::r2(model3)


########################## Herbaceous #######################################

herbdata=subset(data,TissueType=="Herbaceous")

##------------------------- Invasive Range Size -----------------------------##

model1h=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=herbdata,correlation = corStruct,na.action = na.omit)

summary(model1h)
performance::r2(model1h)



##--------------------------- Habitat Breadth -------------------------------##

model2h=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=herbdata,correlation = corStruct,na.action = na.omit)
summary(model2h)
performance::r2(model2h)



##--------------------------- Local Abundance -------------------------------##


model3h=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=herbdata,correlation = corStruct,na.action = na.omit)
summary(model3h)
performance::r2(model3h)


########################## Woody #######################################

woodydata=subset(data,TissueType=="Woody")


##------------------------- Invasive Range Size -----------------------------##

model1w=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=woodydata,correlation = corStruct,na.action = na.omit)

summary(model1w)
performance::r2(model1w)


##--------------------------- Habitat Breadth -------------------------------##

model2w=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=woodydata,correlation = corStruct,na.action = na.omit)
summary(model2w)
performance::r2(model2w)


##--------------------------- Local Abundance -------------------------------##


model3w=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=woodydata,correlation = corStruct,na.action = na.omit)
summary(model3w)
performance::r2(model3w)



##################### Figures ###########################

##### Main results dot-n-whiskerplots ########
# combine tidied models
df <- bind_rows(
  map_df(list(All=model1,Herbaceous=model1h,Woody=model1w),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Invasive range size"),
  map_df(list(All=model2,Herbaceous=model2h,Woody=model2w),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Habitat breadth"),
  map_df(list(All=model3,Herbaceous=model3h,Woody=model3w),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Local abundance")
) %>%
  filter(term!="(Intercept)") %>%
  mutate(term=factor(term,
                     levels=c("scale(ResidenceTime_t)","scale(AGRPC2)","scale(AGRPC1)",
                              "scale(SETPC3)","scale(SETPC2)","scale(SETPC1)",
                              "RegenerativeCapacity",
                              "scale(natrng)","scale(monoploid_t)"),
                     labels=c("Residence time","Acquisitive-conservative gradient","Pioneer-competitor gradient",
                              "Montane affinity","Soil moisture/pH niche","Warm-adaptation",
                              "Regenerative capacity",
                              "Native range size","Genome size"))
  )%>% 
  mutate(sig = p.value <= 0.05)

# plot with facets
pos <- position_dodge(width = 0.4)

A=ggplot(df, aes(x=estimate, y=term, color=group, alpha = sig)) +
  geom_vline(xintercept=0, colour="grey60", linetype=2) +
  geom_point(size=2, position=pos) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0,position=pos) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2), guide = "none") +
  facet_wrap(~ response, scales="fixed") +
  scale_colour_manual(values=c(All="black",Herbaceous="forestgreen",Woody="saddlebrown")) +
  labs(x="Coefficient Estimate", y=NULL) +
  scale_x_continuous(
    breaks = scales::breaks_width(1))+
  theme_bw(base_size=14) +
  theme(
    legend.title=element_blank(),
    legend.position = "bottom"
  )



############################## PGLS Models by Growth Form ##################################

## modifying local abundance to capture probability of high abundance per established hexagon
data$AbundanceRate <- log1p(data$LocalAbundance) -  log1p(data$InvasiveRangeSize)

########################## Graminoid #######################################

gramdata=subset(data,GrowthForm4Model=="Graminoid")


##------------------------- Invasive Range Size -----------------------------##

model1g=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=gramdata,correlation = corStruct,na.action = na.omit)

summary(model1g)
performance::r2(model1g)


##--------------------------- Habitat Breadth -------------------------------##

model2g=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=gramdata,correlation = corStruct,na.action = na.omit)
summary(model2g)
performance::r2(model2g)


##--------------------------- Local Abundance -------------------------------##

model3g=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
           +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=gramdata,correlation = corStruct,na.action = na.omit)
summary(model3g)
performance::r2(model3g)

########################## Forb #######################################

forbdata=subset(data,GrowthForm4Model=="Forb")


##------------------------- Invasive Range Size -----------------------------##

model1f=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=forbdata,correlation = corStruct,na.action = na.omit)

summary(model1f)
performance::r2(model1f)


##--------------------------- Habitat Breadth -------------------------------##

model2f=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=forbdata,correlation = corStruct,na.action = na.omit)
summary(model2f)
performance::r2(model2f)


##--------------------------- Local Abundance -------------------------------##


model3f=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=forbdata,correlation = corStruct,na.action = na.omit)
summary(model3f)
performance::r2(model3f)


########################## Trees & Shrubs #######################################

tsdata=subset(data,GrowthForm4Model=="Tree"|GrowthForm4Model=="Shrub")


##------------------------- Invasive Range Size -----------------------------##

model1ts=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=tsdata,correlation = corStruct,na.action = na.omit)

summary(model1ts)
performance::r2(model1ts)


##--------------------------- Habitat Breadth -------------------------------##

model2ts=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=tsdata,correlation = corStruct,na.action = na.omit)
summary(model2ts)
performance::r2(model2ts)


##--------------------------- Local Abundance -------------------------------##


model3ts=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=tsdata,correlation = corStruct,na.action = na.omit)
summary(model3ts)
performance::r2(model3ts)

########################## Vine #######################################

vinedata=subset(data,GrowthForm4Model=="Vine")

##------------------------- Invasive Range Size -----------------------------##

model1v=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
             +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=vinedata,correlation = corStruct,na.action = na.omit)

summary(model1v)
performance::r2(model1v)


##--------------------------- Habitat Breadth -------------------------------##

model2v=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
             +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=vinedata,correlation = corStruct,na.action = na.omit)
summary(model2v)
performance::r2(model2v)


##--------------------------- Local Abundance -------------------------------##


model3v=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
             +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=vinedata,correlation = corStruct,na.action = na.omit)
summary(model3v)
performance::r2(model3v)


##### Growth Form dot-n-whiskerplots ########

# combine tidied models
df <- bind_rows(
  map_df(list(Forbs=model1f,Graminoids=model1g,TreesShrubs=model1ts,Vines=model1v),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Invasive range size"),
  map_df(list(Forbs=model2f,Graminoids=model2g,TreesShrubs=model2ts,Vines=model2v),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Habitat breadth"),
  map_df(list(Forbs=model3f,Graminoids=model3g,TreesShrubs=model3ts,Vines=model3v),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Local abundance")
) %>%
  filter(term!="(Intercept)") %>%
  mutate(term=factor(term,
                     levels=c("scale(ResidenceTime_t)","scale(AGRPC2)","scale(AGRPC1)",
                              "scale(SETPC3)","scale(SETPC2)","scale(SETPC1)",
                              "RegenerativeCapacity",
                              "scale(natrng)","scale(monoploid_t)"),
                     labels=c("Residence time","Acquisitive-conservative gradient","Pioneer-competitor gradient",
                              "Montane affinity","Soil moisture/pH niche","Warm-adaptation",
                              "Regenerative capacity",
                              "Native range size","Genome size"))
  )%>% 
  mutate(sig = p.value <= 0.05)

# plot with facets
pos <- position_dodge(width = 0.4)

B=ggplot(df, aes(x=estimate, y=term, color=group, alpha = sig)) +
  geom_vline(xintercept=0, colour="grey60", linetype=2) +
  geom_point(size=2, position=pos) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0,position=pos) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2), guide = "none") +
  facet_wrap(~ response, scales="fixed") +
  scale_colour_manual(values=c(Forbs="darkolivegreen",Graminoids="darkseagreen",TreesShrubs="saddlebrown",Vines="mediumorchid4")) +
  labs(x="Coefficient Estimate", y=NULL) +
  scale_x_continuous(
    breaks = scales::breaks_width(1))+
  theme_bw(base_size=14) +
  theme(
    legend.title=element_blank(),
    legend.position = "bottom"
  )


############################## PGLS Models by Lifespan ##################################

## modifying local abundance to capture probability of high abundance per established hexagon
data$AbundanceRate <- log1p(data$LocalAbundance) -  log1p(data$InvasiveRangeSize)

########################## Annuals #######################################

annualdata=subset(data,Lifespan=="Annual"&TissueType=="Herbaceous")

corStruct <- corPagel(0.5, phy = myTree, fixed = FALSE,form = ~species)

##------------------------- Invasive Range Size -----------------------------##

model1a=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=annualdata,correlation = corStruct,na.action = na.omit)

summary(model1a)
performance::r2(model1a)


##--------------------------- Habitat Breadth -------------------------------##

model2a=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=annualdata,correlation = corStruct,na.action = na.omit)
summary(model2a)
performance::r2(model2a)


##--------------------------- Local Abundance -------------------------------##

model3a=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=annualdata,correlation = corStruct,na.action = na.omit)
summary(model3a)
performance::r2(model3a)


########################## Perennials #######################################

pdata=subset(data,Lifespan=="Perennial"&TissueType=="Herbaceous")

corStruct <- corPagel(1, phy = myTree, fixed = FALSE,form = ~species)

##------------------------- Invasive Range Size -----------------------------##

model1p=gls(InvasiveRangeSize_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=pdata,correlation = corStruct,na.action = na.omit)

summary(model1p)
performance::r2(model1p)


##--------------------------- Habitat Breadth -------------------------------##

model2p=gls(HabitatBreadth_t~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=pdata,correlation = corStruct,na.action = na.omit)
summary(model2p)
performance::r2(model2p)


##--------------------------- Local Abundance -------------------------------##

model3p=gls(AbundanceRate~scale(AGRPC1)+scale(AGRPC2)+scale(SETPC1)+scale(SETPC2)+scale(SETPC3)+scale(ResidenceTime_t)
            +RegenerativeCapacity+scale(natrng)+scale(monoploid_t),data=pdata,correlation = corStruct,na.action = na.omit)
summary(model3p)
performance::r2(model3p)



##### Growth Form dot-n-whiskerplots ########

# combine tidied models
df <- bind_rows(
  map_df(list(Annuals=model1a,Perennials=model1p),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Invasive range size"),
  map_df(list(Annuals=model2a,Perennials=model2p),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Habitat breadth"),
  map_df(list(Annuals=model3a,Perennials=model3p),
         ~tidy(.x, conf.int=TRUE), .id="group") %>% mutate(response="Local abundance")
) %>%
  filter(term!="(Intercept)") %>%
  mutate(term=factor(term,
                     levels=c("scale(ResidenceTime_t)","scale(AGRPC2)","scale(AGRPC1)",
                              "scale(SETPC3)","scale(SETPC2)","scale(SETPC1)",
                              "RegenerativeCapacity",
                              "scale(natrng)","scale(monoploid_t)"),
                     labels=c("Residence time","Acquisitive-conservative gradient","Pioneer-competitor gradient",
                              "Montane affinity","Soil moisture/pH niche","Warm-adaptation",
                              "Regenerative capacity",
                              "Native range size","Genome size"))
  )%>% 
  mutate(sig = p.value <= 0.05)

# plot with facets
pos <- position_dodge(width = 0.4)

C=ggplot(df, aes(x=estimate, y=term, color=group, alpha = sig)) +
  geom_vline(xintercept=0, colour="grey60", linetype=2) +
  geom_point(size=2, position=pos) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0,position=pos) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2), guide = "none") +
  facet_wrap(~ response, scales="fixed") +
  scale_colour_manual(values=c(Annuals="darkblue",Perennials="firebrick")) +
  labs(x="Coefficient Estimate", y=NULL) +
  scale_x_continuous(
    breaks = scales::breaks_width(1))+
  theme_bw(base_size=14) +
  theme(
    legend.title=element_blank(),
    legend.position = "bottom"
  )

################## ALL PLOTS TOGETHER #########################

jpeg("C:\\Users\\gnm\\Dropbox\\DISTRIBUTION OF EXOTICS\\Manuscript\\mainfig.jpeg", width = 9.5, height = 12, units = 'in', res = 300)

ggarrange(A,B,C,labels="AUTO",ncol=1)

dev.off()
