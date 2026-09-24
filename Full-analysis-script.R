# R SCRIPT FOR
#
# "Infection induced changes in gene expression linked to foraging choices 
# in self-medicating Lasius platythorax ants"
#


#### PACKAGES ####
library(readxl)
library(tidyverse)
library(MCMC.qpcr)
library(emmeans)
library(car)
library(ggtext)
library(glmmTMB)
library(DHARMa)
library(survival)
library(coxme)
library(survminer)
library(DESeq2)
library(topGO)
library(cowplot)
library(rempsyc)
library(reshape2)
library(ComplexHeatmap)
library(vegan)

## color palettes ##
palet = c("darkred", "burlywood")
surpalet = c("darkred","lightblue4","orange")
piepalet <- c("red", "purple", "lightblue", "darkblue","gold", "black", "brown", "pink","lavender", "darkgreen")
barpalet <- c("darkblue","gold", "brown", "pink")

#### EXPERIMENT A FORAGING ####

## read in data foraging data for experiment a (foraging) ##
smdf <- read_xlsx("Dataset1.xlsx", sheet = "Foraging")

# re-factor columns
smdfcol <- c("OC", "Box", "Day", "Time", "Food_Type", "Treatment")
smdf[smdfcol] <- lapply(smdf[smdfcol], factor)


### Model building ###

# negative binomial model
ovrlfrg_nb <- glmmTMB(Foragers ~ Food_Type * Treatment + Day + (1|OC) + (1|Time/Box), data = smdf, family = "nbinom1")


# model diagnostics

plot(simulateResiduals(ovrlfrg_nb))

#test temporal autocorrelation
res <- simulateResiduals(ovrlfrg_nb)
res2 <- recalculateResiduals(res, group = smdf$Day)
testTemporalAutocorrelation(res2, time = 1:length(res2$scaledResiduals))

testDispersion(ovrlfrg_nb)
testZeroInflation(ovrlfrg_nb)


#ANOVA of fixed effects
Anova(ovrlfrg_nb)


#Pairwise comparisons
emmeans(ovrlfrg_nb, pairwise ~ Food_Type|Day)
emmeans(ovrlfrg_nb, pairwise ~ Treatment|Food_Type)
emmeans(ovrlfrg_nb, pairwise ~ Food_Type|Treatment)
summary(ovrlfrg_nb)


### PLOTTING ###

smOvrlPltP <- ggplot(smdf, aes(x = Treatment, y = Foragers, fill = Food_Type))+
  geom_bar(position = "fill", stat = "summary", alpha = 1, width = 0.7, colour = "black")+
  scale_fill_manual(values = c("white", "red"), name = "Food type", labels=c("Standard", "ROS"))+
  ggtitle("A) Proportional food choice")+
  labs(y=" Food choice ratio\n", x = "Treatment")+
  scale_x_discrete(labels=expression("Control", italic(B.bassiana)))+
  theme_classic()+
  theme(legend.position = "top")+
  geom_text(x=2, y = 0.23, label = "14.6%")+
  geom_text(x=1, y = 0.23, label = "6.9%")


smOvrlPltP <- ggpar(smOvrlPltP,
                   font.main = 14,
                   font.x = 12,
                   font.y = 12,
                   font.legend = 9)

smline2avg <- ggplot(filter(smdf, Food_Type == 2), aes(x=as.factor(Day), y=Foragers, group=Treatment, color=Treatment))+
  stat_summary(aes(group=Treatment), fun = mean, geom="point", size =2.5, position = position_dodge(0.3))+
  stat_summary(aes(group=Treatment), fun = mean, geom="line", size=0.4, position = position_dodge(0.3))+
  stat_summary(aes(group=Treatment), fun.data = mean_cl_boot, geom = "errorbar", width=0.2, size=0.4, position = position_dodge(0.3))+
  ylab("Average forager frequency \n")+
  ggtitle(expression('B) Foraging on ROS food per day'))+
  xlab("Day")+
  scale_color_manual(values = palet, name = "Treatment", labels=expression("Control", italic(B.bassiana)))+
  theme_classic()+
  theme(legend.position = "top")
smline2avg <- ggpar(smline2avg,
                    font.main = 14,
                    font.x = 12,
                    font.y = 12,
                    font.legend = 9)



#### EXPERIMENT A SURVIVAL ####

# read in data and re-factor #
smSur <- read_xlsx("Dataset1.xlsx", sheet = "Survival")

smSurcol <- c("Colony", "OC", "Treatment", "Diet", "TD")
smSur[smSurcol] <- lapply(smSur[smSurcol], factor)


#Exclude day 1 handling mortality
p_excl <- filter(smSur, Day > 1)

#Coxme model building
cxpl1 <- coxme(Surv(Day, Survival) ~ Treatment * Diet + (1|OC) + (1|Colony), data = p_excl)

# ANOVA of fixed effects
Anova(cxpl1)

#Pairwise comparisons
summary(cxpl1)
emmeans(cxpl1, pairwise~Treatment|Diet)
emmeans(cxpl1, pairwise~Diet|Treatment)
emmeans(cxpl1, pairwise~Diet*Treatment)


## Survival plot ##

sur_excl_c <- survfit(Surv(Day, Survival) ~ TD, data = filter(p_excl, Treatment == "C"))

p_1 <- ggsurvplot(
  sur_excl_c,
  data = p_excl,
  conf.int=TRUE,
  pval=TRUE,
  break.time.by=1,
  palette = surpalet,
  xlim = c(0, 9.4),
  ylim = c(0.7, 1),
  legend.title = "Diet",
  legend = "right",
  legend.labs = c("Standard","Food choice","ROS"),
  ggtheme= theme_survminer(),
  title = expression("C) Control treatment"))


p_1 <- ggpar(p_1,
             font.main = 14,
             font.x = 12,
             font.y = 12,
             font.legend = 9)


sur_excl_f <- survfit(Surv(Day, Survival) ~ TD, data = filter(p_excl, Treatment == "F"))

p_2 <- ggsurvplot(
  sur_excl_f,
  data = p_excl,
  conf.int=TRUE,
  pval=TRUE,
  break.time.by=1,
  palette = surpalet,
  xlim = c(0, 9.4),
  ylim = c(0.7, 1),
  linetype = "dashed",
  legend.title = "Diet",
  legend = "bottom",
  legend.labs = c("Standard","Food choice","ROS"),
  ggtheme= theme_survminer(),
  title = expression("D) Exposed to "~italic('B. bassiana')))


p_2 <- ggpar(p_2,
             font.main = 14,
             font.x = 12,
             font.y = 12,
             font.legend = 9)


fplot <- ggarrange(smOvrlPltP,smline2avg, ncol=2, common.legend = F)
splot <- ggarrange(p_1$plot, p_2$plot +rremove("ylab"), ncol=2, common.legend = T, legend = "bottom")
fig1 <- plot_grid(fplot, splot, nrow = 2)


#### DIFFERENTIAL GENE-EXPRESSION ####

#Read in the AA sequences of the whole genome
lplatSeq <- read_xlsx("Dataset2.xlsx", sheet = "lPlatAAseq")

# Sample files and formatting for day 2
sfiled2 <- read_xlsx("Dataset2.xlsx", sheet = "sample_info_d2")

## Re-classify columns
sfilecol <- c("Treatment", "Time", "Replicate", "SampleName")
sfiled2[sfilecol] <- lapply(sfiled2[sfilecol], factor)

# Sample files and formatting for day 4
sfiled4 <- read_xlsx("Dataset2.xlsx", sheet = "sample_info_d4")
sfiled4[sfilecol] <- lapply(sfiled4[sfilecol], factor)



## Analyze the DEG for day 2 ##

# read in reads
tempdf_2 <- read_xlsx("Dataset2.xlsx", sheet = "reads_d2")

# set first column as row names (genes)
tempdf_2 <- column_to_rownames(tempdf_2, var ="...1")

# save dataframe for NMDS visualization later
nmdvis2 <- tempdf_2

# transform zeros to NA and filter out of dataframe
tempdf_2[tempdf_2==0] <- NA
tempdf_2 <- tempdf_2[complete.cases(tempdf_2),]

# Run DESeq2
dds_2 <- DESeqDataSetFromMatrix(countData = tempdf_2,
                                colData = sfiled2,
                                design=~ Replicate + Treatment)

dds_2 <- DESeq(dds_2)

# Summarize results
resd2 <- results(dds_2, name="Treatment_Infected_vs_Control", test = "Wald", alpha = 0.05)
summary(resd2)

# Store full list of DE genes sorted by adjusted pval
gene_info_2 <- resd2[order(resd2$padj),]

## Picking gene lists

#Picking all significantly differentiated genes from first timepoint
#split into up and downregulated genes

resd2 <- resd2[order(resd2$padj),]
resd2 <- as.data.frame(head(resd2, n = 504))
resupd2 <- filter(resd2, resd2$log2FoldChange > 0)
resdownd2 <- filter(resd2, resd2$log2FoldChange < 0)

# All downregulated genes day 2
fastadownd2 <- lplatSeq[lplatSeq$name %in% row.names(resdownd2),]

# All upregulated genes day 2
fastaupd2 <- lplatSeq[lplatSeq$name %in% row.names(resupd2),]

## Analysis for DEG day 4 ##

#read in reads and remove all rows with 0 values
tempdf_4 <- read_xlsx("Dataset2.xlsx", sheet = "reads_d4")
tempdf_4 <- column_to_rownames(tempdf_4, var ="...1")

# save dataframe for nmds visualization later
nmdvis4 <- tempdf_4

# transform zeros to NA and filter out
tempdf_4[tempdf_4==0] <- NA
tempdf_4 <- tempdf_4[complete.cases(tempdf_4),]



# DESeq2 analysis
dds_d4 <- DESeqDataSetFromMatrix(countData = tempdf_4,
                                 colData = sfiled4,
                                 design=~ Replicate + Treatment)

dds_d4 <- DESeq(dds_d4)

# Analysis levels
resultsNames(dds_d4)

# summarize results
resd4 <- results(dds_d4, name="Treatment_Infected_vs_Control", test = "Wald", alpha = 0.05)
summary(resd4)


# Make summary of the dds results ordered by adjusted p-value
gene_info_4 <- resd4[order(resd4$padj),]


## Picking gene lists ##

#Picking all significantly differentiated genes from first timepoint
#split into up and downregulated genes
resd4 <- resd4[order(resd4$padj),]
resd4 <- as.data.frame(head(resd4, n = 956))
resupd4 <- filter(resd4, resd4$log2FoldChange > 0)
resdownd4 <- filter(resd4, resd4$log2FoldChange < 0)

#All downregulated genes day 4
fastadownd4 <- lplatSeq[lplatSeq$name %in% row.names(resdownd4),]

#All upregulated genes day 4
fastaupd4 <- lplatSeq[lplatSeq$name %in% row.names(resupd4),]



## Shared upregulated genes in day 2 and day 4
dsequp_shared <- as.data.frame(c(row.names(resupd4), row.names(resupd2)))
shared_up <- dsequp_shared[duplicated(dsequp_shared), ]
length(shared_up)


dseqdown_shared <- as.data.frame(c(row.names(resdownd4), row.names(resdownd2)))
shared_down <- dseqdown_shared[duplicated(dseqdown_shared), ]
length(shared_down)


## PLOTTING ##

# Upset plot

# Make a list of the up and downregulated genes on each timepoint
listInput <- list("Downregulated Day 2" = row.names(resdownd2), "Upregulated Day 2" = row.names(resupd2), 
                  "Downregulated Day 4" = row.names(resdownd4), "Upregulated Day 4" = row.names(resupd4))

# Make combination matrix from list
mlst <- make_comb_mat(listInput)

# Draw upset plot
upggA <- UpSet(mlst, top_annotation = upset_top_annotation(mlst, add_numbers = T), right_annotation = upset_right_annotation(mlst, add_numbers = T))

# ggplotify the upsetplot for customization
uppgA <- ggplotify::as.ggplot(upggA)
uppgA <- uppgA + ggtitle(" B) Differentially expressed genes distribution")


# NMDS plot

#pick out unique genes for day 2 and 4 to avoid duplicates
bigdf <- as.data.frame(c(row.names(resd2), row.names(resd4)))
tstuniq <- unique(bigdf)

#choose the raw reads for all of the unique genes during both days

bidfg2 <- nmdvis2[row.names(nmdvis2) %in% tstuniq$`c(row.names(resd2), row.names(resd4))`,]
bidfg4 <- nmdvis4[row.names(nmdvis4) %in% tstuniq$`c(row.names(resd2), row.names(resd4))`,]
bigdfg <- cbind(bidfg2, bidfg4)

# transform and format dataframe
totdf <- as.data.frame(t(bigdfg))
totdf$Sample <- row.names(totdf) 
totdf <- relocate(totdf, Sample)
row.names(totdf) <- NULL

# add treatment variable
totdf$Treatment <- rep(rep(c("Control", "B.bassiana"), each = 10), 2)
totdf <- relocate(totdf, Treatment)

#add timepoint variable
totdf$Day <- rep(c(2, 4), each = 20)
totdf <- relocate(totdf, Day)

# set variables
comtot <- as.data.frame(totdf[,-(1:3)])
envtot <- as.data.frame(totdf[,1:3])

#transform into matrix
m_comtot <- as.matrix(comtot)

# make NMDS
nmdstot <- metaMDS(m_comtot, distance = "bray")

# validate with stressplot
stressplot(nmdstot)

# run model
entot <- envfit(nmdstot, envtot, permutations = 999, ra.rm = T)

td <- rep(c("C2", "F2", "C4", "F4"), each = 10)
treat <- envtot$Treatment
dd <- envtot$Day

#extract nmds site scores
dataScores <- as.data.frame(scores(nmdstot)$sites)
dataScores$Treat = totdf$Treatment
dataScores$comb <- as.factor(td)
dataScores$Day <- as.factor(dd)

# draw ggplot
gg <- ggplot(data = dataScores, aes(x = NMDS1, y = NMDS2))+
  geom_point(data = dataScores, aes(colour = Treat, shape = Day), size = 3, alpha = 0.5)+
  stat_ellipse(aes(fill = td, colour = Treat, geom = "ellipse"))+
  theme_classic()+
  scale_color_manual(values = palet, name = "Treatment", labels=expression("Control", italic(B.bassiana)))+
  ggtitle("A) Ordination plot of differential gene expression")

# combine upset plot and nmds plot
fig2 <- ggarrange(gg, uppgA, nrow = 2)




#### GO ANALYSIS ####

#Read in functional annotation of genes
df <- read.table("combined.byGO.nr.tsv", quote = "")


## Data needs re-formatting ##

#Name of vector = gene name, vector = GO terms
#List all the unique genes in the data file
unique(df$V1)

#Create a vector with all names to loop through
unVec <- unique(df$V1)

#Check it's good
unVec[1]

#Create new list
golist <- list()

#Loop variable
x <- 1

#While loop to extract the information and append to new list

while(x < (length(unVec)+1)){
  tempsdf <- filter(df, df$V1 == unVec[x])
  XX <- as.character(unVec[x])
  golist[[paste0(XX)]] <- tempsdf$V2
  x <- x+1
}



#### GO ANALYSIS DAY 2 ###
## DAY 2 UPREGULATED GO TERMS

#Create GeneList object which is a named factor that indicates which
#genes are interesting and which not

#Full gene name list
#geneNames <- names(golist)

#set gene universe based on the list of genes I used for the DESeq2 analysis


geneNames2 <- row.names(tempdf_2)

myInterestingGenes2u <- fastaupd2$name
geneList2u <- factor(as.integer(geneNames2 %in% myInterestingGenes2u))
names(geneList2u) <- geneNames2
str(geneList2u)



#Build GOdata object, with BP ontology (biological process)

GOdata2u <- new("topGOdata", ontology = "BP", allGenes = geneList2u, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = golist)

#test for significance
classic_fisher_result2u = runTest(GOdata2u, algorithm = 'classic', statistic = 'fisher')
weight_fisher_result2u = runTest(GOdata2u, algorithm = 'weight01', statistic = 'fisher')

#generate a table of results
allGO2u = usedGO(GOdata2u)
all_res2u = GenTable(GOdata2u, weightFisher = weight_fisher_result2u, classicFisher = classic_fisher_result2u, orderBy = 'weightFisher', topNodes = length(allGO2u))

sig_all_res2u <- filter(all_res2u, weightFisher < 0.05)


# Store significantly enriched processes and calculate percent enrichment
all_res_sig2u <- filter(all_res2u, weightFisher < 0.5)
PercentEnrichment2u <- all_res_sig2u$Significant/all_res_sig2u$Annotated
all_res_sig2u$PercentEnrichment <- PercentEnrichment2u

### DOWNREGULATED DAY 2 ###

myInterestingGenesD2 <- fastadownd2$name
geneListD2 <- factor(as.integer(geneNames2 %in% myInterestingGenesD2))
names(geneListD2) <- geneNames2
str(geneListD2)

#Build GOdata object with BP ontology (biological process)

GOdataD2 <- new("topGOdata", ontology = "BP", allGenes = geneListD2, nodeSize = 2, annot = annFUN.gene2GO, gene2GO = golist)

#test for significance
classic_fisher_resultD2 = runTest(GOdataD2, algorithm = 'classic', statistic = 'fisher')
weight_fisher_resultD2 = runTest(GOdataD2, algorithm = 'weight01', statistic = 'fisher')

#generate a table of results
allGOD2 = usedGO(GOdataD2)
all_resD2 = GenTable(GOdataD2, weightFisher = weight_fisher_resultD2, classicFisher = classic_fisher_resultD2, orderBy = 'weightFisher', topNodes = length(allGOD2))


sig_all_resD2 <- filter(all_resD2, weightFisher < 0.05)

# Store significantly enriched processes and calculate percent enrichment
PercentEnrichmentD2 <- sig_all_resD2$Significant/sig_all_resD2$Annotated
sig_all_resD2$PercentEnrichment <- PercentEnrichmentD2

#### GO ANALYSIS DAY 4 ###

## UPREGULATED DAY 4 ##

#Create GeneList object which is a named factor that indicates which
#genes are interesting and which not

#Full gene name list
#geneNames <- names(golist)

#set gene universe based on the list of genes I used for the DESeq2 analysis
geneNames4 <- row.names(tempdf_4)

myInterestingGenes4u <- fastaupd4$name
geneList4u <- factor(as.integer(geneNames4 %in% myInterestingGenes4u))
names(geneList4u) <- geneNames4
str(geneList4u)

#Build GOdata object with BP ontology (biological process)

GOdata4u <- new("topGOdata", ontology = "BP", allGenes = geneList4u, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = golist)

#test for significance
classic_fisher_result4u = runTest(GOdata4u, algorithm = 'classic', statistic = 'fisher')
weight_fisher_result4u = runTest(GOdata4u, algorithm = 'weight01', statistic = 'fisher')

#generate a table of results
allGO4u = usedGO(GOdata4u)
all_res4u = GenTable(GOdata4u, weightFisher = weight_fisher_result4u, classicFisher = classic_fisher_result4u, orderBy = 'weightFisher', topNodes = length(allGO4u))


# Store significantly enriched processes and calculate percent enrichment
all_res_sig4u <- filter(all_res4u, weightFisher < 0.05)
PercentEnrichment4u <- all_res_sig4u$Significant/all_res_sig4u$Annotated
all_res_sig4u$PercentEnrichment <- PercentEnrichment4u


# What genes are in Perception of Sweet Taste GO term
invgolist <- inverseList(golist)

sweet_taste_genes <- invgolist["GO:0050916"]
sweet_taste_genes$`GO:0050916`%in%fastaupd4$name
sweet_taste_genes $`GO:0050916`


#g12173.t1 and g12173.t3 are significantly upregulated

#g12173 is malvolio-gene (BLASTED)



### FOR DOWNREGULATED GENES ###

myInterestingGenesD4 <- fastadownd4$name
geneListD4 <- factor(as.integer(geneNames4 %in% myInterestingGenesD4))
names(geneListD4) <- geneNames4
str(geneListD4)

#Build GOdata object with BP ontology (biological process)

GOdataD4 <- new("topGOdata", ontology = "BP", allGenes = geneListD4, nodeSize = 2, annot = annFUN.gene2GO, gene2GO = golist)

#test for significance
classic_fisher_resultD4 = runTest(GOdataD4, algorithm = 'classic', statistic = 'fisher')
weight_fisher_resultD4 = runTest(GOdataD4, algorithm = 'weight01', statistic = 'fisher')

#generate a table of results
allGOD4 = usedGO(GOdataD4)
all_resD4 = GenTable(GOdataD4, weightFisher = weight_fisher_resultD4, classicFisher = classic_fisher_resultD4, orderBy = 'weightFisher', topNodes = length(allGOD4))


# Store significantly enriched processes and calculate percent enrichment
all_res_sig4d <- filter(all_resD4, weightFisher < 0.05)
PercentEnrichment4d <- all_res_sig4d$Significant/all_res_sig4d$Annotated
all_res_sig4d$PercentEnrichment <- PercentEnrichment4d

#### GO Visualizations ###

#read in data
fr_tbl_u4 <- read_xlsx("Dataset2.xlsx", sheet = "Freq_tbl_ancestor")

#aggregate frequency table
piedata <- aggregate(fr~Parent.Process, FUN = "sum", data = fr_tbl_u4)


pieplot <- ggplot(piedata, aes(x="", y = fr, fill=Parent.Process))+
  geom_bar(stat="identity", width = 1, color = "white")+
  coord_polar("y", start = 0)+
  scale_fill_manual(values=piepalet, name = "GO term")+
  ggtitle("A) High level summary of enriched GO terms (Biological process)")+
  theme_void()

#Extract legend
leg <- get_legend(pieplot)
leg <- as_ggplot(leg)

#remove legend from figure
pieplot <- pieplot + theme(legend.position = "none")

# Enrichment top 10 by percent enrichment
all_res_sig4u$Term <- as.factor(all_res_sig4u$Term)
ars10 <- head(all_res_sig4u, 10)
Parent <- as.character(fr_tbl_u4$Parent.Process[1:10])

ars10$Parent <- Parent



#plot
brplt <- ggplot(data = ars10, aes(x=(PercentEnrichment*100), y = Term, fill = Parent))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = barpalet)+
  theme_classic()+
  xlab("Percent enrichment")+
  ylab("GO term (Biological process")+
  ggtitle("B) Top 10 significantly enriched GO terms (Biological process)")+
  guides(fill = "none")

#combine plots and add common legend
goplot <- ggarrange(pieplot, brplt, nrow=2, common.legend = F)
fig3 <- ggarrange(goplot, leg, ncol=2)



#### MALVOLIO FORAGING ####

#Read in datafile, omit colony 7 due to escaped workers
fda <- read_xlsx("Dataset3.xlsx", sheet = "Mlv_foraging") %>% filter(Box != 7)

#re-classify
fdacol <- c("OC", "Box", "Day", "Time", "Food_type", "Treatment")
fda[fdacol] <- lapply(fda[fdacol], factor)

### Model building ###

#poisson model
mlvfrg <- glmmTMB(Foragers ~ Treatment * Food_type * Day + (1|OC) + (1|Box/Time), data = fda, family = "poisson")

#Diagnostics with DHARMa
plot(simulateResiduals(mlvfrg))
testOutliers(mlvfrg, type = "bootstrap")
testZeroInflation(mlvfrg)

#negative binomial model
mlvfrg_nb1 <- glmmTMB(Foragers ~ Treatment * Food_type * Day +  (1|OC) + (1|Box/Time), data = fda, family = "nbinom1")

#Diagnostics
plot(simulateResiduals(mlvfrg_nb1))

#model comparison
anova(mlvfrg, mlvfrg_nb1)
# negative binomial model better fit

# ANOVA of fixed effects
Anova(mlvfrg_nb1)

#pairwise comparisons
#Treatment over food
emmeans(mlvfrg_nb1, pairwise~Treatment|Food_type)

#within treatment comparison of food types per day
emmeans(mlvfrg_nb1, pairwise~Food_type|Treatment|Day)

#within days, comparison of treatments on food types
emmeans(mlvfrg_nb1, pairwise~Treatment|Day|Food_type)

## PLOTTING ##

line2avg <- ggplot(filter(fda, Food_type == 2), aes(x=as.factor(Day), y=Foragers, group=Treatment, color=Treatment))+
  stat_summary(aes(group=Treatment), fun = mean, geom="point", size =2.5, position = position_dodge(0.3))+
  stat_summary(aes(group=Treatment), fun = mean, geom="line", size=0.4, position = position_dodge(0.3))+
  stat_summary(aes(group=Treatment), fun.data = mean_cl_boot, geom = "errorbar", width=0.2, size=0.4, position = position_dodge(0.3))+
  labs(y="Average foraging frequency")+
  ggtitle(expression('B) Foraging on ROS food per day'))+
  scale_color_manual(values = palet, name = "Treatment", labels=expression("Control", italic(B.bassiana)))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  geom_text(x=4, y = 0.5464, label = "**", color = "black", size = 2.5)+
  geom_text(x=5, y = 0.5464, label = "***", color = "black", size = 2.5)+
  geom_text(x=6, y = 0.5464, label = "**", color = "black", size = 2.5)+
  geom_text(x=7, y = 0.5464, label = "***", color = "black", size = 2.5)+
  geom_text(x=8, y = 0.5464, label = "**", color = "black", size = 2.5)+
  geom_text(x=9, y = 0.5464, label = "**", color = "black", size = 2.5)+
  geom_text(x=10, y = 0.5464, label = "**", color = "black", size = 2.5)



line1avg <- ggplot(filter(fda, Food_type == 1), aes(x=as.factor(Day), y=Foragers, group=Treatment, color=Treatment))+
  stat_summary(aes(group=Treatment), fun = mean, geom="point", size =2.5, position = position_dodge(0.3))+
  stat_summary(aes(group=Treatment), fun = mean, geom="line", size=0.4, position = position_dodge(0.3))+
  stat_summary(aes(group=Treatment), fun.data = mean_cl_boot, geom = "errorbar", width=0.2, size=0.4, position = position_dodge(0.3))+
  labs(y="Average foraging frequency")+
  ggtitle(expression('A) Foraging on standard food per day'))+
  scale_color_manual(values = palet, name = "Treatment", labels=expression("Control", italic(B.bassiana)))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  geom_text(x=5, y = 2.09, label = "***", color = "black", size = 2.5)+
  geom_text(x=6, y = 2.09, label = "***", color = "black", size = 2.5)+
  geom_text(x=7, y = 2.09, label = "***", color = "black", size = 2.5)+
  geom_text(x=8, y = 2.09, label = "***", color = "black", size = 2.5)+
  geom_text(x=9, y = 2.09, label = "**", color = "black", size = 2.5)+
  geom_text(x=13, y = 2.09, label = "*", color = "black", size = 2.5)



#### QPCR ANALYSIS ####

# Read in data file
mdf2way <- read_xlsx("Dataset3.xlsx", sheet = "MCMC")

#re-classifying columns
colf <- c("Colony", "Nest", "Day", "Treatment", "sample", "Plate")
mdf2way[colf] <- lapply(mdf2way[colf], factor)
mdf2way <- droplevels(mdf2way)

#read in primer efficiency
meff <- as.data.frame(read_xlsx("Dataset3.xlsx", sheet = "meff"))
meff$efficiency <- as.numeric(meff$efficiency)

#prepare MCMC.qpcr dataframe
qs2wayc <- cq2counts(data = mdf2way, genecols = c(7:8), condcols = c(1:6), effic = meff, Cq1 = 'formula')

#run MCMC.qpcr model
mm2wayc <- mcmc.qpcr(
  fixed = "Day + Treatment + Treatment:Day",
  random = c("Colony", "Nest", "Plate"),
  data = qs2wayc,
  controls = "Act",
  m.fix = 1.2)

summary(mm2wayc)

#extract summary data and make trellisplot
sm3 <- HPDsummary(mm2wayc, data = qs2wayc, relative = F)
trlplt <- trellisByGene(sm3, xFactor = "Day", groupFactor = "Treatment")

#filter out just the MalG expression
trldat <- as.data.frame(trlplt[1]$data) %>%
  filter(gene == "MalG")

#make trellisplot
trlplot <- ggplot(data = trldat, aes(x = Day, y = mean, group = Treatment, colour = Treatment))

# customize trellisplot for figure 4
gpl=ggplot(trldat,
            aes(x = Day, group = Treatment, colour = Treatment, y = mean)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), lwd = 0.4, width = 0.2, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) + geom_point(position = position_dodge(0.3), size = 2.5) +
  theme_bw() +
  ggtitle(expression('C) Expression of' ~italic('Malvolio')))+
  ylab("log(abundance)") +
  theme(legend.position="bottom")+
  scale_color_manual(values = palet, name = "Treatment", labels=expression("Control", italic(B.bassiana)))+
  scale_x_discrete(labels=c("01" = "1", "02" = "2", "03" = "3", "04" = "4", "05" = "5",
                            "06" = "6", "07" = "7", "08" = "8", "09" = "9", "10" = "10",
                            "11" = "11", "12" = "12", "13" = "13"))+
  geom_text(x=4, y = 15.9, label = "***", color = "black", size = 2.5)+
  geom_text(x=5, y = 15.9, label = "***", color = "black", size = 2.5)+
  geom_text(x=6, y = 15.9, label = "***", color = "black", size = 2.5)+
  geom_text(x=7, y = 15.9, label = "***", color = "black", size = 2.5)+
  geom_text(x=8, y = 15.9, label = "**", color = "black", size = 2.5)+
  geom_text(x=9, y = 15.9, label = "**", color = "black", size = 2.5)+
  geom_text(x=10, y = 15.9, label = "**", color = "black", size = 2.5)+
  geom_text(x=11, y = 15.9, label = "**", color = "black", size = 2.5)+
  geom_text(x=13, y = 15.9, label = "***", color = "black", size = 2.5)




# combine plots for figure 4
fig4 <- ggarrange(line1avg,line2avg, gpl, nrow=3, common.legend = T)









#### Render figures and tables ####
#uncomment for rendering

##Main manuscript
#Figure 1 - Experiment a) foraging and survival
#fig1

#figure 2 - Gene expression results
#fig2

# Figure 3 - GO terms ancestry and top 10 significant enriched GO terms
#fig3


#Figure 4 - Foraging and Malvolio expression
#fig4


## Supplemental information
#Tables for GO term analysis

#Upregulated day 2
#all_res2u

#Downregulated day 2
#all_resD2

#upregulated day 4
#all_res4u

#downregulated day 4
#all_resD4


#within treatment comparison of food types per day
#emmeans(mlvfrg_nb1, pairwise~Food_type|Treatment|Day)

#results of MCMC.qpcr analysis
# extract and order output of mcmc model summary
#smr <- summary(mm2wayc)
#sm23 <- as.data.frame(smr$solutions)
#sm23$Combination <- row.names(sm23)
#row.names(sm23) <- NULL
#names(sm23)[2] <- "CI_lower"
#names(sm23)[3] <- "CI_higher"
#names(sm23)[5] <- "p"
#sm23 <- sm23[,c(6,1,2,3,4,5)]

# Extract the combination outputs (effect of treatment in days)
#sm23 <- tail(sm23, 26)

#Write table with significant p values highlighted
#tbl1 <- nice_table(sm23, highlight = T)
#tbl1

