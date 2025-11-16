#####################################################################################
# Gene expression from human brain across development: brainspan
####################################################################################
#BiocManager::install("edgeR")
library(edgeR)

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gplots)


#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

######################################################################
# Basic function to convert mouse to human gene names
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# convertMouseGeneList <- function(x){
#   #for hg38:
#   genesV2 = getLDS(attributes = c("mgi_symbol","ensembl_gene_id"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol","ensembl_gene_id","entrezgene"), martL = human, uniqueRows=T)
#   humanx <- unique(genesV2[,])
#   
#   # Print the first 6 genes found to the screen
#   print(head(humanx))
#   return(humanx)
# }

######################################################################
#Change to your working directory
setwd("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/")
######################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
#Initial data processing, sample cleaning and differential expression analysis 
############################################################################################################################
############################################################################################################################
############################################################################################################################

############################################################################################################################
#make RPKM data
############################################################################################################################

#loop for combining GeneID and counts for each txt file
path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/counts"

setwd(path)

# out.file<-""
# group<-""
# file.names <- dir(path, pattern ='*collapsedUMI.counts.txt')

# for(i in 1:length(file.names)){
#   file <- read.table(file.names[i],skip=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#   counts <- data.frame(file$Geneid,file[,7])
#   #extract file name without count.txt
#   m <- as.data.frame(file.names[i])
#   names(m) <-  c("split")
#   m$split <- as.character(m$split)
#   m <-tidyr::separate(m,split,into= c("name","stuff","stuff2","stuff3"),sep="\\.")
#   
#   print(m$name)
#   #replace column headings with names
#   names(counts) <- c("Geneid",m$name)
#   
#   #write to output file
#   out.file <- counts
#   write.table(out.file, file =paste(m$name, ".out.txt", sep=""), sep="\t",quote=F)
# }


#take files from loop above as new input
input.files <-dir(path, pattern ='*out.txt')
input.files <- as.data.frame(input.files)
names(input.files) <- c("filenames")

#get treatments names
input.files$Sample.ID <- gsub(".out.txt","",input.files$filenames)

#read in info
Info3 <- read.csv("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/MasterExperimentInfo.csv")

#remove PA92FC30 based on MSDS
Info3 <- Info3 %>% filter(Sample.ID !="PA92FC30")


setwd("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/counts")
path <- c("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/counts")
RG <- readDGE(Info3$filenames, group=Info3$Group)
colnames(RG$counts) <- gsub(".out","",colnames(RG$counts))

#save
path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/WGCNA"
setwd(path)
save(RG,Info3, file="RNAseqDEGlist.RData")

# number of unique transcripts
dim(RG)
#53801   177

###########################################################
# Filter for one count per million in at least 1/4 of libraries = 130/4 = 32
allsamples <- table(Info3$Dam.Treatment.Group,Info3$Pup.Treatment.Group,Info3$Sex.x, Info3$Brain.Region)
allsamples

#smallest group is n=6
keep <- rowSums(cpm(RG)>1)>=6
RGcpm <- RG[keep,]
dim(RGcpm)
# 21283   177

geneid <- rownames(RGcpm) #ensemble IDs for biomart

##################################################################################
#read in gene info
##################################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)

#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#us_mart <- useEnsembl(biomart = "ensembl", mirror = "uswest")
us_mart <- useEnsembl(biomart = "ensembl", mirror = "uswest",dataset = "mmusculus_gene_ensembl")

#head(listAttributes(human),100)
#attributes <- listAttributes(us_mart)
#attributes[grep("exon", attributes$name),]

genes <- getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","mgi_symbol","external_gene_name","description"), 
               filter= "ensembl_gene_id",
               values = rownames(RGcpm), 
               mart = us_mart)



genes$genelength <- abs(genes$end_position - genes$start_position) 

#remove duplicates:keeps only first entry
genes <- genes[!duplicated(genes$ensembl_gene_id),]

#match order: df[match(target, df$name),]
genes <- genes[match(rownames(RGcpm),genes$ensembl_gene_id),]

##########################################################################
#add gene info to DEGlist object
##########################################################################

#add into DEGlist
RGcpm$genes <- genes

#save
save(RG,RGcpm,genes, file="DEGlist.RData")

###########################################################
#filtering plot
###########################################################
library(RColorBrewer)
nsamples <- ncol(RGcpm)

colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)

#plot:
pdf('FilteringCPM_plots.pdf')
par(mfrow=c(1,2))

#prefilter:
lcpm <- cpm(RG, log=TRUE, prior.count=2)


plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")

#filtered data
#og-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(RGcpm, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()


###########################################################
###########################################################
#reset library sizes
RGcpm$samples$lib.size <- colSums(RGcpm$counts)


#plot library sizes
pdf('LibrarySizes.pdf',w=30,h=8)
barplot(RGcpm$samples$lib.size,names=colnames(RGcpm),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# Get log2 counts per million
logcounts <- cpm(RGcpm,log=TRUE)
# Check distributions of samples using boxplots
pdf('NonNormalizedLogCPM.pdf',w=30,h=10)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

##########################################################################
#Diagnosis * treatment + sex with subjects correlation removed
##################################################################################
#make design matrix
Info3$Sex <- factor(Info3$Sex.x, levels=c("M","F"))

Info3$Dam.Treatment.Group <- factor(Info3$Dam.Treatment.Group,levels = c("Saline","PolyIC"))

Info3$Pup.Treatment.Group <- factor(Info3$Pup.Treatment.Group,levels = c("Saline","Treg"))

Info3$Brain.Region <- factor(Info3$Brain.Region)

Info3$Group <- factor(Info3$Group)

#set design matrix
design <- model.matrix(~0+Group, data=Info3) #if using a 0 intercept must set up contrasts

colnames(design)


##################################################################################
#Normalization TMM and Voom
##################################################################################
#The two steps refer to different aspects of normalization. CPM "normalization" accounts for library size differences between samples, and produces normalized values that can be compared on an absolute scale (e.g., for filtering). TMM normalization accounts for composition bias, and computes normalization factors for comparing between libraries on a relative scale. CPM normalization doesn't account for composition bias, and TMM normalization doesn't produce normalized values. Thus, you need both steps in the analysis pipeline. This isn't a problem, as the two steps aren't really redundant.
#TMM normalization for library composition
DGE=calcNormFactors(RGcpm,method =c("TMM")) 

pdf('VoomTrend.pdf',w=6,h=4)
v=voom(DGE,design,plot=T)
dev.off()

corfit <- duplicateCorrelation(v, design, block=Info3$Mouse.Number)

fit <- lmFit(v, design, block = Info3$Mouse.Number, correlation = corfit$consensus)


#compute RPKM off adjusted library sizes
RPKM <- rpkm(DGE,gene.length =DGE$genes$genelength,normalized.lib.sizes = TRUE, log=F)

#save
save(RG,RGcpm,Info3,design,genes, RPKM,v,file="DEGlist.RData")


############################################################################################################################
############################################################################################################################
############################################################################################################################
#WGCNA
############################################################################################################################
############################################################################################################################
############################################################################################################################

#make WGCNA table
#remove extra columns:
DF <- RPKM

#make gene matrix with genes as columns and rows as samples
DF2 <- t(DF)
DF2 <- as.data.frame(DF2)
write.csv(DF2,"WGCNA_RPKMinputmatrix.csv")

################## process data ########################
#detect genes with missing values or samples with missing values
gsg = goodSamplesGenes(DF2, verbose = 3);
gsg$allOK # TRUE

#remove the genes with too many missing values
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(DF2)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(DF2)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  DF2 = DF2[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(DF2, verbose = 3);
gsg$allOK # TRUE


#gene with RPKM value of 2 or higher in at least one sample.
#max(col)>=2
#https://stackoverflow.com/questions/24212739/how-to-find-the-highest-value-of-a-column-in-a-data-frame-in-r/24212879
colMax <- function(data) sapply(data, max, na.rm = TRUE)

DF3 <- DF2[, colMax(DF2) >= 2]
dim(DF2)
dim(DF3)

#log2 transform:
Log2DF3 <- log2(DF3+1)

# median absolute deviation 
#remove if = 0 https://support.bioconductor.org/p/65124/
colMad <- function(data) sapply(data, mad, na.rm = TRUE)

DF4 = Log2DF3[, colMad(Log2DF3) != 0]
dim(DF4) #6464 genes left

write.csv(DF4,"WGCNA_log2RPKMinputmatrix.csv")



#=====================================================================================
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(DF4), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)


pdf(file = "sampleClustering_preCut.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
#remove two samples E18 tube 8 and P60 tube 50
abline(h = 160, col = "red");

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep. > all in this case
keepSamples = (clust==1)
datExpr = DF4[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#phenotype data:
#Read in experient information sheet
names(Info3)
info <- Info3

#filter to include only seq files:
info$Sample.ID <- as.character(info$Sample.ID)

infowt <- info[which(info$Sample.ID %in% rownames(datExpr)),]

table(infowt$Sample.ID %in% rownames(datExpr)) #177 samples


#match sample names and trait data:
Samples = rownames(datExpr)
traitRows = match(Samples, infowt$Sample.ID)
datTraits = infowt[traitRows, ]

datTraits<- as.data.frame(datTraits)
datTraits$Sample.ID

#change to numeric: sex #F =2 = red,  M= 1 = white
datTraits$sex2 <- as.numeric(as.factor(datTraits$Sex)) 

#change to numeric: dam treatment 1=saline, 2=polyIC
datTraits$Dam.Treatment.Group2 <- as.numeric(datTraits$Dam.Treatment.Group)

#change to numeric: pup treatment 1=saline, 2=treg
datTraits$Pup.Treatment.Group2 <- as.numeric(datTraits$Pup.Treatment.Group)

#change to numeric: 1=Cerebellum, 2= Frontal, 3= Cortex Hippocampus 
datTraits$Brain.Region2 <- as.numeric(datTraits$Brain.Region)
  
#numeric traits only
numerictraits <- datTraits[,c(20,21,19,22)]
rownames(numerictraits) <- datTraits$Sample.ID

#We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors1 = numbers2colors(numerictraits$Dam.Treatment.Group2,signed = FALSE);
  traitColors2 = numbers2colors(numerictraits$Pup.Treatment.Group2,signed = FALSE);
  traitColors3 = numbers2colors(numerictraits$sex2,signed = FALSE);
  traitColors <-cbind(traitColors1,traitColors2,traitColors3)
  colnames(traitColors) <- c("Dam treatment","Pup treatment","sex")
# Plot the sample dendrogram and the colors underneath.

pdf(file = "sampleClustering_withtraits.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = colnames(traitColors),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(numerictraits,datTraits,datExpr, file = "WTWGCNAnetworkConstruction-inputdata.RData")


#=====================================================================================
#choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency 
# Choose a set of soft-thresholding powers
powers = c(c(1:30))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",corFnc="bicor" ,verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="black")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

write.csv(sft, "SoftThresholding_powercalc.csv")
#lowest power for which the scale-free topology fit index of > 0.8 ->5

#=====================================================================================
#Run this section on an external server, save the network data and load into R studio
#=====================================================================================
#Co-expression similarity and adjacency
library(WGCNA)
allowWGCNAThreads()

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#load data



Allnet = blockwiseModules(datExpr, maxBlockSize = 50000,
                          power = 5, TOMType = "signed", minModuleSize = 20,
                          networkType = "signed",
                          corType = "bicor", #biweight midcorrelation
                          maxPOutliers = 0.05, #forces bicor to never regard more than the specified proportion of samples as outliers.
                          reassignThreshold = 0,
                          numericLabels = TRUE,
                          saveTOMs = FALSE,
                          nThreads = 12,
                          #saveTOMFileBase = "TOM-BilboMGExpression",
                          verbose = 3)

save(Allnet, file = "NetworkConstruction-auto_WTWGCNA.RData")
#=====================================================================================
table(Allnet$colors)

# Convert labels to colors for plotting
AllnetModuleColors = labels2colors(Allnet$colors)


AllnetTree = Allnet$dendrograms[[1]]

#plot the gene dendrogram and the corresponding module colors
sizeGrWindow(8,6);

pdf(file = "AllnetDendrogram_allsamples.pdf", wi = 8, he = 6)
plotDendroAndColors(AllnetTree, AllnetModuleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()


#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = AllnetModuleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

pdf(file = "Pretrim_allnet_dendrogram.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=0.6, col = "red")

dev.off()



#=====================================================================================
#merge down to fewer modules based on plot height
MEDissThres = 0.1

# Call an automatic merging function
merge = mergeCloseModules(datExpr, AllnetModuleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#replot
sizeGrWindow(12, 9)
pdf(file = "Allnet_postmerg.1MEClustering.pdf", wi = 9, he = 6)

plotDendroAndColors(AllnetTree, mColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()

#=====================================================================================
# Rename to moduleColors
moduleColors = mColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, AllnetTree, datTraits, numerictraits,datExpr,Allnet,file = "allsamplesWTMGtimecourse-networkpostmerge.1.RData")

#load("allsamplesbilboMGtimecourse-networkpostmerge.25.RData")
#=====================================================================================
#repeat clustering of new MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
pdf(file = "Pretrim_allnet_dendrogram_postmerge.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes after merging modules",
     xlab = "", sub = "")

dev.off()

#=====================================================================================
#=====================================================================================
#correlate eigengenes with external traits and look for the most signicant associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, numerictraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#FDR correct the pvalues
#sapply(pval,p.adjust,method="fdr") #per column

#FDR correction on entire matrix
FDR <- matrix(p.adjust(as.vector(moduleTraitPvalue), method='fdr'),ncol=ncol(moduleTraitPvalue))
colnames(FDR) <- colnames(moduleTraitPvalue)
rownames(FDR) <- rownames(moduleTraitPvalue)


corout <- cbind(moduleTraitCor,FDR)
write.csv(corout, "PearsonsCorrelations_mod-trait.csv",row.names = T)

#We color code each association by the correlation value:

sizeGrWindow(10,6)

pdf(file = "module-traitrelationshipsFDRcorrected.pdf", wi = 8.5, he = 11)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\ (", #pearson correlation coefficient, space, (FDR corrected pvalue)
                    signif(FDR, 4), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(numerictraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

#=====================================================================================
#ANOVA for trait associations
#=====================================================================================

library(nlme)
library(lsmeans)
datTraits <- as.data.frame(datTraits)
datTraits$sex <- factor(datTraits$sex)
datTraits$Dam.Treatment.Group <- factor(datTraits$Dam.Treatment.Group)
datTraits$Pup.Treatment.Group <- factor(datTraits$Pup.Treatment.Group)
datTraits$Brain.Region <- factor(datTraits$Brain.Region)
datTraits$subject <- datTraits$Sample.ID


MEs$Sample.ID <- rownames(MEs)
datComb <- merge(datTraits,MEs, by = "Sample.ID")
datComb2 <- datComb %>% gather(module, ME, 25:ncol(datComb))

write.csv(datComb2,"Moduel-ME_dataforANOVA.csv")

mod <- unique(datComb2$module)
#remove grey
mod <- head(mod, -1)

anova_out <- NULL
for (i in mod) {
  print(i)
  #select data
  dattmp <- datComb2 %>% dplyr::filter(module == i)
  #model
  tmp <- lme(ME ~ Dam.Treatment.Group*Pup.Treatment.Group*sex*Brain.Region, ~1|Sample.ID, data = dattmp)
  #model indcludes random effects subject ID
  anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
  #write anova to output file
  anova <- as.data.frame(anova)
  anova$condition <- rownames(anova)
  anova <- anova[-1,]
  anova2 <- anova %>% gather(stats,values,1:4) %>%
    mutate(col = paste(condition,stats, sep="_")) %>%
    dplyr::select(col,values) %>%
    spread(col,values)
  anova2$module <- paste(i)
  #save to loop
  anova_out <- rbind(anova_out,anova2)
}


write.csv(anova_out,"MixedModelAnova_moudle-traits.csv")



#include post hocs
anova_out <- NULL
ME_posthoc <- NULL
for (i in mod) {
  print(i)
  #select data
  dattmp <- datComb2 %>% dplyr::filter(module == i)
  #model
  tmp <- lme(ME ~ Dam.Treatment.Group*Pup.Treatment.Group*sex*Brain.Region, ~1|Sample.ID, data = dattmp)
  #model indcludes random effects subject ID
  anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
  #write anova to output file
  anova <- as.data.frame(anova)
  anova$condition <- rownames(anova)
  anova <- anova[-1,]
  anova2 <- anova %>% gather(stats,values,1:4) %>%
    mutate(col = paste(condition,stats, sep="_")) %>%
    dplyr::select(col,values) %>%
    spread(col,values)
  anova2$module <- paste(i)
  #save to loop
  anova_out <- rbind(anova_out,anova2)
  
  ## POSTHOC starts here ##
  #create a reference grid model object
  ME_refgrid <- ref.grid(tmp)
  
  #creates a fitted model using the reference grid, based on treatment and brain region interactions
  ME_lsmeans <- lsmeans(ME_refgrid, ~Dam.Treatment.Group*Pup.Treatment.Group|sex*Brain.Region)
  
  #summarize paired comparisons
  ME_lsmeans_summary <- summary(pairs(ME_lsmeans, adjust = "none"))
  ME_lsmeans_summary <- as.data.frame(ME_lsmeans_summary)
  
  #posthocHSD = contrast(refgrid, method = "pairwise",adjust = "none")
  #outsum <- as.data.frame(summary(posthocHSD))
  
  
  #get only control vs LPS
  # outsum <- outsum[grepl("CONTROL", outsum$contrast),]
  
  #correct for multiple comparisons
  ME_lsmeans_summary$padjust <- p.adjust(ME_lsmeans_summary$p.value, method = "BH")
  
  #paste module information
  ME_lsmeans_summary$module <- paste(i)
  
  #save to output
  ME_posthoc <- rbind(ME_posthoc,ME_lsmeans_summary)

}

write.csv(ME_posthoc,"BHposthocs_moudle-traits.csv")


#graph MEs
pdf(file = "boxplot_MEs.pdf", wi = 22, he = 30)
# 
p <- ggplot(datComb2, aes(x=Dam.Treatment.Group, y=ME), group = Pup.Treatment.Group) +
  facet_wrap(~module*Brain.Region*Sex.x,ncol = 6, scales = "free") +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Pup.Treatment.Group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=Pup.Treatment.Group, colour=Sex.x)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        #  legend.position="none",
        strip.text.x = element_text(size = 20))+
  xlab(label = c("Dam Treatment")) +
  ylab(label = c("ME"))
p
dev.off()
# 


#graph MEs Cerebellum ME blue
pdf(file = "boxplot_MEs_CB_blue.pdf", wi = 6, he = 4)
# 
CBblue <- datComb2 %>% filter(module == "MEblue") %>% filter(Brain.Region == "Cerebellum")

p <- ggplot(CBblue, aes(x=Dam.Treatment.Group, y=ME), group = Pup.Treatment.Group) +
  facet_wrap(~Sex.x,ncol = 2) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Pup.Treatment.Group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=Pup.Treatment.Group)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        #  legend.position="none",
        strip.text.x = element_text(size = 20))+
  xlab(label = c("Dam Treatment")) +
  ylab(label = c("ME"))
p
dev.off()
# 
#graph MEs Cerebellum ME turquoise
pdf(file = "boxplot_MEs_CB_turquoise.pdf", wi = 6, he = 4)
# 
CBturquoise <- datComb2 %>% filter(module == "MEturquoise") %>% filter(Brain.Region == "Cerebellum")

p <- ggplot(CBturquoise, aes(x=Dam.Treatment.Group, y=ME), group = Pup.Treatment.Group) +
  facet_wrap(~Sex.x,ncol = 2) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Pup.Treatment.Group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=Pup.Treatment.Group)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        #  legend.position="none",
        strip.text.x = element_text(size = 20))+
  xlab(label = c("Dam Treatment")) +
  ylab(label = c("ME"))
p
dev.off()
# 

#graph MEs HC ME turquoise
pdf(file = "boxplot_MEs_HC_turquoise.pdf", wi = 6, he = 4)
# 
HCturquoise <- datComb2 %>% filter(module == "MEturquoise") %>% filter(Brain.Region == "Hippocampus")

p <- ggplot(HCturquoise, aes(x=Dam.Treatment.Group, y=ME), group = Pup.Treatment.Group) +
  facet_wrap(~Sex.x,ncol = 2) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Pup.Treatment.Group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=Pup.Treatment.Group)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        #  legend.position="none",
        strip.text.x = element_text(size = 20))+
  xlab(label = c("Dam Treatment")) +
  ylab(label = c("ME"))
p
dev.off()
# 

#graph MEs Cerebellum ME turquoise
pdf(file = "boxplot_MEs_FC_green.pdf",wi = 6, he = 4)
# 
FCgreen <- datComb2 %>% filter(module == "MEgreen") %>% filter(Brain.Region == "Frontal Cortex")

p <- ggplot(FCgreen, aes(x=Dam.Treatment.Group, y=ME), group = Pup.Treatment.Group) +
  facet_wrap(~Sex.x,ncol = 2) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Pup.Treatment.Group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=Pup.Treatment.Group)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        #  legend.position="none",
        strip.text.x = element_text(size = 20))+
  xlab(label = c("Dam Treatment")) +
  ylab(label = c("ME"))
p
dev.off()
# 

pdf(file = "boxplot_MEs_FC_yellow.pdf", wi = 6, he = 4)
# 
FCyellow <- datComb2 %>% filter(module == "MEyellow") %>% filter(Brain.Region == "Frontal Cortex")

p <- ggplot(FCyellow, aes(x=Dam.Treatment.Group, y=ME), group = Pup.Treatment.Group) +
  facet_wrap(~Sex.x,ncol = 2) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Pup.Treatment.Group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=Pup.Treatment.Group)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        #  legend.position="none",
        strip.text.x = element_text(size = 20))+
  xlab(label = c("Dam Treatment")) +
  ylab(label = c("ME"))
p
dev.off()
# 




datComb2 <- read.csv("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/WGCNA/Moduel-ME_dataforANOVA.csv")

#graph MEs
datComb2$Pup.Treatment.Group <- gsub("Treg","Tregs",datComb2$Pup.Treatment.Group)
datComb2$Treatment <- paste(datComb2$Dam.Treatment.Group, datComb2$Pup.Treatment.Group, sep="-")
datComb2$Treatment <- factor(datComb2$Treatment,levels=c("Saline-Saline","Saline-Tregs","PolyIC-Saline","PolyIC-Tregs"))

#Saline–Saline (dark blue) — #777787

#Saline–Treg (light blue) — #747788

#Poly I:C–Saline (red) — #b15b53

#Poly I:C–Treg (light red/pink) — #b2817f
cbPalette <- c("#5756f9","#92bffa","#fe0100","#fe8080") 

#graph MEs Cerebellum ME blue
pdf(file = "boxplot_MEs_CB_blue_wide.pdf", wi = 8, he = 6)
# 
CBblue <- datComb2 %>% filter(module == "MEblue") %>% filter(Brain.Region == "Cerebellum")

p <- ggplot(CBblue, aes(x = Treatment, y = ME, fill = Treatment)) +
  facet_wrap(~ Sex.x, ncol = 2) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      qs <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      setNames(qs, c("ymin", "lower", "middle", "upper", "ymax"))
    },
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  geom_point(position = position_dodge(width = 0.9), aes(group = Treatment), alpha = 0.8) +
  scale_fill_manual(values = cbPalette) +
  # two-line x labels: keep "-" on the top line
  scale_x_discrete(labels = function(x) {
    parts <- strsplit(x, "\\s*[\\-–—]\\s*", perl = TRUE)
    sapply(seq_along(parts), function(i) {
      if (length(parts[[i]]) >= 2) paste0(parts[[i]][1], " -\n", parts[[i]][2]) else x[i]
    })
  }) +
  theme_cowplot(font_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.background = element_rect(fill = "#B0C4DE"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  xlab("Treatment") +
  ylab("ME")

p

dev.off()
# 
#graph MEs Cerebellum ME turquoise
pdf(file = "boxplot_MEs_CB_turquoise_wide.pdf", wi = 8, he = 6)
# 
CBturquoise <- datComb2 %>% filter(module == "MEturquoise") %>% filter(Brain.Region == "Cerebellum")

p <- ggplot(CBturquoise, aes(x = Treatment, y = ME, fill = Treatment)) +
  facet_wrap(~ Sex.x, ncol = 2) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      qs <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      setNames(qs, c("ymin", "lower", "middle", "upper", "ymax"))
    },
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  geom_point(position = position_dodge(width = 0.9), aes(group = Treatment), alpha = 0.8) +
  scale_fill_manual(values = cbPalette) +
  # two-line x labels: keep "-" on the top line
  scale_x_discrete(labels = function(x) {
    parts <- strsplit(x, "\\s*[\\-–—]\\s*", perl = TRUE)
    sapply(seq_along(parts), function(i) {
      if (length(parts[[i]]) >= 2) paste0(parts[[i]][1], " -\n", parts[[i]][2]) else x[i]
    })
  }) +
  theme_cowplot(font_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.background = element_rect(fill = "#B0C4DE"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  xlab("Treatment") +
  ylab("ME")

p

dev.off()
# 

#graph MEs HC ME turquoise
pdf(file = "boxplot_MEs_HC_turquoise_wide.pdf", wi = 8, he = 6)
# 
HCturquoise <- datComb2 %>% filter(module == "MEturquoise") %>% filter(Brain.Region == "Hippocampus")

p <- ggplot(HCturquoise, aes(x = Treatment, y = ME, fill = Treatment)) +
  facet_wrap(~ Sex.x, ncol = 2) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      qs <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      setNames(qs, c("ymin", "lower", "middle", "upper", "ymax"))
    },
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  geom_point(position = position_dodge(width = 0.9), aes(group = Treatment), alpha = 0.8) +
  scale_fill_manual(values = cbPalette) +
  # two-line x labels: keep "-" on the top line
  scale_x_discrete(labels = function(x) {
    parts <- strsplit(x, "\\s*[\\-–—]\\s*", perl = TRUE)
    sapply(seq_along(parts), function(i) {
      if (length(parts[[i]]) >= 2) paste0(parts[[i]][1], " -\n", parts[[i]][2]) else x[i]
    })
  }) +
  theme_cowplot(font_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.background = element_rect(fill = "#B0C4DE"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  xlab("Treatment") +
  ylab("ME")

p
dev.off()
# 

#graph MEs Cerebellum ME turquoise
pdf(file = "boxplot_MEs_FC_green_wide.pdf",wi = 8, he = 6)
# 
FCgreen <- datComb2 %>% filter(module == "MEgreen") %>% filter(Brain.Region == "Frontal Cortex")

p <- ggplot(FCgreen, aes(x = Treatment, y = ME, fill = Treatment)) +
  facet_wrap(~ Sex.x, ncol = 2) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      qs <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      setNames(qs, c("ymin", "lower", "middle", "upper", "ymax"))
    },
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  geom_point(position = position_dodge(width = 0.9), aes(group = Treatment), alpha = 0.8) +
  scale_fill_manual(values = cbPalette) +
  # two-line x labels: keep "-" on the top line
  scale_x_discrete(labels = function(x) {
    parts <- strsplit(x, "\\s*[\\-–—]\\s*", perl = TRUE)
    sapply(seq_along(parts), function(i) {
      if (length(parts[[i]]) >= 2) paste0(parts[[i]][1], " -\n", parts[[i]][2]) else x[i]
    })
  }) +
  theme_cowplot(font_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.background = element_rect(fill = "#B0C4DE"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  xlab("Treatment") +
  ylab("ME")
p
dev.off()
# 

pdf(file = "boxplot_MEs_FC_yellow_wide.pdf", wi = 8, he = 6)
# 
FCyellow <- datComb2 %>% filter(module == "MEyellow") %>% filter(Brain.Region == "Frontal Cortex")


p <- ggplot(FCyellow, aes(x = Treatment, y = ME, fill = Treatment)) +
  facet_wrap(~ Sex.x, ncol = 2) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      qs <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      setNames(qs, c("ymin", "lower", "middle", "upper", "ymax"))
    },
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  geom_point(position = position_dodge(width = 0.9), aes(group = Treatment), alpha = 0.8) +
  scale_fill_manual(values = cbPalette) +
  # two-line x labels: keep "-" on the top line
  scale_x_discrete(labels = function(x) {
    parts <- strsplit(x, "\\s*[\\-–—]\\s*", perl = TRUE)
    sapply(seq_along(parts), function(i) {
      if (length(parts[[i]]) >= 2) paste0(parts[[i]][1], " -\n", parts[[i]][2]) else x[i]
    })
  }) +
  theme_cowplot(font_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.background = element_rect(fill = "#B0C4DE"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  xlab("Treatment") +
  ylab("ME")

p
dev.off()
# 


pdf(file = "boxplot_MEs_FC_brown_wide.pdf", wi = 8, he = 6)
# 
FCbrown <- datComb2 %>% filter(module == "MEbrown") %>% filter(Brain.Region == "Frontal Cortex")


p <- ggplot(FCbrown, aes(x = Treatment, y = ME, fill = Treatment)) +
  facet_wrap(~ Sex.x, ncol = 2) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      qs <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      setNames(qs, c("ymin", "lower", "middle", "upper", "ymax"))
    },
    position = position_dodge2(width = 0.9, preserve = "single")
  ) +
  geom_point(position = position_dodge(width = 0.9), aes(group = Treatment), alpha = 0.8) +
  scale_fill_manual(values = cbPalette) +
  # two-line x labels: keep "-" on the top line
  scale_x_discrete(labels = function(x) {
    parts <- strsplit(x, "\\s*[\\-–—]\\s*", perl = TRUE)
    sapply(seq_along(parts), function(i) {
      if (length(parts[[i]]) >= 2) paste0(parts[[i]][1], " -\n", parts[[i]][2]) else x[i]
    })
  }) +
  theme_cowplot(font_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.background = element_rect(fill = "#B0C4DE"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  xlab("Treatment") +
  ylab("ME")

p
dev.off()
# 






#save.image(file = "AllFiles.RData")
#load(file = "/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/WGCNA/AllFiles.RData")

#quantify associations of individual genes with our trait of interest (years) by dening Gene Signicance GS as
#(the absolute value of) the correlation between the gene and the trait

#For each module, we also dene a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression prole. This
#allows us to quantify the similarity of all genes on the array to every module.


##module membership MM 
#correlation between expression of each gene and MEs (uncorrected p values)
MEs <- MEs %>% dplyr::select(-Sample.ID)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#or
# calculate the module membership values (aka. module eigengene based
# connectivity kME):
#datKME = signedKME(PFCdatExpr, MEs)

#correlate each gene expression with variable of interest
# Define variable age containing the years column of datTrait
#dam
Dam.Treatment = as.data.frame(numerictraits$Dam.Treatment.Group2);
names(Dam.Treatment) = "Dam.Treatment"

Dam_geneTraitSignificance = as.data.frame(cor(datExpr, Dam.Treatment, use = "p"));
Dam_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Dam_geneTraitSignificance), nSamples));

names(Dam_geneTraitSignificance) = paste("GS.", names(Dam.Treatment), sep="");
names(Dam_GSPvalue) = paste("p.GS.", names(Dam.Treatment), sep="")

#pups
Pup.Treatment = as.data.frame(numerictraits$Pup.Treatment.Group2);
names(Pup.Treatment) = "Pup.Treatment"

Pup_geneTraitSignificance = as.data.frame(cor(datExpr, Pup.Treatment, use = "p"));
Pup_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Pup_geneTraitSignificance), nSamples));

names(Pup_geneTraitSignificance) = paste("GS.", names(Pup.Treatment), sep="");
names(Pup_GSPvalue) = paste("p.GS.", names(Pup.Treatment), sep="")

#Regions
Region = as.data.frame(numerictraits$Brain.Region2);
names(Region) = "Brain Region"

Region_geneTraitSignificance = as.data.frame(cor(datExpr, Region, use = "p"));
Region_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Region_geneTraitSignificance), nSamples));

names(Region_geneTraitSignificance) = paste("GS.", names(Region), sep="");
names(Region_GSPvalue) = paste("p.GS.", names(Region), sep="")

#=====================================================================================
#Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high signicance for weight as well as high module
#membership in interesting modules.


#As an example, we look at the brown module that has the highest association with age.

pdf(file = "Damtreatment_geneModuleMembershipvsGeneSignificanceAGE.pdf", wi = 10, he = 10)

par(mfrow = c(2,2));

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Dam_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Dam Treatment Gene significance for age",
                   main = paste("Module membership vs. Dam treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Dam_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Dam Treatment Gene significance for age",
                   main = paste("Module membership vs. Dam treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Dam_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Dam Treatment Gene significance for age",
                   main = paste("Module membership vs. Dam treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Dam_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Dam Treatment Gene significance for age",
                   main = paste("Module membership vs. Dam treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Dam_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Dam Treatment Gene significance for age",
                   main = paste("Module membership vs. Dam treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "grey"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Dam_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Dam Treatment Gene significance for age",
                   main = paste("Module membership vs. Dam treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()



pdf(file = "Puptreatment_geneModuleMembershipvsGeneSignificanceAGE.pdf", wi = 10, he = 10)

par(mfrow = c(2,2));

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Pup_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pup Treatment Gene significance for age",
                   main = paste("Module membership vs. Pup treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Pup_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pup Treatment Gene significance for age",
                   main = paste("Module membership vs. Pup treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Pup_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pup Treatment Gene significance for age",
                   main = paste("Module membership vs. Pup treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Pup_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pup Treatment Gene significance for age",
                   main = paste("Module membership vs. Pup treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Pup_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pup Treatment Gene significance for age",
                   main = paste("Module membership vs. Pup treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "grey"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Pup_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pup Treatment Gene significance for age",
                   main = paste("Module membership vs. Pup treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()


pdf(file = "Region_geneModuleMembershipvsGeneSignificance.pdf", wi = 10, he = 10)

par(mfrow = c(2,2));

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Region_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Region Treatment Gene significance for Region",
                   main = paste("Module membership vs. Region treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Region_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Region Treatment Gene significance for Region",
                   main = paste("Module membership vs. Region treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Region_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Region Treatment Gene significance for Region",
                   main = paste("Module membership vs. Region treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Region_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Region Treatment Gene significance for Region",
                   main = paste("Module membership vs. Region treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Region_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Region Treatment Gene significance for Region",
                   main = paste("Module membership vs. Region treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "grey"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(Region_geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Region Treatment Gene significance for Region",
                   main = paste("Module membership vs. Region treatment gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()


#Clearly, GS and MM are highly correlated, illustrating that genes highly significantly associated with a trait are 
#often also the most important (central) elements of modules associated with the trait.

#=====================================================================================

#We have found modules with high association with our trait of interest, and have identified their central players by the Module Membership measure.
#We now merge this statistical information with gene annotation and write out a file that summarizes the most important results and can be inspected in standard spreadsheet software 
# Create the starting data frame
Ensemble = colnames(datExpr)

#df <- nc %>% dplyr::select(ensembl_gene_id,mgi_symbol,description,genelength) %>% distinct()

#get matching data for each ensemble id from the row data from Brainspain input
ids <- genes[which(genes$ensembl_gene_id %in% Ensemble),]

#match order

idRows = match(Ensemble, ids$ensembl_gene_id)
ids2 = ids[idRows, ]


geneInfo0 = data.frame(ensembl_gene_id = colnames(datExpr),
                       geneSymbol = ids2$mgi_symbol,
                       description = ids2$description,
                       genelength = ids2$genelength,
                       moduleColor = moduleColors,
                       Dam_geneTraitSignificance,
                       Dam_GSPvalue,
                       Pup_geneTraitSignificance,
                       Pup_GSPvalue,
                       Region_geneTraitSignificance,
                       Region_GSPvalue,
                       t(datExpr))
# Order modules by their significance for age
modOrder = order(-abs(cor(MEs, Dam.Treatment, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}



write.csv(geneInfo0, file = "geneInfo_Master_WGCNA.csv")

#geneInfo0 <- read.csv("geneInfo_BilboMGtimecourse.csv")

#=====================================================================================
#output a character vector of genes, where the genes are the hub gene picked for each module, 
#and the names correspond to the module in which each gene is a hub.

hubs <- chooseTopHubInEachModule(datExpr,moduleColors, type = "signed",power =10) #https://support.bioconductor.org/p/46342/
#power of 2 for unsigned, 4 for signed
write.csv(hubs,"TopHubInEachModule.csv")

#These genes represent the top 10 genes per module based on kME  
topGenesKME = NULL
for (i in 1:length(colnames(geneModuleMembership))){
  KMErank = geneModuleMembership[(order(-geneModuleMembership[,i])),] #order by column, - decreasing
  KMErank$Ensemble <- rownames(KMErank)
  GenesKME = KMErank[c(1:10),c(i,7)] #where column 7 is the ensemble id
  
  #get gene info
  topKMEinfo <- ids2[which(ids2$ensembl_gene_id %in% rownames(GenesKME)),]
  
  #get kME
  merge <- merge(topKMEinfo,GenesKME,by.x="ensembl_gene_id", by.y = "Ensemble")
  colnames(merge)[9] <- c("kME") #module column
  
  merge$module <- substr(colnames(geneModuleMembership[i]),3,nchar(colnames(geneModuleMembership[i])))
  
  topGenesKME = rbind(topGenesKME,merge)
}

write.csv(topGenesKME,"Top10HubInEachModule_KME.csv")


#topGenesKME <- read.csv("TopHubInEachModule.csv")

#=====================================================================================
#Extract modules
module_colors = setdiff(unique(moduleColors), "grey")
# for (color in module_colors){
#   module=datExpr[,moduleColors==module_colors]
#   write.table(module, paste("module_",color, "log2RPKM+1expression_DLPFC.txt",sep=""), sep="\t", row.names=T, col.names=T,quote=FALSE)
#   
# }

#Look at expression patterns of these genes, as they are clustered
#heatmap colors:
#install.packages("RColorBrewer")
library("RColorBrewer")
library(gplots)

#gene-module info

heatDF <- geneInfo0 %>% filter(moduleColor != "grey") %>% arrange(moduleColor)
heatDF$moduleColor <- factor(heatDF$moduleColor)

data <- as.matrix(heatDF[,12:188]) #RPKM values
rownames(data) <- heatDF$ensembl_gene_id


#data <- t(m)
library(pheatmap)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
annotation_col <- data.frame(
  #Sample = factor(d$donor_name),
  Dam.Treatment = factor(numerictraits$Dam.Treatment.Group2, levels = unique(numerictraits$Dam.Treatment.Group2)),
  Pup.Treatment = factor(numerictraits$Pup.Treatment.Group2, levels = unique(numerictraits$Pup.Treatment.Group2)),
  Region = factor(numerictraits$Brain.Region2, levels = unique(numerictraits$Brain.Region2)),
  Sex = factor(numerictraits$sex2))

rownames(annotation_col) = colnames(data)

head(annotation_col)

annotation_row <- data.frame(heatDF$moduleColor)

rownames(annotation_row) = rownames(data)

colnames(annotation_row) <- c("module")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 2)
names(Var1) <- unique(annotation_col$Dam.Treatment)

Var2        <- sample(col_vector, 2)
names(Var2) <- unique(annotation_col$Pup.Treatment)

Var3        <-  c("cornflowerblue","darksalmon")
names(Var3) <- unique(annotation_col$Sex)

Var4        <-  sample(col_vector, 3)
names(Var4) <- unique(annotation_col$Region)

Var5 <- as.character(unique(heatDF$moduleColor))
names(Var5) <-  unique(heatDF$moduleColor)

anno_colors <- list(Dam.Treatment = Var1, Pup.Treatment = Var2, Sex = Var3, Brain.Region = Var4, moduleColor = Var5)

pdf(file = "module-Heatmap_log2RPKM_zscore.pdf", wi = 20, he = 16)
pheatmap::pheatmap(data, 
                   cluster_row = F,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   fontsize = 10,
                   fontsize_row=6, 
                   show_rownames = F,
                   fontsize_col = 10,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Gene Expression Modules")

dev.off()

#=====================================================================================
#average heatmap by condition
heatDF2 <- heatDF[,c(1,5,12:188)] #ensembl ID, modulecolour, log2PRKM values

heatDF2 <- heatDF2 %>% gather(Sample.ID, log2RPKM, 3:ncol(heatDF2))
headDF3 <- merge(heatDF2,datTraits, by = "Sample.ID", all.x=T)


heatDF4 <- headDF3 %>% group_by(moduleColor,ensembl_gene_id,Group) %>% summarize(meanlog2RPKM= mean(log2RPKM)) %>%
  spread(Group,meanlog2RPKM)%>% arrange(moduleColor)

matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
rownames(matrix) <- heatDF4$ensembl_gene_id

#order columns
colnames(matrix)
#newcolorder <- c( "7 F" , "7 M" ,"15 F", "15 M", "35 F" ,"35 M" , "84 F", "84 M")

#matrix <- matrix[,newcolorder]
#colnames(matrix)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
Dam.Treatment <- substr(colnames(matrix),1,6)
tmp <- as.data.frame(colnames(matrix))
colnames(tmp) <- c("group")
tmp <- tmp %>% separate(group,into=c("Dam","Pup","sex","region"),sep="_")

sex <- substr(colnames(matrix),(nchar(colnames(matrix))+1)-1,nchar(colnames(matrix)))

annotation_col <- tmp

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

annotation_row <- data.frame(heatDF4$moduleColor)

rownames(annotation_row) = rownames(matrix)

colnames(annotation_row) <- c("moduleColor")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 2)
names(Var1) <- unique(annotation_col$Dam)

Var2       <- sample(col_vector, 2)
names(Var2) <- unique(annotation_col$Pup)

Var3        <-  c("cornflowerblue","darksalmon")
names(Var3) <- unique(annotation_col$sex)

Var4        <-  sample(col_vector, 3)
names(Var4) <- unique(annotation_col$region)

Var5 <- as.character(unique(heatDF4$moduleColor))
names(Var5) <-  unique(heatDF4$moduleColor)

anno_colors <- list(DamTreatment = Var1, PupTreatment = Var2, Sex = Var3, BrainRegion = Var4, moduleColor = Var5)

pdf(file = "Heatmap_log2RPKM_zscore_noReplicates.pdf", wi = 11, he = 8.5)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   show_rownames = F,
                   fontsize = 12,
                   fontsize_row=6, 
                   fontsize_col = 12,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Gene Expression Modules")

dev.off()


#=====================================================================================
#split by brain region
# #=====================================================================================

##############Cerebellum###################
CB <- headDF3 %>% filter(Brain.Region == "Cerebellum")

heatDF4 <- CB %>% group_by(moduleColor,ensembl_gene_id,Group) %>% summarize(meanlog2RPKM= mean(log2RPKM)) %>%
  spread(Group,meanlog2RPKM)%>% arrange(moduleColor)

matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
rownames(matrix) <- heatDF4$ensembl_gene_id

#order columns
colnames(matrix)
newcolorder <- c("Saline_Saline_M_Cerebellum", "Saline_Saline_F_Cerebellum",
                 "Saline_Treg_M_Cerebellum","Saline_Treg_F_Cerebellum",
                 "PolyIC_Saline_M_Cerebellum","PolyIC_Saline_F_Cerebellum",
                 "PolyIC_Treg_M_Cerebellum", "PolyIC_Treg_F_Cerebellum")
  
matrix <- matrix[,newcolorder]
colnames(matrix)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
tmp <- as.data.frame(colnames(matrix))
colnames(tmp) <- c("group")
tmp <- tmp %>% separate(group,into=c("Dam","Pup","sex","region"),sep="_")

sex <- substr(colnames(matrix),(nchar(colnames(matrix))+1)-1,nchar(colnames(matrix)))

annotation_col <- tmp

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

annotation_row <- data.frame(heatDF4$moduleColor)

rownames(annotation_row) = rownames(matrix)

colnames(annotation_row) <- c("moduleColor")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 2)
names(Var1) <- unique(annotation_col$Dam)

Var2       <- sample(col_vector, 2)
names(Var2) <- unique(annotation_col$Pup)

Var3        <-  c("cornflowerblue","darksalmon")
names(Var3) <- unique(annotation_col$sex)

Var4        <-  sample(col_vector, 3)
names(Var4) <- unique(annotation_col$region)

Var5 <- as.character(unique(heatDF4$moduleColor))
names(Var5) <-  unique(heatDF4$moduleColor)

anno_colors <- list(DamTreatment = Var1, PupTreatment = Var2, Sex = Var3, BrainRegion = Var4, moduleColor = Var5)

pdf(file = "CB_Heatmap_log2RPKM_zscore_noReplicates.pdf", wi = 11, he = 8.5)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   show_rownames = F,
                   fontsize = 12,
                   fontsize_row=6, 
                   fontsize_col = 12,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Gene Expression Modules")

dev.off()

##############Frontal Cortex###################
FC <- headDF3 %>% filter(Brain.Region == "Frontal Cortex")

heatDF4 <- FC %>% group_by(moduleColor,ensembl_gene_id,Group) %>% summarize(meanlog2RPKM= mean(log2RPKM)) %>%
  spread(Group,meanlog2RPKM)%>% arrange(moduleColor)

matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
rownames(matrix) <- heatDF4$ensembl_gene_id

#order columns
#order columns
colnames(matrix)
newcolorder <- c("Saline_Saline_M_Frontal Cortex", "Saline_Saline_F_Frontal Cortex",
                 "Saline_Treg_M_Frontal Cortex","Saline_Treg_F_Frontal Cortex",
                 "PolyIC_Saline_M_Frontal Cortex","PolyIC_Saline_F_Frontal Cortex",
                 "PolyIC_Treg_M_Frontal Cortex", "PolyIC_Treg_F_Frontal Cortex")

matrix <- matrix[,newcolorder]
colnames(matrix)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
tmp <- as.data.frame(colnames(matrix))
colnames(tmp) <- c("group")
tmp <- tmp %>% separate(group,into=c("Dam","Pup","sex","region"),sep="_")

sex <- substr(colnames(matrix),(nchar(colnames(matrix))+1)-1,nchar(colnames(matrix)))

annotation_col <- tmp

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

annotation_row <- data.frame(heatDF4$moduleColor)

rownames(annotation_row) = rownames(matrix)

colnames(annotation_row) <- c("moduleColor")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 2)
names(Var1) <- unique(annotation_col$Dam)

Var2       <- sample(col_vector, 2)
names(Var2) <- unique(annotation_col$Pup)

Var3        <-  c("cornflowerblue","darksalmon")
names(Var3) <- unique(annotation_col$sex)

Var4        <-  sample(col_vector, 3)
names(Var4) <- unique(annotation_col$region)

Var5 <- as.character(unique(heatDF4$moduleColor))
names(Var5) <-  unique(heatDF4$moduleColor)

anno_colors <- list(DamTreatment = Var1, PupTreatment = Var2, Sex = Var3, BrainRegion = Var4, moduleColor = Var5)

pdf(file = "FC_Heatmap_log2RPKM_zscore_noReplicates.pdf", wi = 11, he = 8.5)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   show_rownames = F,
                   fontsize = 12,
                   fontsize_row=6, 
                   fontsize_col = 12,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Gene Expression Modules")

dev.off()

##############FHippocampus###################
HC <- headDF3 %>% filter(Brain.Region == "Hippocampus")

heatDF4 <- HC %>% group_by(moduleColor,ensembl_gene_id,Group) %>% summarize(meanlog2RPKM= mean(log2RPKM)) %>%
  spread(Group,meanlog2RPKM)%>% arrange(moduleColor)

matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
rownames(matrix) <- heatDF4$ensembl_gene_id

#order columns
colnames(matrix)
newcolorder <- c("Saline_Saline_M_Hippocampus", "Saline_Saline_F_Hippocampus",
                 "Saline_Treg_M_Hippocampus","Saline_Treg_F_Hippocampus",
                 "PolyIC_Saline_M_Hippocampus","PolyIC_Saline_F_Hippocampus",
                 "PolyIC_Treg_M_Hippocampus", "PolyIC_Treg_F_Hippocampus")

matrix <- matrix[,newcolorder]
colnames(matrix)


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
tmp <- as.data.frame(colnames(matrix))
colnames(tmp) <- c("group")
tmp <- tmp %>% separate(group,into=c("Dam","Pup","sex","region"),sep="_")


annotation_col <- tmp

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

annotation_row <- data.frame(heatDF4$moduleColor)

rownames(annotation_row) = rownames(matrix)

colnames(annotation_row) <- c("moduleColor")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 2)
names(Var1) <- unique(annotation_col$Dam)

Var2       <- sample(col_vector, 2)
names(Var2) <- unique(annotation_col$Pup)

Var3        <-  c("cornflowerblue","darksalmon")
names(Var3) <- unique(annotation_col$sex)

Var4        <-  sample(col_vector, 3)
names(Var4) <- unique(annotation_col$region)

Var5 <- as.character(unique(heatDF4$moduleColor))
names(Var5) <-  unique(heatDF4$moduleColor)

anno_colors <- list(DamTreatment = Var1, PupTreatment = Var2, Sex = Var3, BrainRegion = Var4, moduleColor = Var5)

pdf(file = "HC_Heatmap_log2RPKM_zscore_noReplicates.pdf", wi = 11, he = 8.5)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   show_rownames = F,
                   fontsize = 12,
                   fontsize_row=6, 
                   fontsize_col = 12,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Gene Expression Modules")

dev.off()


#=====================================================================================
#line plots
# #=====================================================================================

# # #scale data same as pheatmap: subtract mean and divide by std. dev
# # #https://www.biostars.org/p/223532/
# scale_rows = function(x){
#    m = apply(x, 1, mean, na.rm = T)
#    s = apply(x, 1, sd, na.rm = T)
#    return((x - m) / s)
#  }
# # 
# scale_mat = function(mat, scale){
#    if(!(scale %in% c("none", "row", "column"))){
#      stop("scale argument shoud take values: 'none', 'row' or 'column'")
#    }
#    mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
#    return(mat)
#  }
# # 
# # #data:
#  heatDF4 <- headDF3 %>% group_by(moduleColor,ensembl_gene_id,Group) %>% summarize(meanlog2exp = mean(log2RPKM)) %>%
#    spread(Group,meanlog2exp)%>% arrange(moduleColor)
# # 
#  matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
#  rownames(matrix) <- heatDF4$ensembl_gene_id
# # 
# # 
# # #scale by row:
# merPFCmatrix <- scale_mat(matrix, scale = "row")
# merPFCmatrix <- as.data.frame(merPFCmatrix)
# # 
# merPFCmatrix$module <- heatDF4$moduleColor
# # 
# merPFCDF <- merPFCmatrix %>% gather(group, scaled_expression, 1:(ncol(merPFCmatrix)-1)) %>% 
#   group_by(module,group) 
  #%>%
#  # summarize(meanscaled_expression = mean(scaled_expression)) 
# # 
#  merPFCDF <- as.data.frame(merPFCDF)
# # 
# # #arrange 
# # newcolorder <- c( "7 F" , "7 M" ,"15 F", "15 M", "35 F" ,"35 M" , "84 F", "84 M")
# # 
# # merPFCDF$age <- factor(merPFCDF$age,levels = newcolorder)
# # 
# # #set colors:
#  my_palette2 <- unique(merPFCDF$module)
# # 
#  my_palette2 <- col2hex(my_palette2)
# # 
# library(cowplot)
#  pdf(file = "lineplot_log2RPKM_scaled_noReplicates.pdf", wi = 11, he = 8.5)
# # 
#  p <- ggplot(merPFCDF, aes(x=group, y=meanscaled_expression, group=module,colour=module)) +
#    geom_line(size=2) +
#    scale_color_manual(values = my_palette2)+
#    theme(axis.text.x  = element_text(angle=45, vjust=0.5))+
#    xlab(label = c("Treatment")) +
#    ylab(label = c("Scaled log2(RPKM+1)"))
# # 
#  p
# # 
#  dev.off()
# 
# #facet by module
 
 merPFCDF2 <- merPFCDF %>% separate(group, into=c("DamTreatment","PupTreatment","sex","region"), sep="_")
 merPFCDF2$DamTreatment <- factor( merPFCDF2$DamTreatment,levels =c("Saline","PolyIC" ))

  theme_set(theme_classic(base_size = 24)) 
 cbPalette <- c("lightblue","pink")
 
 pdf(file = "boxplot_log2RPKM+1_scaled_noReplicates_facet.pdf", wi = 18, he = 30)
# 
 p <- ggplot(merPFCDF2, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
  facet_wrap(~region*module*sex,ncol = 5) +
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
 # geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
      #  legend.position="none",
        strip.text.x = element_text(size = 30))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()
# 

 pdf(file = "boxplot_log2RPKM+1_scaled_noReplicates_facet2.pdf", wi = 22, he = 30)
 # 
 p <- ggplot(merPFCDF2, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
   facet_wrap(~module*region*sex,ncol = 6) +
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
   #geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
         strip.background = element_rect(fill="#B0C4DE"),
         #  legend.position="none",
         strip.text.x = element_text(size = 20))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()
 # 
 
 pdf(file = "boxplot_log2RPKM+1_scaled_CB_bluemodule.pdf", wi = 6, he = 4)
 # 
 CBblueRPKM <- merPFCDF2 %>% filter(module == "blue") %>% filter(region == "Cerebellum")
 p <- ggplot( CBblueRPKM, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
   facet_wrap(~sex,ncol = 2) +
   
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
   #geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
         strip.background = element_rect(fill="#B0C4DE"),
         #  legend.position="none",
         strip.text.x = element_text(size = 20))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()
 
 pdf(file = "boxplot_log2RPKM+1_scaled_CB_turquoisemodule.pdf", wi = 6, he = 4)
 # 
 CBturquoiseRPKM <- merPFCDF2 %>% filter(module == "turquoise") %>% filter(region == "Cerebellum")
 p <- ggplot( CBturquoiseRPKM, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
   facet_wrap(~sex,ncol = 2) +
   
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
   #geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
         strip.background = element_rect(fill="#B0C4DE"),
         #  legend.position="none",
         strip.text.x = element_text(size = 20))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()
 
 pdf(file = "boxplot_log2RPKM+1_scaled_HC_turquoisemodule.pdf", wi = 6, he = 4)
 # 
HCturquoiseRPKM <- merPFCDF2 %>% filter(module == "turquoise") %>% filter(region == "Hippocampus")
 p <- ggplot( HCturquoiseRPKM, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
   facet_wrap(~sex,ncol = 2) +
   
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
   #geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
         strip.background = element_rect(fill="#B0C4DE"),
         #  legend.position="none",
         strip.text.x = element_text(size = 20))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()

 pdf(file = "boxplot_log2RPKM+1_scaled_FC_Greenmodule.pdf", wi = 6, he = 4)
 # 
FCgreenRPKM <- merPFCDF2 %>% filter(module == "green") %>% filter(region == "Frontal Cortex")
 p <- ggplot(FCgreenRPKM, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
   facet_wrap(~sex,ncol = 2) +
   
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
   #geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
         strip.background = element_rect(fill="#B0C4DE"),
         #  legend.position="none",
         strip.text.x = element_text(size = 20))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()
 
 pdf(file = "boxplot_log2RPKM+1_scaled_FC_yellowmodule.pdf", wi = 6, he = 4)
 # 
 FCyellowRPKM <- merPFCDF2 %>% filter(module == "yellow") %>% filter(region == "Frontal Cortex")
 p <- ggplot(FCyellowRPKM, aes(x=DamTreatment, y=scaled_expression), group = PupTreatment) +
   facet_wrap(~sex,ncol = 2) +
   
   stat_summary(geom = "boxplot", 
                fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                position = "dodge", aes(fill=PupTreatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
   #geom_point(position = position_dodge(width = 0.90),aes(group=PupTreatment, colour=sex)) + 
   scale_fill_manual(values = cbPalette) +
   theme_cowplot(font_size = 15)+
   theme(axis.text.x  = element_text(angle=45, vjust=0.5),
         strip.background = element_rect(fill="#B0C4DE"),
         #  legend.position="none",
         strip.text.x = element_text(size = 20))+
   xlab(label = c("Dam Treatment")) +
   ylab(label = c("Scaled log2(RPKM)"))
 p
 dev.off()
 
 
#=====================================================================================
#GO term enrichment:https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/
#source("http://bioconductor.org/biocLite.R");
#biocLite(c("impute", "AnnotationDBI", "GO.db", "org.Hs.eg.db", "org.Mm.eg.db"))

#install.packages("anRichment_1.01-2.tar.gz")
#source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
# installAnRichment();

library(anRichment)
options(stringsAsFactors = FALSE)


###################################################################
GOcollection = buildGOcollection(organism = "mouse")

#get entrez IDs:
entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
               filter= "ensembl_gene_id",
               values = geneInfo0$ensembl_gene_id,
               mart = mouse)

geneInfo <- merge(geneInfo0, entrez, by="ensembl_gene_id")

# evaluates the enrichment of the gene modules in the collection of GO terms

GOenrichment = enrichmentAnalysis(
   classLabels = geneInfo$moduleColor, identifiers = geneInfo$entrezgene_id,
   refCollection = GOcollection,
   useBackground = "intersection", #intersection of genes in identifiers and the organism database
  threshold = 0.05,
  thresholdType = "FDR",
   getOverlapEntrez = FALSE,
   getOverlapSymbols = TRUE,
   entrySeparator = ",",
  maxReportedOverlapGenes = 1000,
   ignoreLabels = "grey") #ignore genes not assigned to modules (grey)
# 
 collectGarbage()
# 
 names(GOenrichment)
# 
 #enrichment results are summarized in the component enrichmentTable:
 names(GOenrichment$enrichmentTable)
# 

 write.csv(GOenrichment$enrichmentTable, file = "GOenrichment-enrichmentTable.csv")
# 
# 
