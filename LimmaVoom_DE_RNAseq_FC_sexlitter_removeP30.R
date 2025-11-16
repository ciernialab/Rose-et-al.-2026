## code for taking counts.txt file outputs from FeatureCounts and importing into R for EdgeR and/or limmavoom analysis
## loops through each txt file and takes the GeneID and accepted_hits.bam counts columns (counts per gene)
## writes a new text file samplename.out.txt that contains only these two columns
## these files are then read into EdgeR 
## a targets dataframe is made that contains a list of the .out.txt files and the group they belong to
## the targets dataframe is then fed into EdgeR using readDGE function of EdgeR
##Oct 2022 AVC

#remove outlier from PCA: PA92FC30.out.txt

library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(tidyr)
library(readxl)
library(edgeR)
library(limma)
library(dplyr)

#loop for combining GeneID and counts for each txt file
path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/counts"

setwd(path)

# out.file<-""
# group<-""
# file.names <- dir(path, pattern ='*collapsedUMI.counts.txt')
# 
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

##The following code is for analysis of RNASeq in Limma Voom package
##########################################################################################

#take files from loop above as new input
input.files <-dir(path, pattern ='*out.txt')
input.files <- as.data.frame(input.files)
names(input.files) <- c("filenames")

#get treatments names
input.files$Sample.ID <- gsub(".out.txt","",input.files$filenames)

#remove outlier from PCA:
input.files <- input.files %>% filter(Sample.ID != "PA92FC30")

##########################################################################################
#read in subject info
##########################################################################################

Info3 <- read.csv("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/MasterExperimentInfo.csv")
##########################################################################################
#filter for only FC samples
FCsamples <- Info3 %>% filter(Brain.Region == "Frontal Cortex")

FCsamples$Group <- gsub("_Frontal Cortex","",FCsamples$Group)

#remove PA92FC30"
FCsamples <- FCsamples %>% filter(Sample.ID != "PA92FC30")

table(FCsamples$Group, FCsamples$Dam..)
##########################################################################################
# Make DEGList


setwd("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/counts")
path <- c("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/counts")
RG <- readDGE(FCsamples$filenames, group=FCsamples$Group)
colnames(RG$counts) <- gsub(".out","",colnames(RG$counts))

#save
path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/FrontalCortex/sexlitter_removePA92FC30"
setwd(path)
save(RG,Info3, file="RNAseqDEGlist.FC.RData")

# number of unique transcripts
dim(RG)
#53801    58

##################################################################################
#read in gene info
##################################################################################

# #source("http://bioconductor.org/biocLite.R")
# #biocLite("biomaRt")
# library(biomaRt)
# 
# #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# #us_mart <- useEnsembl(biomart = "ensembl", mirror = "uswest")
# us_mart <- useEnsembl(biomart = "ensembl", mirror = "uswest",dataset = "mmusculus_gene_ensembl")
# 
# #head(listAttributes(human),100)
# #attributes <- listAttributes(us_mart)
# #attributes[grep("exon", attributes$name),]
# 
# genes <- getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","mgi_symbol","external_gene_name","description"), 
#                filter= "ensembl_gene_id",
#                values = rownames(RGcpm), 
#                mart = us_mart)

genes <- read.csv("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/Allmm10Genes.csv")

genes$genelength <- abs(genes$end_position - genes$start_position) 

#remove duplicates:keeps only first entry
genes <- genes[!duplicated(genes$ensembl_gene_id),]

#filter for genes in list
values = rownames(RG)
values <- as.data.frame(values)

tmp <- merge(values, genes, by.x="values", by.y = "ensembl_gene_id", all.x=T )

tmp <- tmp %>% dplyr::select(-X)

#match order: df[match(target, df$name),]
genes <- genes[match(rownames(RG),genes$ensembl_gene_id),]


##########################################################################
#add gene info to DEGlist object
##########################################################################

#add into DEGlist
RG$genes <- genes


###########################################################
# Filter for one count per million 
###########################################################
allsamples <- table(FCsamples$Dam.Treatment.Group,FCsamples$Pup.Treatment.Group,FCsamples$Sex.x)
allsamples

#smallest group is n=6
keep <- rowSums(cpm(RG)>1)>=6
RGcpm <- RG[keep,]
dim(RGcpm)
# 16680    59

#save
save(RG,RGcpm,file="RNAseqDEGlist.FC.RData")
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
#treatment * sex with subjects correlation removed
##################################################################################
# Make litter a random effect, since limma warns "coefficients not estimable" for some litters
# Ref: https://support.bioconductor.org/p/11956/
# Obstacle: Cannot do this properly with surrogtate variables, since there's an error when including litter in null model
# Duplicate correlations alternative for other scenarios:
# https://support.bioconductor.org/p/68916/
# https://support.bioconductor.org/p/110987/
# https://support.bioconductor.org/p/131179/


#make design matrix
FCsamples$Sex <- factor(FCsamples$Sex.x, levels=c("M","F"))

FCsamples$Dam.Treatment.Group <- factor(FCsamples$Dam.Treatment.Group,levels = c("Saline","PolyIC"))

FCsamples$Pup.Treatment.Group <- factor(FCsamples$Pup.Treatment.Group,levels = c("Saline","Treg"))

FCsamples$DamID <- factor(FCsamples$Dam..)

#dam treatment, pup treatment, and sex
FCsamples$Group <- paste(FCsamples$Dam.Treatment.Group,FCsamples$Pup.Treatment.Group,FCsamples$Sex, sep="_")
table(FCsamples$Group)

#set design matrix: full sex * group interactions
design <- model.matrix(~0+Group, data=FCsamples) #if using a 0 intercept must set up contrasts

colnames(design)

#make contrasts
#https://support.bioconductor.org/p/91718/
#test whether the average treatment effect across all cell lines and treatment methods is significantly different from zero
#contrast compares the average of all treatment samples to the average of all control samples.

#https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
#Contrasts within each sex.
#a positive log2-fold-change (logFC) will indicate a gene up-regulated in PolyIC relative to saline, whereas a negative logFC will indicate a gene more highly expressed in saline.
#https://support.bioconductor.org/p/74466/

my.contrasts <- makeContrasts(
  #inidividual comparisons
  M_PolyIC_salvsTreg = GroupPolyIC_Treg_M -GroupPolyIC_Saline_M,
  M_Saline_salvsTreg =  GroupSaline_Treg_M - GroupSaline_Saline_M,
  M_PolyICvsSaline_sal = GroupPolyIC_Saline_M -GroupSaline_Saline_M,
  M_PolyICvsSaline_Treg =  GroupPolyIC_Treg_M - GroupSaline_Treg_M,
  F_PolyIC_salvsTreg = GroupPolyIC_Treg_F -GroupPolyIC_Saline_F,
  F_Saline_salvsTreg =  GroupSaline_Treg_F - GroupSaline_Saline_F,
  F_PolyICvsSaline_sal = GroupPolyIC_Saline_F -GroupSaline_Saline_F,
  F_PolyICvsSaline_Treg =  GroupPolyIC_Treg_F - GroupSaline_Treg_F,
  
  
  #average across sex treatment
  PolyIC_salvsTreg = ((GroupPolyIC_Treg_M + GroupPolyIC_Treg_F)/2 - (GroupPolyIC_Saline_M + GroupPolyIC_Saline_F)/2),
  Saline_salvsTreg =  ((GroupSaline_Treg_M + GroupSaline_Treg_F)/2 - (GroupSaline_Saline_M + GroupSaline_Saline_F)/2),
  PolyICvsSaline_sal = ((GroupPolyIC_Saline_M + GroupPolyIC_Saline_F)/2 - (GroupSaline_Saline_M + GroupSaline_Saline_F)/2),
  PolyICvsSaline_Treg =  ((GroupPolyIC_Treg_M + GroupPolyIC_Treg_F)/2 - (GroupSaline_Treg_M + GroupSaline_Treg_F)/2),
  
  levels=design)
  #interactions: NT vs LPS for ASD vs Typical:
  #ASDNTvsLPS_TypicalNTvsLPS = (GroupASD.LPS - GroupASD.NT) - (GroupTypical.LPS-GroupTypical.NT),
  #interactions: NT vs LPS for PDDNOS vs Typical:
  #ASDNTvsLPS_TypicalNTvsLPS = (GroupPDDNOS.LPS - GroupPDDNOS.NT) - (GroupTypical.LPS-GroupTypical.NT),
  #interactions: NT vs LPS for PDDNOS vs ASD:
  #ASDNTvsLPS_TypicalNTvsLPS = (GroupPDDNOS.LPS - GroupPDDNOS.NT) - (GroupASD.LPS-GroupASD.NT),
  #interactions: NT vs LTA for ASD vs Typical:
 # ASDNTvsLTA_TypicalNTvsLTA = (GroupASD.LTA - GroupASD.NT) - (GroupTypical.LTA-GroupTypical.NT),
  #interactions: NT vs LTA for PDDNOS vs Typical:
 # ASDNTvsLTA_TypicalNTvsLTA = (GroupPDDNOS.LTA - GroupPDDNOS.NT) - (GroupTypical.LTA-GroupTypical.NT),
  #interactions: NT vs LTA for PDDNOS vs ASD:
 # ASDNTvsLTA_TypicalNTvsLTA = (GroupPDDNOS.LTA - GroupPDDNOS.NT) - (GroupASD.LTA-GroupASD.NT),


##################################################################################
#Normalization TMM and Voom
##################################################################################
#The two steps refer to different aspects of normalization. CPM "normalization" accounts for library size differences between samples, and produces normalized values that can be compared on an absolute scale (e.g., for filtering). TMM normalization accounts for composition bias, and computes normalization factors for comparing between libraries on a relative scale. CPM normalization doesn't account for composition bias, and TMM normalization doesn't produce normalized values. Thus, you need both steps in the analysis pipeline. This isn't a problem, as the two steps aren't really redundant.
#TMM normalization for library composition
DGE=calcNormFactors(RGcpm,method =c("TMM")) 

pdf('VoomTrend.pdf',w=6,h=4)
#v=voom(DGE,design,plot=T)


#remove varience due to litter
#corfit <- duplicateCorrelation(v, design, block=FCsamples$DamID)

# fit <- lmFit(v, design)
# 
# fit2 <- contrasts.fit(fit,my.contrasts)
# fit2 <- eBayes(fit2,robust=TRUE )

#https://support.bioconductor.org/p/131179/
#The function is analogous to calling voom followed by duplicateCorrelation and lmFit except for the modified residual df values and residual standard deviation sigma values. This function returns df.residual values that are less than or equal to those from lmFit and sigma values that are greater than or equal to those from lmFit. voomLmFit is more robust to zero counts than calling voom, duplicateCorrelation and lmFit separately and provides more rigorous error rate control.

#lmFit computes coefficients, residual variances and standard errors.

#contrasts.fit converts the coefficients and standard errors to reflect the contrasts rather than the original design matrix, but does not compute t-statistics or p-values.

#eBayes computes t-statistics and p-values from the coefficients and standard errors.


fit <- voomLmFit(DGE, design, block=FCsamples$DamID, plot=TRUE)
fit2 <- contrasts.fit(fit,my.contrasts)
fit2 <- eBayes(fit2,robust=TRUE )

dev.off()

#What is voom doing?

# Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
# A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
# A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
# The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.
# More details at https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29
# 

pdf('PlotSA_VoomTrend.pdf',w=6,h=4)
plotSA(fit2, main="Final model: Mean variance trend")
dev.off()


#log2 CPM
log2cpm <- cpm(DGE, prior.count=2, log=TRUE)


#plot only normalized data
pdf('TMM_VOOM_NormalizedLogCPM.pdf',w=30,h=10)
boxplot(log2cpm, xlab="", ylab="Log2 counts per million",las=2,main="TMM & Voom transformed logCPM")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
dev.off()

#save
save(RG,RGcpm,DGE,log2cpm, FCsamples,design,my.contrasts,fit,fit2, file="RNAseqDEGlist.FC.RData")

##########################################################################
#MAA plots
##########################################################################

#MAA plot for each sample
FCsamples$Name <- paste(FCsamples$Mouse.Number,FCsamples$Group,sep=".")
FCsamples$Name <- gsub("\\."," ", FCsamples$Name)

C <- unique(FCsamples$Name)
#Mean-Difference Plot of Expression Data
#8 per page
pdf('Samples_MAAplots.pdf')
par(mfrow=c(4,2),mar=c(5,6,4,2)+0.1)
for (i in 1:length(C)) {
  print(paste(i))
  name <- C[i]
  #pdf(file = paste(name,"MAAplot.pdf", sep ="_"), wi = 3, he = 3)
  #MAA: should center on zero
  plotMD(DGE,column=i,
         ylab = "Expression log-ratio\n(this sample vs others)")
  abline(h=0, col="red", lty=2, lwd=2)
  #dev.off()
  
}

dev.off()

##################################################################################
#MDS plots
##################################################################################

#plot MDS for group
FCsamples$Group <- as.factor(FCsamples$Group)

pdf(file = "MDSplot_group.pdf", wi = 12, he = 10, useDingbats=F)

levels(FCsamples$Group)
col.cell2 <- c("blue","orange","purple","red","pink","yellow","green","brown")[FCsamples$Group]
data.frame(FCsamples$Group,col.cell2)


# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(DGE,col=col.cell2)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1.5, 1,
       fill=c("blue","orange","purple","red","pink","yellow","green","brown"),
       legend=unique(FCsamples$Group),
       cex = 0.8)
# Add a title
title("Group MDS Plot")

dev.off()





#plot MDS for sex
pdf(file = "MDSplot_Sex.pdf", wi = 12, he = 10, useDingbats=F)

levels(FCsamples$Sex)
col.cell2 <- c("blue","orange")[FCsamples$Sex]
data.frame(FCsamples$Sex,col.cell2)


# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(DGE,col=col.cell2)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1.5, 1,
       fill=c("blue","orange"),
       legend=levels(FCsamples$Sex),
       cex = 0.8)
# Add a title
title("Sex MDS Plot")

dev.off()


#by dam Treatment:
pdf(file = "MDSplot_DamTreatment.pdf", wi = 12, he = 10, useDingbats=F)

#par(mfrow=c(1,2))
#plot MDS for treatment
levels(FCsamples$Dam.Treatment.Group)
col.cell <- c("blue","orange")[FCsamples$Dam.Treatment.Group]
data.frame(FCsamples$Dam.Treatment.Group,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(DGE,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1.5, 1,
       fill=c("blue","orange"),
       legend=levels(FCsamples$Dam.Treatment.Group),
       cex = 0.8)
# Add a title
title("Dam Treatment MDS Plot")
dev.off()

#pup treatment
pdf(file = "MDSplot_PupTreatment.pdf", wi = 12, he = 10, useDingbats=F)

#par(mfrow=c(1,2))
#plot MDS for treatment
levels(FCsamples$Pup.Treatment.Group)
col.cell <- c("blue","purple")[FCsamples$Pup.Treatment.Group]
data.frame(FCsamples$Pup.Treatment.Group,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(DGE,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1.5, 1,
       fill=c("blue","purple"),
       legend=levels(FCsamples$Pup.Treatment.Group),
       cex = 0.8)
# Add a title
title("Pup Treatment MDS Plot")
dev.off()

#library size

size <- as.data.frame(RGcpm$samples)
tmp <- merge(FCsamples, size, by.x = "filenames",by.y = "files")

tmp <- tmp %>% arrange(lib.size)

colfunc <- colorRampPalette(c("pink", "red"))
colfunc(length(tmp$lib.size))

pdf(file = "MDSplot_librarysize.pdf", wi = 12, he = 10, useDingbats=F)

#par(mfrow=c(1,2))
#plot MDS for treatment

col.cell <- colfunc(length(tmp$lib.size))
data.frame(tmp$lib.size,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(DGE,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(3, 1,
       fill=col.cell,
       legend=unique(tmp$lib.size),
       cex = 0.8)
# Add a title
title("Library Size MDS Plot")
dev.off()


#litter

library(RColorBrewer)
n <- length(unique(FCsamples$DamID))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

pdf(file = "MDSplot_Litter.pdf", wi = 12, he = 10, useDingbats=F)


#colors
#par(mfrow=c(1,2))
#plot MDS for treatment
levels(FCsamples$DamID)
col.cell <- col_vector[FCsamples$DamID]
data.frame(FCsamples$DamID,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(DGE,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(3, 1,
       fill=col_vector,
       legend=levels(FCsamples$DamID),
       cex = 0.8)
# Add a title
title("Litter MDS Plot")
dev.off()


##################################################################################
#Average CPM for each condition
##################################################################################

#get log2CPM counts from voom and put in dataframe:
library(plotrix)
library(xlsx)
#average log2 CPM and sem
countdf <- as.data.frame(log2cpm)
countdf$GeneID <- rownames(log2cpm)
#targets <- v$targets
#targets$sample <- gsub(".out","",rownames(targets))

#DF <- merge(countdf,comp_out, by.x ="GeneID",by.y="genes")

#summarize 
countdf2 <- countdf %>% group_by(GeneID) %>% gather(Sample.ID,log2CPM, 1:(ncol(countdf)-1)) 
countdf2 <- as.data.frame(countdf2)
countdf3 <-merge(countdf2,FCsamples,by="Sample.ID")

GeneSummary <- countdf3 %>% group_by(GeneID,Group) %>% 
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
 # dplyr::select(GeneID,Group,meanlog2CPM) %>%
  spread(Group,meanlog2CPM)

#add in gene symbol information
geneids <- DGE$genes
GeneSummary2 <- merge(GeneSummary,geneids, by.x="GeneID",by.y="ensembl_gene_id")
 
write.csv(GeneSummary2,file= "MeanLog2CPM_pergroup_GeneExpression.csv")

#log2CPM values for individuals
countdf3$Individual <- paste(countdf3$Group,countdf3$Mouse.Number, sep=".")
countdf3$Individual <- factor(countdf3$Individual)

IndGeneSummary <- countdf3 %>% dplyr::group_by(GeneID,Individual) %>% 
  summarise(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  # dplyr::select(GeneID,Group,meanlog2CPM) %>%
  spread(Individual,meanlog2CPM)

#add in gene symbol information
IndGeneSummary2 <- merge(IndGeneSummary,geneids ,by.x="GeneID",by.y="ensembl_gene_id")

write.csv(IndGeneSummary2,file= "IndividualSamplesLog2CPM_GeneExpression.csv")


#go back to CPM from log2: 2^log2CPM = CPM
GeneSummaryCPM <- countdf3 %>% group_by(GeneID,Group) %>% 
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  mutate(meanCPM = 2^meanlog2CPM) %>%
  dplyr::select(GeneID,Group,meanCPM) %>%
  spread(Group,meanCPM)

#add in gene symbol information
GeneSummaryCPM2 <- merge(GeneSummaryCPM,geneids,by.x="GeneID",by.y="ensembl_gene_id")

write.csv(GeneSummaryCPM2,file= "CPM_GeneExpression.csv")

save(GeneSummary2,GeneSummaryCPM2,file="AverageCPMandLog2CPM.RData")

##################################################################################
#DE analysis:
##################################################################################
str(fit2)
#https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#differential-expression-analysis
#more weight to fold-changes in the ranking
#treat computes empirical Bayes moderated-t p-values relative to a minimum meaningful fold-change threshold. 
#Instead of testing for genes that have true log-fold-changes different from zero, it tests whether the true 
#log2-fold-change is greater than lfc in absolute value (McCarthy and Smyth, 2009). 
#In other words, it uses an interval null hypothesis, where the interval is [-lfc,lfc]. 
#When the number of DE genes is large, treat is often useful for giving preference to larger fold-changes and 
#for prioritizing genes that are biologically important. 
#treat is concerned with p-values rather than posterior odds, so it does not compute the B-statistic lods. 
#The idea of thresholding doesn't apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed. 
#When lfc=0, treat is identical to eBayes, except that F-statistics and B-statistics are not computed. 
#The lfc threshold is usually chosen relatively small, because significantly DE genes must all have fold changes substantially greater than the testing threshold. 
#Typical values for lfc are log2(1.1), log2(2) or log2(2). The top genes chosen by treat can be examined using topTreat.

# Treat relative to a ~2 fold-change
#requires genes to have a log-FC that is significantly greater than 1 (equivalent to a 2-fold difference between cell types on the original scale).
tfit <- treat(fit2,lfc=log2(1))
dt <- decideTests(tfit,adjust.method="fdr", method="separate")
sum <- summary(dt)
sum
write.csv(sum,"SummaryCount_DEGs.csv")

dt <- as.data.frame(dt)

write.fit(tfit, file="results.txt")

#######################################################################
#get out DE lists for each contrast:
#######################################################################
#if memory error: https://www.biostars.org/p/432389/
options(future.globals.maxSize = 4000 * 1024^5)
library(calibrate)

comparisons=(coef(tfit))
comparisons=colnames(comparisons)

comp_out <- as.data.frame(DGE$genes$ensembl_gene_id)
names(comp_out) <- c("GeneID")
nrowkeep <- nrow(comp_out)

SumTableOut <- NULL


for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  print(comp)
  #make comparisons 
  
  tmp=topTreat(tfit,coef=i,number=nrowkeep,adjust.method="fdr")
  #nrow(tmp[(tmp$adj.P.Val<0.05),]) # number of DE genes
  
  #LogFC values:https://support.bioconductor.org/p/82478/
  tmp$direction <- c("none")
  tmp$direction[which(tmp$logFC > 0)] = c("Increase")
  tmp$direction[which(tmp$logFC < 0)] = c("Decrease")
  
  tmp$significance <- c("nonDE")
  tmp$significance[which(tmp$adj.P.Val <0.05)] = c("DE")
  
  #summary counts table based on Ensemble Gene ID counts:
  SumTable <- table(tmp$significance,tmp$direction)
  SumTable <- as.data.frame(SumTable)
  SumTable$comparison <- paste(comp)
  SumTableOut <- rbind(SumTable,SumTableOut)
  
  #get geneids  
  tmp$GeneID <- rownames(tmp)
  
  #gene gene names and expression levels
  tmp2 <- tmp
  
  tmp2$comparison <- paste(comp)
  
  write.csv(tmp2,file = paste(comp,"_DEgenes.csv"))
  
  #save to output:
  #merge <- merge(comp_out,tmp2, by= "GeneID")
 # merge2 <- tmp2 %>% dplyr::select(ensembl_gene_id,logFC,t,P.Value,adj.P.Val,direction,significance)
  #colnames(merge2) <- paste(colnames(merge2),comp,sep=".")
  #colnames(merge2)[1] <- c("GeneID")
  #comp_out <- merge(comp_out, merge2, by="GeneID")

  
  #data for plot with gene names:
  genenames <- tmp2 %>% dplyr::select(adj.P.Val,logFC,mgi_symbol) %>% distinct()
  
  #names for plots
  plotname <- gsub("\\."," ",comp)
  plotname <- gsub("vs"," vs ",plotname)
  
  #volcano plot
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  
  with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(genenames, logFC < -2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(genenames, logFC > 2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-2,2), col = "black", lty = 2, lwd = 1)
  
  #Label points with the textxy function from the calibrate plot

  with(subset(genenames, adj.P.Val<0.05 & abs(logFC)>2), textxy(logFC, -log10(adj.P.Val), labs=mgi_symbol, cex=.4))
  
  dev.off()
  
}


write.csv(SumTableOut,"SummaryTableDEgenes.csv")

#master outfile to get log2CPM values
#mout <- merge(comp_out, GeneSummary2, by="GeneID")
#write.csv(mout,"AllDEG_AllConditions_log2CPM.csv")

save(tfit,dt,SumTableOut,GeneSummaryCPM2,GeneSummary2,countdf,RG,genes,file = "RNAseqDEGlist.FC.RData")
#load("RNAseqDEGlist.FC.RData")


#######################################################################
#plot only FC_M_PolyICvsSaline_Treg_Volcano
#######################################################################
i="M_PolyICvsSaline_Treg"
tmp=topTreat(tfit,coef=i,number=nrowkeep,adjust.method="fdr")
#nrow(tmp[(tmp$adj.P.Val<0.05),]) # number of DE genes

#LogFC values:https://support.bioconductor.org/p/82478/
tmp$direction <- c("none")
tmp$direction[which(tmp$logFC > 0)] = c("Increase")
tmp$direction[which(tmp$logFC < 0)] = c("Decrease")

tmp$significance <- c("nonDE")
tmp$significance[which(tmp$adj.P.Val <0.05)] = c("DE")


#get geneids  
tmp$GeneID <- rownames(tmp)

#gene gene names and expression levels
tmp2 <- tmp

#data for plot with gene names:
genenames <- tmp2 %>% dplyr::select(adj.P.Val,logFC,mgi_symbol) %>% distinct()

comp <- i
plotname <- i
#names for plots
plotname <- gsub("\\."," ",comp)
plotname <- gsub("vs"," vs ",plotname)

#volcano plot
pdf(file = paste(comp,"_VolcanoFC.pdf", sep=""), wi = 6, he = 6, useDingbats=F)

with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))

#color points red when sig and log2 FC > .5 and blue if log2 FC < -.5 
with(subset(genenames, logFC < -.1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
with(subset(genenames, logFC > .1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))

#add lines
abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
#abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)

#Label points with the textxy function from the calibrate plot

with(subset(genenames, adj.P.Val<0.05 & abs(logFC)>.5), textxy(logFC, -log10(adj.P.Val), labs=mgi_symbol, cex=.6))

dev.off()

#######################################################################
#Venn overlaps for M_PolyICvsSaline_Treg and PolyICvsSaline_Treg
#######################################################################

#venns:
#The number of genes that are not DE in either comparison are marked in the bottom-right.

pdf("Venn_PolyICvsSaline_Treg_down.pdf",height=10,width = 10)
vennDiagram(dt[,c("PolyICvsSaline_Treg","M_PolyICvsSaline_Treg")],
            circle.col=c("salmon","purple","grey"), include="down",
            names=c("All PolyICvsSaline_Treg","Male PolyICvsSaline_Treg"))
dev.off()

pdf("Venn_PolyICvsSaline_Treg_up.pdf",height=10,width = 10)
vennDiagram(dt[,c("PolyICvsSaline_Treg","M_PolyICvsSaline_Treg")],
            circle.col=c("salmon","purple","grey"), include="up",
            names=c("All PolyICvsSaline_Treg","Male PolyICvsSaline_Treg"))
dev.off()


#######euler#######
library(eulerr)
t <- euler(c("A" = 18, "B" = 1149, "A&B" = 319))
pdf("euler_PolyICvsSaline_Treg_up.pdf",height=10,width = 10)
plot(t, quantities = TRUE, edges = FALSE,  fills = c("cadetblue4", "darkgoldenrod1"), labels = c("All PolyICvsSaline_Treg","Male PolyICvsSaline_Treg")) 
dev.off()

t <- euler(c("A" = 17, "B" = 864, "A&B" = 264))
pdf("euler_PolyICvsSaline_Treg_down.pdf",height=10,width = 10)
plot(t, quantities = TRUE, edges = FALSE,  fills = c("cadetblue4", "darkgoldenrod1"), labels = c("All PolyICvsSaline_Treg","Male PolyICvsSaline_Treg")) 
dev.off()


#list of overlapping genes #######
#Genes that are DE in multiple comparisons can be extracted using the results from decideTests, where 0s represent genes that are not DE, 1s represent genes that are up-regulated, and -1s represent genes that are down-regulated. 
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#useful-graphical-representations-of-differential-expression-results

#up in both conditions
de.common.up <- which(dt[,c("PolyICvsSaline_Treg")]==1 & dt[,"M_PolyICvsSaline_Treg"]==1)
de.common.up <- tfit$genes$mgi_symbol[de.common.up]
de.common.up <- unique(de.common.up)
length(de.common.up)

#down in both conditions
de.common.down <- which(dt[,c("PolyICvsSaline_Treg")]==-1 & dt[,"M_PolyICvsSaline_Treg"]==-1)
de.common.down <- tfit$genes$mgi_symbol[de.common.down]
de.common.down <- unique(de.common.down)
length(de.common.down)

#only up in All condition
de.allonly.up <- which(dt[,c("PolyICvsSaline_Treg")]==1 & dt[,"M_PolyICvsSaline_Treg"]==0)
de.allonly.up <- tfit$genes$mgi_symbol[de.allonly.up]
de.allonly.up <- unique(de.allonly.up)
length(de.allonly.up)

#only down in All condition (some repeat symbols so 176 instead of 181)
de.allonly.down <- which(dt[,c("PolyICvsSaline_Treg")]==-1 & dt[,"M_PolyICvsSaline_Treg"]==0)
de.allonly.down <- tfit$genes$mgi_symbol[de.allonly.down]
de.allonly.down <- unique(de.allonly.down)
length(de.allonly.down)

#only up in Male condition
de.Maleonly.up <- which(dt[,c("PolyICvsSaline_Treg")]==0 & dt[,"M_PolyICvsSaline_Treg"]==1)
de.Maleonly.up <- tfit$genes$mgi_symbol[de.Maleonly.up]
de.Maleonly.up <- unique(de.Maleonly.up)
length(de.Maleonly.up)

#only down in Male condition 
de.Maleonly.up <- which(dt[,c("PolyICvsSaline_Treg")]==0 & dt[,"M_PolyICvsSaline_Treg"]==-1)
de.Maleonly.down <- tfit$genes$mgi_symbol[de.Maleonly.down]
de.Maleonly.down <- unique(de.Maleonly.down)
length(de.Maleonly.down)


#combine into list
x = list(de.common.up,
         de.common.down,
         de.allonly.up,
         de.allonly.down,
         de.Maleonly.up,
         de.Maleonly.up)


rapply(x, length, how="list")


##################################################################
#Fisher's exact test for list overlaps LPS
##################################################################
# library(GeneOverlap)
# #https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf
# 
# #background = all genes expressed in the experiment
# BG <- length(GeneSummaryCPM$GeneID)
# 
# #overlap list all vs all
# go.obj <-  newGOM(x,genome.size=BG)
# print(go.obj)
# drawHeatmap(go.obj)
# 
# #get p values
# options(digits=10)
# pval <- getMatrix(go.obj, name="pval")
# 
# pval <- as.data.frame(pval)
# pval_df <- pval %>% mutate(List1 = rownames(pval)) %>%
#   gather(List2,FisherExactpvalue,1:3)
# 
# #get odds ratio
# OR <- getMatrix(go.obj, "odds.ratio")
# OR <- as.data.frame(OR)
# 
# OR_df <- OR %>% mutate(List1 = rownames(OR)) %>%
#   gather(List2,OddsRatio,1:3)
# 
# Overlaps <- merge(pval_df,OR_df,by=c("List1","List2"))
# 
# #remove List1 == List2
# Overlaps <- Overlaps[Overlaps$List1!=Overlaps$List2,]
# 
# #adjust by fdr
# Overlaps$fdr <- p.adjust(as.numeric(Overlaps$FisherExactpvalue),method="fdr")
# 
# write.csv(Overlaps,"LPSVennListOverlapsFishersExactTest.csv")
# 
# ##################################################################
# #plot:
# ##################################################################
# 
# library(nVennR)
# 
# #get overlapping genes:
# myV <- plotVenn(x,showPlot = T)
# listregions <- listVennRegions(myV)
# 
# #fix names:
# namesregionlist <- names(listregions)
# 
# library(qdapRegex)
# names(listregions) <- rm_between(namesregionlist, "(", ")", extract=TRUE)
# 
# library(tidyverse)
# l_tib <- listregions %>% 
#   #unlist(recursive = T) %>% 
#   enframe() %>% 
#   unnest()
# 
# l_tib <- as.data.frame(l_tib)
# names(l_tib) <- c("LPSVennGroup","GeneID")
# l_tib$GeneID <- as.character(l_tib$GeneID)
# 
# mout2 <- merge(mout,l_tib, by = "GeneID",all.x=T)
# 
# write.csv(mout2,"AllDEG_AllConditions.csv")
# 



#######################################################################
# GO enrichment of each list using ClusterCompare
#######################################################################
#library(devtools)
#devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
library("clusterProfiler")
#The input for geneCluster parameter should be a named list of gene IDs. 

#all
list1.up <- which(dt[,c("PolyICvsSaline_Treg")]==1)
list1.up <- unique(tfit$genes$mgi_symbol[list1.up])

list1.down <- which(dt[,c("PolyICvsSaline_Treg")]==-1)
list1.down <- unique(tfit$genes$mgi_symbol[list1.down])


#male only 
list2.up <- which(dt[,c("M_PolyICvsSaline_Treg")]==1)
list2.up <- unique(tfit$genes$mgi_symbol[list2.up])

list2.down <- which(dt[,c("M_PolyICvsSaline_Treg")]==-1)
list2.down <- unique(tfit$genes$mgi_symbol[list2.down])

#combine into list
x2 = list(list1.up,
          list1.down,
          list2.up,
          list2.down)



names(x2) <- c("DamPolyIC>Saline_PupTreg", "DamPolyIC<Saline_PupTreg", "DamPolyIC>Saline_MalePupTreg", "DamPolyIC<Saline_MalePupTreg")

rapply(x2, length, how="list")

#######################################################################
#for KEGG convert to entreze
#######################################################################

#BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
library(xlsx)

listDEGs2  <- NULL

for (i in 1:length(names(x2))){

  e <- bitr(x2[[i]],fromType='SYMBOL',toType="ENTREZID", OrgDb="org.Mm.eg.db")
  head(e$ENTREZID)
  
  #save to output list:
  listDEGs2[[i]] <- unique(e$ENTREZID)
  
}

names(listDEGs2) <- names(x2)

#universe = all detected genes in experiment
universe_enterz <- bitr(GeneSummaryCPM$GeneID,fromType='ENSEMBL',toType="ENTREZID", OrgDb="org.Mm.eg.db")
length(universe_enterz$ENTREZID)

#clusters
ck <- compareCluster(geneCluster = listDEGs2,
                     organism="mmu", 
                     universe      = unique(universe_enterz$ENTREZID),
                     fun = "enrichKEGG",
                     qvalueCutoff  = 0.05)




pdf('KeggEnrichment_allDEGs.pdf',w=9,h=11)
d <- dotplot(ck,showCategory = 20)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

ck <- as.data.frame(ck)
write.csv(ck,file="KEGGenrichments_allDEGs.csv")


#######################################################################
#GO term enrichment
#######################################################################

#######################################################################
#cellular component
#######################################################################

cc <- compareCluster(geneCluster = x2,
                     universe      = unique(GeneSummary2$mgi_symbol),
                     OrgDb         = org.Mm.eg.db,
                     ont           = "CC",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     fun = "enrichGO")

#cc <- compareCluster(geneCluster = ck, use_internal_data = TRUE, fun = "enrichKEGG")

head(cc)
dim(cc)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
cc2 <- simplify(cc, cutoff=0.7, by="p.adjust", select_fun=min)
dim(cc2)


pdf('GO_CCEnrichment_allDEGs.pdf',w=9,h=11)
d <- dotplot(cc2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

cc2 <- as.data.frame(cc2)
write.csv(cc2,file="GO_CCenrichments_allDEGs.csv")

#######################################################################
#biological process
#######################################################################

bp <-  compareCluster(geneCluster = x2,
                      universe      = unique(GeneSummary2$mgi_symbol),
                      OrgDb         = org.Mm.eg.db,
                      ont           = "BP",
                      keyType       = 'SYMBOL',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      fun = "enrichGO")
head(bp)
dim(bp)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
dim(bp2)


pdf('GO_BPEnrichment_allDEGs.pdf',w=9,h=14)
d <- dotplot(bp2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

bp2 <- as.data.frame(bp2)
write.csv(bp2,file="GO_BPenrichments_allDEGs.csv")

#######################################################################
#molecular function
#######################################################################

mf <- compareCluster(geneCluster = x2,
                     universe      = unique(GeneSummary2$mgi_symbol),
                     OrgDb         = org.Mm.eg.db,
                     ont           = "mf",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     fun = "enrichGO")
head(mf)
dim(mf)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
mf2 <- simplify(mf, cutoff=0.7, by="p.adjust", select_fun=min)
dim(mf2)


pdf('GO_MFEnrichment_allDEGs.pdf',w=9,h=14)
d <- dotplot(mf2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

mf2 <- as.data.frame(mf2)
write.csv(mf2,file="GO_MFenrichments_allDEGs.csv")


################################################################################################
#GLIMMA interactive plot building
################################################################################################
#BiocManager::install("Glimma", version = "3.8")
#http://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf
#https://bioconductor.org/packages/devel/bioc/vignettes/Glimma/inst/doc/limma_edger.html
library(Glimma)


#MDS plots
samplefiles <- DGE$samples$files

#match ordre of df to v
group.df <- FCsamples[match(samplefiles, FCsamples$filenames),]
group.df <- group.df %>% dplyr::select(Dam.Treatment.Group,Pup.Treatment.Group,Sex,DamID,Cage,Group)


#writes out html file:
glMDSPlot(DGE, groups=group.df)

#DE genes
#tfit <- treat(fit2,lfc=log2(1))
#dt <- decideTests(tfit,adjust.method="fdr", method="separate")

summary(dt)

#glMDPlot(tfit, counts=v$E,transform=FALSE,status = results,coef=1,
#         anno=tfit$genes, groups = group.df$Group,
#          main=paste("MD plot:",colnames(results)[1]),
 #       display.columns=c("hgnc_symbol", "ensembl_gene_id"),
  #       folder=paste(colnames(results)[1]),side.main="hgnc_symbol")

group.df$Group <- as.factor(group.df$Group)
sample.cols <- c("purple","orange","green","blue","pink",
                 "brown","red","black")[group.df$Group]

for (COEF in 1:12) {
  glMDPlot(tfit, dge = DGE,
           anno=tfit$genes, groups = group.df$Group,
           coef=COEF, main=colnames(tfit)[COEF],
           side.ylab="Log2CPM Expression", side.main="mgi_symbol",
           folder="glimma_results",status = dt,
           sample.cols =sample.cols,
           html = paste("MD-Plot",colnames(tfit)[COEF]))
}

 

################################################################################################
#Heatmaps
################################################################################################

#heatmap for DEGs from both sexes; average Log2CPM per group:
library(pheatmap)

DEGs_nonsex1 <- x2$`DamPolyIC>Saline_PupTreg`

DEGs_nonsex2 <-x2$`DamPolyIC<Saline_PupTreg`

tmp1 <- subset(GeneSummary2, mgi_symbol %in% DEGs_nonsex1)  
tmp1$DEG <- c("DamPolyIC>Saline_PupTreg")

tmp2 <- subset(GeneSummary2, mgi_symbol %in% DEGs_nonsex2)  
tmp2$DEG <- c("DamPolyIC<Saline_PupTreg")

heatDF <- rbind(tmp1,tmp2)

heatDF2 <- heatDF[,2:9]

rownames(heatDF2) <- heatDF$mgi_symbol

#matrix
matrix <- as.matrix(heatDF2)


#define row groups:
annotation_row <- heatDF$DEG
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- c("DEG Condition")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)

#column annotations:
annotation_col <- as.data.frame(colnames(matrix))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Dam Treatment", "Pup Treatment", "Pup Sex"),sep="_") 


rownames(annotation_col) = colnames(matrix)

head(annotation_col)

# Specify colors
Var1        <- c("royalblue3","orange2")
names(Var1) <- unique(annotation_col$`Dam Treatment`)

Var2 <- c("gold","darkred")
names(Var2) <- unique(annotation_col$`Pup Treatment`)

Var3 <- c("darkorchid4","darkorange")
names(Var3) <- unique(annotation_col$`Pup Sex`)

Var4 <- c("violet","darkgreen")
names(Var4) <- unique(annotation_row$`DEG Condition`)


#combined for heatmap labels
anno_colors <- list(`Dam Treatment`=Var1,
                    `Pup Treatment`=Var2,
                    `Pup Sex` = Var3,
                    `DEG Condition` = Var4)


#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_DEGs_bothsexes.pdf", wi = 8, he = 50)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                  annotation_row = annotation_row,
                   show_rownames = T,
                  annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()

pdf(file = "Heatmap_DEGs_bothsexes_nonames.pdf", wi = 8, he = 16)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()

########################################################################
#heatmap for DEGs from both sexes; individual animal Log2CPM:
########################################################################


DEGs_nonsex1 <- x2$`DamPolyIC>Saline_PupTreg`

DEGs_nonsex2 <-x2$`DamPolyIC<Saline_PupTreg`

tmp1 <- subset(IndGeneSummary2, mgi_symbol %in% DEGs_nonsex1)  
tmp1$DEG <- c("DamPolyIC>Saline_PupTreg")

tmp2 <- subset(IndGeneSummary2, mgi_symbol %in% DEGs_nonsex2)  
tmp2$DEG <- c("DamPolyIC<Saline_PupTreg")

heatDF <- rbind(tmp1,tmp2)

heatDF2 <- heatDF[,2:59]

rownames(heatDF2) <- heatDF$mgi_symbol

#matrix
matrix <- as.matrix(heatDF2)


#define row groups:
annotation_row <- heatDF$DEG
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- c("DEG Condition")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)

#column annotations:
annotation_col <- as.data.frame(colnames(matrix))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Dam Treatment", "Pup Treatment", "PupSex"),sep="_") 

annotation_col <- tidyr::separate(annotation_col,PupSex,into=c("Pup Sex","Mouse.Number"),sep="\\.") 

#add in litter info
tmp <- merge(annotation_col, FCsamples, by="Mouse.Number")

annotation_col <- tmp %>% dplyr::select(`Dam Treatment` , `Pup Treatment`, `Pup Sex`,DamID)

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

# Specify colors
Var1        <- c("royalblue3","orange2")
names(Var1) <- unique(annotation_col$`Dam Treatment`)

Var2 <- c("gold","darkred")
names(Var2) <- unique(annotation_col$`Pup Treatment`)

Var3 <- c("darkorchid4","darkorange")
names(Var3) <- unique(annotation_col$`Pup Sex`)

Var4 <- c("violet","darkgreen")
names(Var4) <- unique(annotation_row$`DEG Condition`)

library(RColorBrewer)
n <- length(unique(FCsamples$DamID))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Var5 <- col_vector
names(Var5) <- unique(annotation_col$DamID)

#combined for heatmap labels
anno_colors <- list(`Dam Treatment`=Var1,
                    `Pup Treatment`=Var2,
                    `Pup Sex` = Var3,
                    `DEG Condition` = Var4,
                     Litter = Var5)


#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_DEGs_bothsexes_IndividualAnimals.pdf", wi = 16, he = 50)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = T,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()

pdf(file = "Heatmap_DEGs_bothsexes_nonames_IndividualAnimals.pdf", wi = 16, he = 16)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()


################################################################################################
# Male DEGs only
################################################################################################
#heatmap for DEGs from males only; average Log2CPM per group:

DEGs_nonsex1 <- x2$`DamPolyIC>Saline_MalePupTreg`

DEGs_nonsex2 <-x2$`DamPolyIC<Saline_MalePupTreg`

tmp1 <- subset(GeneSummary2, mgi_symbol %in% DEGs_nonsex1)  
tmp1$DEG <- c("DamPolyIC>Saline_MalePupTreg")

tmp2 <- subset(GeneSummary2, mgi_symbol %in% DEGs_nonsex2)  
tmp2$DEG <- c("DamPolyIC<Saline_MalePupTreg")

heatDF <- rbind(tmp1,tmp2)

heatDF2 <- heatDF[,2:9]

rownames(heatDF2) <- heatDF$mgi_symbol

#matrix
matrix <- as.matrix(heatDF2)


#define row groups:
annotation_row <- heatDF$DEG
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- c("DEG Condition")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)

#column annotations:
annotation_col <- as.data.frame(colnames(matrix))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Dam Treatment", "Pup Treatment", "Pup Sex"),sep="_") 


rownames(annotation_col) = colnames(matrix)

head(annotation_col)

# Specify colors
Var1        <- c("royalblue3","orange2")
names(Var1) <- unique(annotation_col$`Dam Treatment`)

Var2 <- c("gold","darkred")
names(Var2) <- unique(annotation_col$`Pup Treatment`)

Var3 <- c("darkorchid4","darkorange")
names(Var3) <- unique(annotation_col$`Pup Sex`)

Var4 <- c("violet","darkgreen")
names(Var4) <- unique(annotation_row$`DEG Condition`)


#combined for heatmap labels
anno_colors <- list(`Dam Treatment`=Var1,
                    `Pup Treatment`=Var2,
                    `Pup Sex` = Var3,
                    `DEG Condition` = Var4)


#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_DEGs_Males.pdf", wi = 8, he = 50)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = T,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()

pdf(file = "Heatmap_DEGs_Males_nonames.pdf", wi = 8, he = 16)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()

########################################################################
#heatmap for DEGs males only; individual animal Log2CPM:
########################################################################

DEGs_nonsex1 <- x2$`DamPolyIC>Saline_MalePupTreg`

DEGs_nonsex2 <-x2$`DamPolyIC<Saline_MalePupTreg`


tmp1 <- subset(IndGeneSummary2, mgi_symbol %in% DEGs_nonsex1)  
tmp1$DEG <- c("DamPolyIC>Saline_MalePupTreg")

tmp2 <- subset(IndGeneSummary2, mgi_symbol %in% DEGs_nonsex2)  
tmp2$DEG <- c("DamPolyIC<Saline_MalePupTreg")

heatDF <- rbind(tmp1,tmp2)

heatDF2 <- heatDF[,2:59]

rownames(heatDF2) <- heatDF$mgi_symbol

#matrix
matrix <- as.matrix(heatDF2)


#define row groups:
annotation_row <- heatDF$DEG
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- c("DEG Condition")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)

#column annotations:
annotation_col <- as.data.frame(colnames(matrix))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Dam Treatment", "Pup Treatment", "PupSex"),sep="_") 

annotation_col <- tidyr::separate(annotation_col,PupSex,into=c("Pup Sex","Mouse.Number"),sep="\\.") 

#add in litter info
tmp <- merge(annotation_col, FCsamples, by="Mouse.Number")

annotation_col <- tmp %>% dplyr::select(`Dam Treatment` , `Pup Treatment`, `Pup Sex`,DamID)

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

# Specify colors
Var1        <- c("royalblue3","orange2")
names(Var1) <- unique(annotation_col$`Dam Treatment`)

Var2 <- c("gold","darkred")
names(Var2) <- unique(annotation_col$`Pup Treatment`)

Var3 <- c("darkorchid4","darkorange")
names(Var3) <- unique(annotation_col$`Pup Sex`)

Var4 <- c("violet","darkgreen")
names(Var4) <- unique(annotation_row$`DEG Condition`)

Var5 <- col_vector
names(Var5) <- unique(annotation_col$DamID)

#combined for heatmap labels
anno_colors <- list(`Dam Treatment`=Var1,
                    `Pup Treatment`=Var2,
                    `Pup Sex` = Var3,
                    `DEG Condition` = Var4,
                    Litter = Var5)



#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_DEGs_Males_IndividualAnimals.pdf", wi = 16, he = 50)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = T,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()

pdf(file = "Heatmap_DEGs_Males_nonames_IndividualAnimals.pdf", wi = 16, he = 16)
pheatmap::pheatmap(matrix, 
                   cluster_row = T,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                   fontsize = 16,
                   fontsize_row=6, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score log2CPM mRNA Expression")

dev.off()


