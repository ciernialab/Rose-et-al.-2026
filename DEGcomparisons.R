## code for comparing DEGs across brain regions
## AVC Oct 2022


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(readxl)
library(tidyverse)

#load DEGs from Cerebellum:
path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/Cerebellum/sexlitter"
setwd(path)

#load CB example data to get gene info:
CB <- read.csv("M_PolyIC_salvsTreg _DEgenes.csv")
CB_genes <- as.data.frame(CB[,3:10])
#CB_genes_data <- as.data.frame(CB[,c(11,13,14,15,16,17)])
#colnames(CB_genes_data) <- paste(colnames(CB_genes_data), unique(CB$comparison), sep="_" )


out.file <- data.frame(CB_genes$ensembl_gene_id)
out.file <- out.file[complete.cases(out.file),]

out.file <- data.frame(CB_genes$ensembl_gene_id)
names(out.file) <- c("ensembl_gene_id")


file.names <- dir(path, pattern ='*_DEgenes.csv')
file.names


for(i in 1:length(file.names)){
  file <- read.csv(file.names[i])
  
  #exract stats and rename columns to include comparison performed
  print(unique(file$comparison))
  
  CB_genes_data <- as.data.frame(file[,c(11,13,14,15,16,17)])
  colnames(CB_genes_data) <- paste(colnames(CB_genes_data), unique(file$comparison), sep="_" )
  CB_genes_data$ensembl_gene_id <- file$ensembl_gene_id
  CB_genes_data <- CB_genes_data[complete.cases(CB_genes_data),]
    
  #merge into datafram by ensembl id
  out.file <- merge(out.file, CB_genes_data, by="ensembl_gene_id")
  
  print("done")

}

#add in log2CPM per animal
CB_CPM <- read.csv("IndividualSamplesLog2CPM_GeneExpression.csv")
CB_CPM$X.1 <- NULL

CBmaster <- merge(out.file, CB_CPM, by.x="ensembl_gene_id",by.y="GeneID")

#write to file
write.csv(CBmaster, file = "CB_masterDEGs_log2CPM.csv")

#DEGs

#fix column names to contain M, F or All for comparisons
CBmaster2 <- CBmaster
colnames(CBmaster)<- gsub("_PolyIC_","_All_PolyIC_", colnames(CBmaster))
colnames(CBmaster)<- gsub("F_All_PolyIC_","F_PolyIC_", colnames(CBmaster))
colnames(CBmaster) <- gsub("M_All_PolyIC_","M_PolyIC_", colnames(CBmaster))

colnames(CBmaster)<- gsub("_Saline_","_All_Saline_", colnames(CBmaster))
colnames(CBmaster)<- gsub("F_All_Saline_","F_Saline_", colnames(CBmaster))
colnames(CBmaster) <- gsub("M_All_Saline_","M_Saline_", colnames(CBmaster))

colnames(CBmaster)<- gsub("_PolyICvsSaline_","_All_PolyICvsSaline_", colnames(CBmaster))
colnames(CBmaster)<- gsub("F_All_PolyICvsSaline_","F_PolyICvsSaline_", colnames(CBmaster))
colnames(CBmaster) <- gsub("M_All_PolyICvsSaline_","M_PolyICvsSaline_", colnames(CBmaster))

colnames(CBmaster) <- gsub("_All_All_","_All_", colnames(CBmaster))


#identify significance for each DEG column
#pull all significance columns
SigCol <- CBmaster2[ , grepl( "significance" , names( CBmaster2 ) ) ]
rownames(SigCol) <- CBmaster2$ensembl_gene_id

#detect DE at beginning of string
#sum(str_detect(CBmaster$significance_F_PolyIC_salvsTreg, '^DE')) > 0

#apply over columns. returns vector with TRUE if contains a significant DEG
Siglist <- sapply(SigCol, function(x) sum(str_detect(x, '^DE')) > 0 ) 

#filter for only true value containing lists of DEGs
SiglistTrue <- names(Siglist[which(Siglist)])

SiglistTrue <- gsub("significance_","",SiglistTrue)

#build master list of ensembl IDs for each DEG list:

masterDEGlist <- list()

for (i in 1:length(SiglistTrue)){
 
  df.subset <- CBmaster[ , grepl(SiglistTrue[i], names( CBmaster ) ) ]

  df.subset$ensembl_gene_id <- CBmaster$ensembl_gene_id
  df.subset$mgi_symbol <- CBmaster$mgi_symbol
  df.subset$description <- CBmaster$description
 
  #change column names
  colnames(df.subset)[6] <- substr(colnames(df.subset[6]), 1, 3) 
  colnames(df.subset)[5] <- substr(colnames(df.subset[5]), 1, 3)
  
  df.sig.up <- df.subset %>% filter(sig == "DE") %>%
    filter(dir == "Increase") %>% select(ensembl_gene_id)
  
  names(df.sig.up) <- paste(SiglistTrue[i],"Increase",sep="_")
  
  df.sig.down <- df.subset %>% filter(sig == "DE") %>%
    filter(dir == "Decrease") %>% select(ensembl_gene_id)
  
  names(df.sig.down) <- paste(SiglistTrue[i],"Decrease",sep="_")
  
  #make to list
  list1 <- append(df.sig.up,df.sig.down)
  
  #save to a master list
  masterDEGlist <- append(masterDEGlist, list1)
}

names(masterDEGlist)
head(masterDEGlist)

CBmasterlist <- masterDEGlist

#############################################################################################################################load DEGs from Hippocampus:
##########################################################################################################################

path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/Hippocampus/sexlitter"
setwd(path)

#load HC example data to get gene info:
HC <- read.csv("M_PolyIC_salvsTreg _DEgenes.csv")
HC_genes <- as.data.frame(HC[,3:10])
#HC_genes_data <- as.data.frame(HC[,c(11,13,14,15,16,17)])
#colnames(HC_genes_data) <- paste(colnames(HC_genes_data), unique(HC$comparison), sep="_" )


out.file <- data.frame(HC_genes$ensembl_gene_id)
out.file <- out.file[complete.cases(out.file),]

out.file <- data.frame(HC_genes$ensembl_gene_id)
names(out.file) <- c("ensembl_gene_id")


file.names <- dir(path, pattern ='*_DEgenes.csv')
file.names


for(i in 1:length(file.names)){
  file <- read.csv(file.names[i])
  
  #exract stats and rename columns to include comparison performed
  print(unique(file$comparison))
  
  HC_genes_data <- as.data.frame(file[,c(11,13,14,15,16,17)])
  colnames(HC_genes_data) <- paste(colnames(HC_genes_data), unique(file$comparison), sep="_" )
  HC_genes_data$ensembl_gene_id <- file$ensembl_gene_id
  HC_genes_data <- HC_genes_data[complete.cases(HC_genes_data),]
  
  #merge into datafram by ensembl id
  out.file <- merge(out.file, HC_genes_data, by="ensembl_gene_id")
  
  print("done")
  
}

#add in log2CPM per animal
HC_CPM <- read.csv("IndividualSamplesLog2CPM_GeneExpression.csv")
HC_CPM$X.1 <- NULL

HCmaster <- merge(out.file, HC_CPM, by.x="ensembl_gene_id",by.y="GeneID")

#write to file
write.csv(HCmaster, file = "HC_masterDEGs_log2CPM.csv")

#DEGs

#fix column names to contain M, F or All for comparisons

colnames(HCmaster)<- gsub("_PolyIC_","_All_PolyIC_", colnames(HCmaster))
colnames(HCmaster)<- gsub("F_All_PolyIC_","F_PolyIC_", colnames(HCmaster))
colnames(HCmaster) <- gsub("M_All_PolyIC_","M_PolyIC_", colnames(HCmaster))

colnames(HCmaster)<- gsub("_Saline_","_All_Saline_", colnames(HCmaster))
colnames(HCmaster)<- gsub("F_All_Saline_","F_Saline_", colnames(HCmaster))
colnames(HCmaster) <- gsub("M_All_Saline_","M_Saline_", colnames(HCmaster))

colnames(HCmaster)<- gsub("_PolyICvsSaline_","_All_PolyICvsSaline_", colnames(HCmaster))
colnames(HCmaster)<- gsub("F_All_PolyICvsSaline_","F_PolyICvsSaline_", colnames(HCmaster))
colnames(HCmaster) <- gsub("M_All_PolyICvsSaline_","M_PolyICvsSaline_", colnames(HCmaster))

colnames(HCmaster) <- gsub("_All_All_","_All_", colnames(HCmaster))
colnames(HCmaster)

HCmaster2 <- HCmaster

#identify significance for each DEG column
#pull all significance columns
SigCol <- HCmaster2[ , grepl( "significance" , names( HCmaster2 ) ) ]
rownames(SigCol) <- HCmaster2$ensembl_gene_id

#detect DE at beginning of string
#sum(str_detect(HCmaster$significance_F_PolyIC_salvsTreg, '^DE')) > 0

#apply over columns. returns vector with TRUE if contains a significant DEG
Siglist <- sapply(SigCol, function(x) sum(str_detect(x, '^DE')) > 0 ) 

#filter for only true value containing lists of DEGs
SiglistTrue <- names(Siglist[which(Siglist)])

SiglistTrue <- gsub("significance_","",SiglistTrue)
SiglistTrue

#build master list of ensembl IDs for each DEG list:

masterDEGlist <- list()

for (i in 1:length(SiglistTrue)){
  
  df.subset <- HCmaster[ , grepl(SiglistTrue[i], names( HCmaster ) ) ]
  
  df.subset$ensembl_gene_id <- HCmaster$ensembl_gene_id
  df.subset$mgi_symbol <- HCmaster$mgi_symbol
  df.subset$description <- HCmaster$description
  
  #change column names
  colnames(df.subset)[6] <- substr(colnames(df.subset[6]), 1, 3) 
  colnames(df.subset)[5] <- substr(colnames(df.subset[5]), 1, 3)
  
  df.sig.up <- df.subset %>% filter(sig == "DE") %>%
    filter(dir == "Increase") %>% select(ensembl_gene_id)
  
  names(df.sig.up) <- paste(SiglistTrue[i],"Increase",sep="_")
  
  df.sig.down <- df.subset %>% filter(sig == "DE") %>%
    filter(dir == "Decrease") %>% select(ensembl_gene_id)
  
  names(df.sig.down) <- paste(SiglistTrue[i],"Decrease",sep="_")
  
  #make to list
  list1 <- append(df.sig.up,df.sig.down)
  
  #save to a master list
  masterDEGlist <- append(masterDEGlist, list1)
}

names(masterDEGlist)
head(masterDEGlist)

HCmasterlist <- masterDEGlist

#############################################################################################################################load DEGs from FrontalCortex:/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/FrontalCortex/sexlitter_removePA92FC30
##########################################################################################################################

path ="/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/FrontalCortex/sexlitter_removePA92FC30"
setwd(path)

#load FC example data to get gene info:
FC <- read.csv("M_PolyIC_salvsTreg _DEgenes.csv")
FC_genes <- as.data.frame(FC[,3:10])
#FC_genes_data <- as.data.frame(FC[,c(11,13,14,15,16,17)])
#colnames(FC_genes_data) <- paste(colnames(FC_genes_data), unique(FC$comparison), sep="_" )


out.file <- data.frame(FC_genes$ensembl_gene_id)
out.file <- out.file[complete.cases(out.file),]

out.file <- data.frame(FC_genes$ensembl_gene_id)
names(out.file) <- c("ensembl_gene_id")


file.names <- dir(path, pattern ='*_DEgenes.csv')
file.names


for(i in 1:length(file.names)){
  file <- read.csv(file.names[i])
  
  #exract stats and rename columns to include comparison performed
  print(unique(file$comparison))
  
  FC_genes_data <- as.data.frame(file[,c(11,13,14,15,16,17)])
  colnames(FC_genes_data) <- paste(colnames(FC_genes_data), unique(file$comparison), sep="_" )
  FC_genes_data$ensembl_gene_id <- file$ensembl_gene_id
  FC_genes_data <- FC_genes_data[complete.cases(FC_genes_data),]
  
  #merge into datafram by ensembl id
  out.file <- merge(out.file, FC_genes_data, by="ensembl_gene_id")
  
  print("done")
  
}

#add in log2CPM per animal
FC_CPM <- read.csv("IndividualSamplesLog2CPM_GeneExpression.csv")
FC_CPM$X.1 <- NULL

FCmaster <- merge(out.file, FC_CPM, by.x="ensembl_gene_id",by.y="GeneID")

#write to file
write.csv(FCmaster, file = "FC_masterDEGs_log2CPM.csv")

#DEGs

#fix column names to contain M, F or All for comparisons

colnames(FCmaster)<- gsub("_PolyIC_","_All_PolyIC_", colnames(FCmaster))
colnames(FCmaster)<- gsub("F_All_PolyIC_","F_PolyIC_", colnames(FCmaster))
colnames(FCmaster) <- gsub("M_All_PolyIC_","M_PolyIC_", colnames(FCmaster))

colnames(FCmaster)<- gsub("_Saline_","_All_Saline_", colnames(FCmaster))
colnames(FCmaster)<- gsub("F_All_Saline_","F_Saline_", colnames(FCmaster))
colnames(FCmaster) <- gsub("M_All_Saline_","M_Saline_", colnames(FCmaster))

colnames(FCmaster)<- gsub("_PolyICvsSaline_","_All_PolyICvsSaline_", colnames(FCmaster))
colnames(FCmaster)<- gsub("F_All_PolyICvsSaline_","F_PolyICvsSaline_", colnames(FCmaster))
colnames(FCmaster) <- gsub("M_All_PolyICvsSaline_","M_PolyICvsSaline_", colnames(FCmaster))

colnames(FCmaster) <- gsub("_All_All_","_All_", colnames(FCmaster))
colnames(FCmaster)

FCmaster2 <- FCmaster

#identify significance for each DEG column
#pull all significance columns
SigCol <- FCmaster2[ , grepl( "significance" , names( FCmaster2 ) ) ]
rownames(SigCol) <- FCmaster2$ensembl_gene_id

#detect DE at beginning of string
#sum(str_detect(FCmaster$significance_F_PolyIC_salvsTreg, '^DE')) > 0

#apply over columns. returns vector with TRUE if contains a significant DEG
Siglist <- sapply(SigCol, function(x) sum(str_detect(x, '^DE')) > 0 ) 

#filter for only true value containing lists of DEGs
SiglistTrue <- names(Siglist[which(Siglist)])

SiglistTrue <- gsub("significance_","",SiglistTrue)
SiglistTrue

#build master list of ensembl IDs for each DEG list:

masterDEGlist <- list()

for (i in 1:length(SiglistTrue)){
  
  df.subset <- FCmaster[ , grepl(SiglistTrue[i], names( FCmaster ) ) ]
  
  df.subset$ensembl_gene_id <- FCmaster$ensembl_gene_id
  df.subset$mgi_symbol <- FCmaster$mgi_symbol
  df.subset$description <- FCmaster$description
  
  #change column names
  colnames(df.subset)[6] <- substr(colnames(df.subset[6]), 1, 3) 
  colnames(df.subset)[5] <- substr(colnames(df.subset[5]), 1, 3)
  
  df.sig.up <- df.subset %>% filter(sig == "DE") %>%
    filter(dir == "Increase") %>% select(ensembl_gene_id)
  
  names(df.sig.up) <- paste(SiglistTrue[i],"Increase",sep="_")
  
  df.sig.down <- df.subset %>% filter(sig == "DE") %>%
    filter(dir == "Decrease") %>% select(ensembl_gene_id)
  
  names(df.sig.down) <- paste(SiglistTrue[i],"Decrease",sep="_")
  
  #make to list
  list1 <- append(df.sig.up,df.sig.down)
  
  #save to a master list
  masterDEGlist <- append(masterDEGlist, list1)
}

names(masterDEGlist)
head(masterDEGlist)

FCmasterlist <- masterDEGlist

#############################################################################################################################combined DEG liss across tissues
##########################################################################################################################
names(FCmasterlist) <- paste("FC_",names(FCmasterlist), sep="")
names(HCmasterlist) <- paste("HC_",names(HCmasterlist), sep="")
names(CBmasterlist) <- paste("CB_",names(CBmasterlist), sep="")

BrainList <- c(FCmasterlist, HCmasterlist, CBmasterlist)

#sex specific lists
SexBrainList <- BrainList[!grepl("_All_",names(BrainList))]


#universe of all possible genes: genes detected in at least 1 of teh 3 tissues
CBuniverse <- CBmaster$ensembl_gene_id
FCuniverse <- FCmaster$ensembl_gene_id
HCuniverse <- HCmaster$ensembl_gene_id

Universe <- c(CBuniverse,FCuniverse,HCuniverse)

#unique IDs : 19669
Universe <- unique(Universe)

#save
save(BrainList,SexBrainList, Universe,file="BrainDEGlists.Rdata")
#############################################################################################################################Upset plot
##########################################################################################################################
setwd("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/DEGcomparisons")
library(UpSetR)

pdf(file = "UpSetPlot_AllDEGs.pdf", wi = 14, he = 10, useDingbats=F)

upset(fromList(BrainList),
      order.by = "freq",
      #number.angles =10,
     # keep.order = TRUE,
      text.scale = c(2),
     nsets = length(names(BrainList)))

dev.off()

#repeat with Sex specific lists only

names(SexBrainList)

pdf(file = "UpSetPlot_SexonlyDEGs.pdf", wi = 14, he = 10, useDingbats=F)

upset(fromList(SexBrainList),
      order.by = "freq",
      #number.angles =10,
      # keep.order = TRUE,
      text.scale = c(2),
      nsets = length(names(SexBrainList)))

dev.off()


#######################################################################
# GO enrichment of each list using ClusterCompare
#######################################################################
#library(devtools)
#devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
library("clusterProfiler")
#The input for geneCluster parameter should be a named list of gene IDs. 


l <- rapply(SexBrainList, length, how="list")
l <- as.data.frame(l)
write.csv(l, file="DEGcounts.csv")

x2 <- SexBrainList
#######################################################################
#for KEGG convert to entreze
#######################################################################

#BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
library(xlsx)

listDEGs2  <- NULL

for (i in 1:length(names(x2))){
  
  e <- bitr(x2[[i]],fromType='ENSEMBL',toType="ENTREZID", OrgDb="org.Mm.eg.db")
  head(e$ENTREZID)
  
  #save to output list:
  listDEGs2[[i]] <- unique(e$ENTREZID)
  
}

names(listDEGs2) <- names(x2)

#universe = all detected genes in experiment
universe_enterz <- bitr(Universe,fromType='ENSEMBL',toType="ENTREZID", OrgDb="org.Mm.eg.db")
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
                     universe      = unique(Universe),
                     OrgDb         = org.Mm.eg.db,
                     ont           = "CC",
                     keyType       = 'ENSEMBL',
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
                      universe      = unique(Universe),
                      OrgDb         = org.Mm.eg.db,
                      ont           = "BP",
                      keyType       = 'ENSEMBL',
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
                     universe      = unique(Universe),
                     OrgDb         = org.Mm.eg.db,
                     ont           = "mf",
                     keyType       = 'ENSEMBL',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     fun = "enrichGO")
head(mf)
dim(mf)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
mf2 <- simplify(mf, cutoff=0.7, by="p.adjust", select_fun=min)
dim(mf2)


pdf('GO_MFEnrichment_allDEGs.pdf',w=12,h=14)
d <- dotplot(mf2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

mf2 <- as.data.frame(mf2)
write.csv(mf2,file="GO_MFenrichments_allDEGs.csv")


#######################################################################
# MGEnrichment database Fisher's Exact Tests
#######################################################################


setwd("/Users/aciernia/Sync/collaborations/Ashwood/MIATcellRNAseq/listenrichments")

load("Mouse_Human_GenelistDatabaseJuly2023.RData")

#make to genelist
msdb <- mouse.master3 %>% dplyr::select(ensembl_gene_id,listname)
DatabaseList <- split(msdb$ensembl_gene_id, f = msdb$listname )

##########################################################################################################
#Overlap Function: Fisher's one tailed exact test
#takes a targetlistname > list of ensembl gene ids of interest
#genelists > list of ensembl ids for each gene list in the database
#genomesize > single numeric value that is a count of all unique genes either all mm10 genes or all genes in the database
##########################################################################################################
#BiocManager::install("GeneOverlap")

library("GeneOverlap")

Overlap_fxn <- function(targetlistname,genelists,genomesize){
  out <- NULL
  target <- targetlistname[] #list 
  inputlistname <- names(targetlistname)
  
  for (i in 1:length(genelists)) { 
    
    #call gene overlaps
    go.obj <- newGeneOverlap(target[[1]],
                             genelists[[i]],
                             genome.size=genomesize)
    
    #perform test
    go.obj <- testGeneOverlap(go.obj) #returns onetailed pvalue
    #return odds ratio:
    OR <- getOddsRatio(go.obj)
    pvalue <- getPval(go.obj)
    
    #extract contingency table
    CT <- getContbl(go.obj)
    notAnotB <- CT[1,1]
    inAnotB <- CT[1,2]
    inBnotA <- CT[2,1]
    inBinA <- CT[2,2]
    
    CTlist <- cbind(notAnotB,inAnotB,inBnotA,inBinA)
    
    
    #two sided fisher's exact test
    #test <- fisher.test(CT,alternative='two.sided')
    
    #get gene list B:
    intersection <- go.obj@intersection
    intersection_ensembl <- paste(as.character(intersection),collapse=", ",sep="")
    
    #get intersection gene names
    intersection_genenames <- mouse_genes$mgi_symbol[which(mouse_genes$ensembl_gene_id %in% intersection)]
    intersection_genenames <- paste(as.character(intersection_genenames),collapse=", ",sep="")
    
    #get listname
    listname <- paste(names(genelists[i]))
    
    results <- cbind(listname,pvalue,OR, CTlist,intersection_ensembl,intersection_genenames )
    
    names(results) <- c("listname","pvalue","OR","notAnotB","inAnotB","inBnotA","inBinA","ensembl","geneID")
    out <- rbind(out,results) 
    
  }
  
  #remove first row as overlap with self:
  out <- as.data.frame(out)
  # out2 <- out[- grep(inputlistname, out$listname),]
  
  #add in targetlist name (assumes first list is the input)
  out$targetlist <- paste(inputlistname)
  
  
  rownames(out) <- NULL
  
  #return results
  return(out)
}



##########################################################################################################
#make lists for each brain region
##########################################################################################################


#sex specific lists
SexBrainList <- BrainList[!grepl("_All_",names(BrainList))]
CBBrainlist <- SexBrainList[grepl("CB_",names(SexBrainList))]
FCBrainlist <- SexBrainList[grepl("FC_",names(SexBrainList))]
HCBrainlist <- SexBrainList[grepl("HC_",names(SexBrainList))]

##########################################################################################################
#CB enrichments
##########################################################################################################

genomesize <- length(Universe)
targetsList <- CBBrainlist

dfoverlaps <- NULL
for (i in 1:length(targetsList)) { 
  
  tmp <- Overlap_fxn(targetsList[i], DatabaseList, genomesize)
  dfoverlaps <- rbind(tmp,dfoverlaps)
}

#remove list overlap with self if targetlist is from the database
dfoverlaps2 <- subset(dfoverlaps, listname != targetlist)

#adjust pvalue
dfoverlaps2$pvalue <- as.numeric(as.character(dfoverlaps2$pvalue))
dfoverlaps2$FDR <- p.adjust(dfoverlaps2$pvalue, method='fdr')

#add back in list info
masterlist <- mouse.master3 %>% dplyr::select(-ensembl_gene_id,-mgi_symbol,-entrezgene_id) %>% distinct()

m <- merge(dfoverlaps2,masterlist, by=c("listname"), all.x=T)
m <- unique(m)

write.csv(m,"CBDEGs_enrichments_genelists_July2023.csv")

##########################################################################################################
#FC enrichments
##########################################################################################################

genomesize <- length(Universe)
targetsList <- FCBrainlist

dfoverlaps <- NULL
for (i in 1:length(targetsList)) { 
  
  tmp <- Overlap_fxn(targetsList[i], DatabaseList, genomesize)
  dfoverlaps <- rbind(tmp,dfoverlaps)
}

#remove list overlap with self if targetlist is from the database
dfoverlaps2 <- subset(dfoverlaps, listname != targetlist)

#adjust pvalue
dfoverlaps2$pvalue <- as.numeric(as.character(dfoverlaps2$pvalue))
dfoverlaps2$FDR <- p.adjust(dfoverlaps2$pvalue, method='fdr')

#add back in list info
masterlist <- mouse.master3 %>% dplyr::select(-ensembl_gene_id,-mgi_symbol,-entrezgene_id) %>% distinct()

m <- merge(dfoverlaps2,masterlist, by=c("listname"), all.x=T)
m <- unique(m)

write.csv(m,"FCDEGs_enrichments_genelists_July2023.csv")

##########################################################################################################
#HC enrichments
##########################################################################################################

genomesize <- length(Universe)
targetsList <- HCBrainlist

dfoverlaps <- NULL
for (i in 1:length(targetsList)) { 
  
  tmp <- Overlap_fxn(targetsList[i], DatabaseList, genomesize)
  dfoverlaps <- rbind(tmp,dfoverlaps)
}

#remove list overlap with self if targetlist is from the database
dfoverlaps2 <- subset(dfoverlaps, listname != targetlist)

#adjust pvalue
dfoverlaps2$pvalue <- as.numeric(as.character(dfoverlaps2$pvalue))
dfoverlaps2$FDR <- p.adjust(dfoverlaps2$pvalue, method='fdr')

#add back in list info
masterlist <- mouse.master3 %>% dplyr::select(-ensembl_gene_id,-mgi_symbol,-entrezgene_id) %>% distinct()

m <- merge(dfoverlaps2,masterlist, by=c("listname"), all.x=T)
m <- unique(m)

write.csv(m,"HCDEGs_enrichments_genelists_July2023.csv")

