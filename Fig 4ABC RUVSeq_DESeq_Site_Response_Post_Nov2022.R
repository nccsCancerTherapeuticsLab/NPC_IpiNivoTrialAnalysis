####Differential Expression Analysis####
#Method is adapted from https://github.com/bhattacharya-a-bt/CBCS_normalization

require(NanoStringQCPro)
require(ggplot2)
require(EnvStats)
require(readxl)
require(dplyr)
require(tidyverse)
require(RUVSeq)
require(DESeq2)
require(limma)
require(matrixStats)
require(MASS)

#Functions were adapted from https://github.com/bhattacharya-a-bt/CBCS_normalization/blob/master/nanostring_RUV_functions.R
RUV_total <- function(raw,pData,fData,k,hkgenes = NULL,exclude = NULL){
  
  ### INPUT: raw - p x n raw expressions with p genes and n samples
  ###        pData - phenotype metadata across samples
  ###        fData - feature metadata across genes
  ###        k - number of dimensions of unwanted variation estimated
  ###        exclude - vector of gene names to exclude
  
  library(RUVSeq)
  library(DESeq2)
  library(limma)
  library(matrixStats)
  
  if (!is.null(hkgenes)){
    
    fData(set)$Class[rownames(set) %in% hkgenes] = 'Housekeeping'
    
  }
  
  fData = fData[rownames(raw),]
  int = intersect(rownames(raw),rownames(fData))
  fData = fData[int,]
  raw = raw[int,]
  
  set <- newSeqExpressionSet(as.matrix(round(raw)),
                             phenoData=pData,
                             featureData=fData)
  
  cIdx <- rownames(set)[fData(set)$Class == "Housekeeping"]
  cIdx = cIdx[!(cIdx %in% exclude)]
  #x <- as.factor(pData$Group)
  set <- betweenLaneNormalization(set, which="upper")
  set <- RUVg(set, cIdx, k=k)
  dds <- DESeqDataSetFromMatrix(counts(set),colData=pData(set),design=~1)
  rowData(dds) <- fData
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)$dispGeneEst <- disp
  dds <- estimateDispersionsFit(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat <- removeBatchEffect(mat, covariates=covars)
  assay(vsd) <- mat
  return(list(set = set,vsd = vsd))
  
  
  
}

imagingQC <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for imaging quality
  
  fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
  if (!(fovRatio > .75)) {return('Flag')}
  if (fovRatio > .75) {return('No flag')}
  
}

bindingDensityQC <- function(rcc,low,high){
  
  
  #### INPUT: rcc - input from rcc
  ####         low, high - the lower and upper limits for binding density
  #### OUTPUT: flag for binding density
  
  bd = as.numeric(rcc$Lane_Attributes[6])
  if(!(bd < high & bd > low)) {return('Flag')}
  if (bd < high & bd > low) {return('No flag')}
  
  
}

limitOfDetectionQC <- function(rcc,numSD = 0){
  
  #### INPUT: rcc - input from rcc
  ####         numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: flag for limit of detection
  
  counts = rcc$Code_Summary
  posE = as.numeric(counts$Count[counts$Name == 'POS_E'])
  negControls = as.numeric(counts$Count[grepl('NEG',counts$Name)])
  if(!(posE > mean(negControls) + numSD*sd(negControls))) {return('Flag')}
  if (posE > mean(negControls) + numSD*sd(negControls)) {return('No flag')}
  
}

positiveLinQC <- function(rcc){
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for linearity for positive controls
  
  
  counts = rcc$Code_Summary
  posControls = as.numeric(counts$Count[grepl('POS_',counts$Name)])
  known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
  r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
  if(!(r2 > .95) | is.na(r2)) {return('Flag')}
  if(r2 > .95) {return('No flag')}
  
}

makeRLEplot <- function(data,metadata,id){
  
  #### INPUT: data - matrix of expressions with genes on rows and samples on columns
  ####        metadata - matrix of metadata with a column that corresponds to the colnames of data
  ####        id - colname of sample ids
  #### OUTPUT: ggplot2 RLE plot
  
  data = data - apply(data,1,median)
  stack = stack(data)
  colnames(stack)[1] = id
  stackPlot = merge(stack,metadata,by=id)
  colnames(stackPlot)[1:2] = c('Sample','values')
  rle_plots = ggplot(data = stackPlot,aes(x = Sample,y = values, color = ER_status)) +
    geom_boxplot(coef = 6) + theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=20),
          legend.text=element_text(size=20),
          strip.text = element_text(size=24),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey", 
                                      fill = NA, size = .1),
          legend.position = 'bottom',
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + xlab('Sample') +
    ylab('Median deviation of log expression') + ylim(c(-4,4))
  return(rle_plots)
  
}

#Differential Expression analysis is done on all samples and by removing one sample at time to find the most consistent differential genes

n=0 #Loop count
#Import Annotation file
nanostring_annot <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/Nanotring Annotation All.xlsx")
nanostring_annot <- nanostring_annot[,-1]
colnames(nanostring_annot)[1] <- "SampleID"
nsd <- nanostring_annot$SampleID[nanostring_annot$Site_Response!="SD"]

#Read how many files are there in total
nfiles.RCC = list.files(path = "C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC RCC/")
nfiles.RCC = nfiles.RCC[grepl('RCC',nfiles.RCC)] #Choose only RCC files
nfiles.RCC = nfiles.RCC[grepl(paste(nsd, collapse="|"),nfiles.RCC)] #Choose only PDPR
nfiles.RCC = nfiles.RCC[grepl('Post',nfiles.RCC)] #Choose only Post-treatment File
nfiles.RCC

#Loop to eliminate 1 sample at a time
for (n in n:length(nfiles.RCC)){
setwd("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC RCC/")
files.RCC = list.files(path = "C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC RCC/")
files.RCC = files.RCC[grepl('RCC',files.RCC)] #Choose only RCC files
files.RCC = files.RCC[grepl(paste(nsd, collapse="|"),files.RCC)] #Choose only PDPR
files.RCC = files.RCC[grepl('Post',files.RCC)] #Choose only Pre-treatment File

#Remove sample 
if (n==0) { nm <- "all"} else { 
  nm <- sub("P.*", "", files.RCC[n]) #extract removed sample name
  files.RCC = files.RCC[-c(n)]}
files.RCC

#Creating Raw Expression Matrix
raw_expression = as.data.frame(matrix(nrow = as.integer("784"),ncol = length(files.RCC)+2))
colnames(raw_expression)[1:2] = c('Gene','Class')
pData = as.data.frame(matrix(nrow = length(files.RCC),ncol = 11))
colnames(pData) = c('BCAC_ID','SampleID','Owner','Comments','Date','GeneRLF','SystemAPF','imagingQC',
                    'bindingDensityQC','limitOfDetectionQC','positiveLinearityQC')
raw_expression[,1:2] = readRcc(files.RCC[1])$Code_Summary[,c(2,1)]


for (i in 1:length(files.RCC)){
  
  print(i)
  rcc = readRcc(files.RCC[i])
  raw = rcc$Code_Summary
  
  raw_expression[,i+2] = as.numeric(raw$Count)
  colnames(raw_expression)[i+2] = strsplit(files.RCC[i],'_')[[1]][1]
  pData[i,2:7] = as.vector(rcc$Sample_Attributes)
  pData$imagingQC[i] = imagingQC(rcc)
  pData$bindingDensityQC[i] = bindingDensityQC(rcc,.05,2.25)
  pData$limitOfDetectionQC[i] = limitOfDetectionQC(rcc)
  pData$positiveLinearityQC[i] = positiveLinQC(rcc)
}

setwd("C:/Users/chris/Dropbox/R/New Clean Environment/NanostringDiff")

raw_expression <- raw_expression%>% rename_at(vars(contains(".RCC")), funs(str_replace(.,".RCC", ""))) #Remove .RCC from colnames
s_id <- colnames(raw_expression)[-c(1,2)]
pData$SampleID <- s_id

#Add annotation to pData, remove unwanted columns
pData <- left_join(pData, nanostring_annot)
pData <- pData[,-c(1,3:11)] #Remove redundant column such as QC flags
pData$Treatment <- factor(pData$Treatment, levels = c("Pre","On"))
pData$Site_Response <- factor(pData$Site_Response, levels = c("PD","PR","SD"))
pData$Sample <- factor(pData$Sample)

#remove POS & NEG control from raw data
raw_expression <- raw_expression %>%filter(Class %in% c("Endogenous","Housekeeping") )

raw = raw_expression[,-c(1:2)]
fData = raw_expression[,c(1:2)]
rownames(raw) = fData$Gene
cIdx <- fData$Gene[fData$Class == "Housekeeping"]
pData$HK_Gene_Miss = colSums(raw[cIdx,] == 0)
rownames(fData) = fData$Gene
rownames(raw) = fData$Gene
rownames(pData) = colnames(raw)

#### CHECK IF Housekeeping genes are associated with Response
hk_raw = raw[cIdx,]
pval = vector(length = nrow(hk_raw))

for (i in 1:nrow(hk_raw)){
  reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(pData$Site_Response)) 
  pval[i] = coef(summary(reg))[2,4]
}
#making sure factor HK genes does not differ between different response
sum(pval <= .05) 

k = 1
vsd = RUV_total(raw,pData,fData,k = k)$vsd
set = RUV_total(raw,pData,fData,k = k)$set
save(vsd,file = paste0("deangelo_ruv_vsd_k_",k,".rda"))
save(set,file = paste0("deangelo_ruv_set_k_",k,".rda"))

plotRLE(set@assayData[["counts"]])
plotRLE(set@assayData[["normalizedCounts"]])

i = 1

results = paste0("deangelo_deg_k",i,".csv")
load(paste0("deangelo_ruv_set_k_",i,".rda"))

#DE Comparision
dds <- DESeqDataSetFromMatrix(countData = counts(set)[1:750,],
                              colData = pData(set),
                              design = ~ Site_Response)
dds <- DESeq(dds)
resultsNames(dds)

res <- as.data.frame(results(dds, name="Site_Response_PR_vs_PD"))

#save as table
write.table(res, paste0("result/Response/Post_Site/DE_result", "_no", nm ,".csv"), sep =",", row.names = TRUE, col.names = NA)
}


#Merge Different DE tables to find overlapping genes
library(dplyr)
library(reshape)
library(tidyverse)
library(plyr)

directory<- "C:/Users/chris/Dropbox/R/New Clean Environment/NanostringDiff/result/Response/Post_Site"
file_list <- list.files(directory)

load_file <- function (filename) {
  setwd(directory)
  db <- read.csv(filename) 
  colnames(db)[1] <- "Gene"
  db <- db %>% dplyr::select(Gene, log2FoldChange, pvalue) %>% filter(pvalue <= 0.05)
  file_label <- sub("DE_result_", "", filename)
  file_label <- sub(".csv", "", file_label)
  db$label <- file_label
  return(db)
}

allfile <- lapply(file_list, load_file) #Apply function to all the files in the list
file_list_name <- sub("DE_result_", "", file_list) #Remove Prefix from file name
file_list_name <- sub(".csv", "", file_list_name) #Remove suffix from file name
names(allfile) <- file_list_name
merged<- merge_all(allfile) #combine all list into a single file
merged_summary <- merged %>% group_by(Gene) %>% count("Gene") #count the frequency of the genes in each table

setwd("C:/Users/chris/Dropbox/R/New Clean Environment/NanostringDiff/result")

#Compiling all the files together
df_final <- allfile %>% purrr::reduce(full_join, by = "Gene") %>% left_join(merged_summary) %>% arrange(desc(freq))
write.table(df_final, 'full_list_Post_Site_PDPR_with_fc_pvalue.csv', sep = ',', row.names = FALSE)

#Filtering for genes that appeared in all DE result
df_filtered <- df_final %>% filter(freq == max(freq)) %>% arrange(log2FoldChange.x)
write.table(df_filtered, 'final_list_Post_Site_PDPR_with_fc_pvalue.csv', sep = ',', row.names = FALSE)

#Plot the Number of Genes & Freq
db_final_summary <- df_final %>% dplyr::select(Gene, freq, contains("label")) 
db_final_summary$freq <- factor(db_final_summary$freq)
count_freq<- db_final_summary %>% dplyr::select(freq) %>% group_by(freq) %>% dplyr::summarise(Gene = n())
ggplot(count_freq, aes(x=freq, y=Gene)) + 
  geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=0.4) +theme_classic()

#Plot Number of Genes and Database
merged$label <- factor(merged$label)
count_freq1<- merged %>% dplyr::select(label) %>% group_by(label) %>% dplyr::summarise(Count = n())
ggplot(count_freq1, aes(x=label, y=Count)) + 
  geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=0.4) +theme_classic()

####Plot Volcano Plot####
library(ggrepel)
library(ggplot2)
DE_ontreat_volcano <- read.csv("C:/Users/chris/Dropbox/R/New Clean Environment/NanostringDiff/result/Response/Post_Site/DE_result_noall.csv", header = TRUE, sep = ',')
DE_ontreat_volcano$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DE_ontreat_volcano$diffexpressed[DE_ontreat_volcano$log2FoldChange > 0 & DE_ontreat_volcano$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DE_ontreat_volcano$diffexpressed[DE_ontreat_volcano$log2FoldChange < 0 & DE_ontreat_volcano$pvalue < 0.05] <- "DOWN"

v <- ggplot(data=DE_ontreat_volcano, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=X)) +
  geom_point(size = 3) + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "none", aspect.ratio = 1) +
  geom_text_repel(size = 5) +
  scale_color_manual(values=c("#1e90ff", "black", "#f768a1"))
v

while (!is.null(dev.list()))  dev.off() #clear device
pdf("post_Site_PRPD_volcano.pdf", width = 8, height = 7)
plot(v)
dev.off()

####Plot Heatmap####
library("pheatmap")
normalised_data <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/NPC All Normalised Data.xlsx", 
                                      sheet = "Log2")
annotation_file <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/Nanotring Annotation All.xlsx")
hm_annot <- annotation_file %>% dplyr::select(SampleID, Site_Response, EBV, Smoking_Status, Sex, Treatment) %>% as.data.frame() %>% tibble::column_to_rownames(var="SampleID")
hm_annot_PDPR <- hm_annot %>% filter(Site_Response != "SD" & Treatment == "On")
hm_annot_PDPR <- within(hm_annot_PDPR, rm(Treatment))

#Extract Genes that has max frequency (appear in all the iteration)
deg <- df_final$Gene[df_final$freq == max(df_final$freq)]
norm_mat <- as.matrix(normalised_data[,-1])
rownames(norm_mat) <- normalised_data$Gene

#Create Annotation color
ebv_col30k = c("#37a124","#e5c064")
names(ebv_col30k) = levels(as.factor(hm_annot_PDPR$EBV))
annot_color = list(
  Site_Response = c(PR= "#f768a1", PD="#1e90ff", SD ="#e6e600"), 
  EBV = ebv_col30k,
  Smoking_Status = c(Smoker= "#282828", `Ex-smoker` ="#767676", `Non-smoker` = "#f6f6f6"), Sex = c(Female ='plum', Male ='lightskyblue'))

#Subset matrix to select only DE genes (row) and PDPR patients (column)
hm_db <- norm_mat[deg,rownames(hm_annot_PDPR)]

hm <- pheatmap(hm_db,
         scale = "row", 
         show_rownames = FALSE, show_colnames = TRUE,
         border_color = NA,
         color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
         treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
         fontsize_row = 9,
         annotation_col = hm_annot_PDPR,
         annotation_colors = annot_color)

while (!is.null(dev.list()))  dev.off() #clear device
pdf("post_Site_PRPD_heatmap_unlabelled.pdf", width = 4.5, height = 6)
hm
dev.off()

####Pathway####
library(DOSE)
library(clusterProfiler)
library("AnnotationHub")
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)

res <- df_filtered[,c(1:2)]
colnames(res)[1] <- "SYMBOL"

gene.df <- bitr(res$SYMBOL, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
DE_result <- inner_join(gene.df, res) 
geneList = DE_result$log2FoldChange.x #create vector with log2FC
names(geneList) = as.character(DE_result$ENTREZID) #give names to each value in vector
geneList = sort(geneList, decreasing = TRUE) #order the genelist in decreasing order

#Plotting Upreg & downreg genes
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf$group <- "High in PR"
mydf$group[mydf$FC < 0] <- "High in PD"
formula_res_go <- compareCluster(Entrez~group, data=mydf, fun="enrichGO", OrgDb= org.Hs.eg.db, ont = "all", readable = TRUE)
#plot
g<- dotplot(formula_res_go, x="group",showCategory = 10, font.size = 11) #Change Show Category to change the number of pathway shown
g

write.table(data.frame(formula_res_go), "Site_Response_post_pathway.csv", sep =",", row.names = FALSE)

#Manually choose 10 significant pathways, removing those that falls under the same category
chosen <- c("GO:0009615",
  "GO:0019221",
  "GO:0009314",
  "GO:0044772",
  "GO:0045069",
  "GO:0097193",
  "GO:2001234",
  "GO:0032608",
  "GO:0034340",
  "GO:0140374",
  "GO:0072678",
  "GO:0031295",
  "GO:0032609",
  "GO:0050852",
  "GO:0050851",
  "GO:0042113",
  "GO:0046651",
  "GO:0030217",
  "GO:1903037",
  "GO:0030098",
  "GO:0042110")

pathway_table <- data.frame(formula_res_go)
pathway_table <- pathway_table %>%filter(ID %in% chosen == TRUE)

pathway_ordered <- pathway_table[order(sapply(pathway_table$ID, function(x) which(x == chosen))), ]

#Plot Pathways
adj_res <- formula_res_go
adj_res@compareClusterResult <- pathway_ordered
h<- dotplot(adj_res, x="group",showCategory = 10, font.size = 11)
h

while (!is.null(dev.list()))  dev.off() #clear device
pdf("Pathway_Post_Site_PDPR.pdf", width = 4, height = 6)
h
dev.off()
