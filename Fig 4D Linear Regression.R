library(tidyverse)
library(broom)
library(readxl)
library(gridExtra)

#Load Normalised data
data <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/NPC All Normalised Data.xlsx", 
                   sheet = "Log2")
data <- data %>% tibble::column_to_rownames("Gene")
data <- t(data) %>% as.data.frame() %>% tibble::rownames_to_column(var = "SampleID")
annotationfile <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/Nanotring Annotation All.xlsx")
annotationfile <- select(annotationfile, c("SampleID", "Treatment","Site_Response", "TTP")) %>% na.omit()

normalised_data_annot <- left_join(annotationfile, data)
normalised_data_annot$Site_Response <- factor(normalised_data_annot$Site_Response, levels = c("PR","SD","PD"))

#Filtering for Pre or On treatment
norm_pre <- dplyr::filter(normalised_data_annot, Treatment == "Pre") 
pre_ttp <- norm_pre[-c(1:3)]
norm_post <- dplyr::filter(normalised_data_annot, Treatment == "On")
post_ttp <- norm_post[-c(1:3)]

#Pearson Correlation for Table
pcor <- function (db){
  require(Hmisc)
  dbname = deparse(substitute(db))
  # cormat : matrix of the correlation coefficients
  # pmat : matrix of the correlation p-values
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  res2<-rcorr(as.matrix(db)) #calculate r & p value (Pearson)
  final_cor <- flattenCorrMatrix(res2$r, res2$P) #create summary table containing r and p value
  #Filter for TTP Only -> 
  final_cor <- filter(final_cor, row == "TTP") %>% filter(p < 0.05) %>% arrange(p)
  final_cor <- final_cor[,-1]
  colnames(final_cor) <- c("Gene", "r", "P-Value")
  write.table(final_cor, paste0("pearson_cor_",dbname, ".csv"), sep =',', row.names = FALSE)
  return(final_cor)
}

out_pre_ttp <- pcor(pre_ttp)
out_post_ttp <- pcor(post_ttp)

####Plot Linear Curve####
fil.col <- c('#f768a1','#e6e600', '#1e90ff')
library(ggrepel)
require(ggplot2)
ggplotRegression <- function (db,dep_val,gene) {
  eqn <- formula(paste0(dep_val," ~ ",gene))
  t <- cor(db[,dep_val], db[,gene], method = "pearson") #calculate r value
  t<- round(t, digits = 3) #3 decimal place
  fit <- lm(eqn, data = db)
  ggplot(db, aes(x = !!ensym(dep_val), y = !!ensym(gene))) + 
    geom_point(aes(fill = Site_Response), colour="black",pch=21, size = 5) +  guides(fill="legend") +
    stat_smooth(method = "lm", col = "red") + theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          aspect.ratio = 1,
          axis.title.x = element_blank(), #remove x-axis title
          axis.title.y = element_blank(), #remove y-axis title
          legend.position = "none") + #remove legend
    scale_fill_manual(values= fil.col)+ labs(title = gene, subtitle = paste(" r = ", t ,", P =",signif(summary(fit)$coef[2,4], 3))) +
    ylab("log2 normalised expression") # +
    #geom_label_repel(aes(label = SampleID))
}

ggplotRegression(norm_post, "TTP","CD4")

ggpost_ttp <- function(gene){
  ggplotRegression(norm_post, "TTP",gene)
}

while (!is.null(dev.list()))  dev.off() #prevent pdf write problem
pdf(paste0('result/correlationplot_post_ttp.pdf'))
lapply(out_post_ttp[[1]], ggpost_ttp)
dev.off()


#Plot chosen genes
ttp_post_chosen <- c("PDCD1", "HAVCR2", "CCL5", "CD244","MSH6", "LDHB", "BIRC5", "EGFR")

#Print them one by one
while (!is.null(dev.list()))  dev.off() #prevent pdf write problem
pdf('result/correlationplot_post_ttp_chosen.pdf')
lapply(ttp_post_chosen, ggpost_ttp)
dev.off()

library(gridExtra)
plist <- lapply(ttp_post_chosen, ggpost_ttp)
f<- grid.arrange(grobs = plist, ncol = 4) ## display plot
ggsave(file = "plot.pdf", arrangeGrob(grobs = plist, ncol = 4), width = 25, height = 15,units = "cm")  ## save plot

