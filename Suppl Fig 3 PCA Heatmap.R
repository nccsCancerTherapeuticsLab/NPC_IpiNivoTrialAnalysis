####PCA Plot####
library("factoextra")
library("ggfortify")
library("cluster")
library("readxl")
normalised_annot_log2 <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/normalised_annot_log2.xlsx")
row.names(normalised_annot_log2) <- unlist(normalised_annot_log2[,1])
normalised_annot_log2$Response<- factor(normalised_annot_log2$Response, levels = c("Partial_Response", "Stable_Disease", "Progressive_Disease"))
train <- normalised_annot_log2
train <- train[,-c(1:4)]
row.names(train) <-unlist(normalised_annot_log2[,1])
response_group <- factor(normalised_annot_log2$Response, levels = c("Partial_Response", "Stable_Disease", "Progressive_Disease"))
res.pca <- prcomp(train, scale = TRUE) #compute PCA
fviz_eig(res.pca) #visualize PCs
fviz_pca_ind(res.pca,
             col.ind = response_group, # Color by the quality of representation
             palette = c("#FF6666","cyan",  "gray"),
             repel = TRUE)    # Avoid text overlapping

#Create PCA Plot with ggplot2
df_out <- as.data.frame(res.pca$x)
df_out$Response <- normalised_annot_log2$Response
df_out$Treatment <- normalised_annot_log2$Treatment
head(df_out)

#goem point; color is the border, fill is the inside
fil.col <- c('#f768a1','#e6e600', '#1e90ff')
p<- ggplot(df_out,aes(x=PC1,y=PC2, label = rownames(df_out))) + 
  geom_point(aes(shape = Treatment, fill = Response), color = "black", size = 5, stroke = 1)+ theme_classic() + xlab("PC1 (24.1%)") + ylab("PC2 (13.9%)")+
  scale_fill_manual(values= fil.col)+ #change color of dots
  scale_shape_manual(values = c(21, 24)) + #need shape 21 and 24 circle and triangle with border
  guides(fill = guide_legend(override.aes = list(color = fil.col))) + #change legend
  #geom_text_repel(aes(label=rownames(df_out))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, axis.text=element_text(size=14),axis.title=element_text(size=14))

p
while (!is.null(dev.list()))  dev.off() #prevent pdf write problem
pdf("result/PCA.pdf", width = 7, height = 6)
p
dev.off()

####Plot Heatmap####
library(dplyr)
library(pheatmap)
library(readxl)

#Load normalised database
normalised_linear <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/NPC All Normalised Data.xlsx", 
    sheet = "Linear")
#Load annotations
annot_file <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/Nanotring Annotation All.xlsx")
annot_file <- annot_file %>% dplyr::select(SampleID, Site_Response, EBV, Smoking_Status, Sex) %>% as.data.frame() %>% tibble::column_to_rownames(var="SampleID")

#Creating separate database for Pre and Post
normalised_pre <- normalised_linear %>% select(Gene, contains("Pre")) %>% tibble::column_to_rownames(var="Gene") %>% as.matrix() %>% log2()
normalised_post <- normalised_linear %>% select(Gene, contains("Post")) %>% tibble::column_to_rownames(var="Gene") %>% as.matrix() %>% log2()

#Create Annotation color
ebv_col30k = c("#37a124","#e5c064")
names(ebv_col30k) = levels(as.factor(annot_file$EBV))
annot_color = list(
  Site_Response = c(PR= "#f768a1", PD="#1e90ff", SD ="#e6e600"), 
  EBV = ebv_col30k,
  Smoking_Status = c(Smoker= "#282828", `Ex-smoker` ="#767676", `Non-smoker` = "#f6f6f6"), Sex = c(Female ='plum', Male ='lightskyblue'))

hm <- pheatmap(normalised_pre,
               scale = "row", 
               show_rownames = FALSE, show_colnames = TRUE,
               border_color = NA,
               color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
               treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
               fontsize_row = 9,
               annotation_col = annot_file,
               annotation_colors = annot_color)

while (!is.null(dev.list()))  dev.off() #clear device
pdf("result/Pre_WholeGene_Heatmap.pdf", width = 5, height = 6)
hm
dev.off()

hmp <- pheatmap(normalised_post,
               scale = "row", 
               show_rownames = FALSE, show_colnames = TRUE,
               border_color = NA,
               color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
               treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
               fontsize_row = 9,
               annotation_col = annot_file,
               annotation_colors = annot_color)

while (!is.null(dev.list()))  dev.off() #clear device
pdf("result/Post_WholeGene_Heatmap.pdf", width = 5, height = 6)
hmp
dev.off()
