library(readxl)
library(dplyr)
require(ggpubr)
require(rstatix)
#Cell Type Profiling has been done using Nanostring nSolver Software
#Data is provided as Source Data
log2_cellscore <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/Relative Cell Type Score.xlsx", 
                                       sheet = "log2")
log2_cellscore <- log2_cellscore %>% tibble::column_to_rownames("Gene")
log2_cellscore <- t(log2_cellscore) %>% as.data.frame() %>% tibble::rownames_to_column(var = "SampleID")
annotationfile <- read_excel("C:/Users/chris/Dropbox/NCCS Personal/NPC IpiNivo Project/Nanostring & Other Analysis/NPC Batch 5/data table/Nanotring Annotation All.xlsx")
annotationfile <- select(annotationfile, c("SampleID", "Treatment","Site_Response")) %>% na.omit()

data_annot <- left_join(annotationfile, log2_cellscore)
data_annot$Site_Response <- factor(data_annot$Site_Response, levels = c("PR","SD","PD"))

#Plot
library(ggplot2)
fil.col <- c('#f768a1','#e6e600', '#1e90ff')
plot_bar <- function(cell_score, prepost) {
  if(prepost =="Pre"){ sub_label ="Pre-treatment"} else {sub_label ="On-treatment"}
  filtered <- dplyr::filter(data_annot, Treatment == prepost) 
  stat_form <- formula(paste0(cell_score, "~ Site_Response"))
  stat.test <- filtered %>% wilcox_test(stat_form, ref.group = "PR") #Perform t.test against reference group
  stat.test <- stat.test %>% add_xy_position(x = "Site_Response")
  print(stat.test)
  stat.test$y.position <- c(3.8,4.1)
  p <- ggplot(filtered, aes(x = Site_Response, y = !!ensym(cell_score))) + geom_boxplot() +
    geom_point(aes(fill = Site_Response), colour="black",pch=21, size = 3) +
    theme_classic() + labs(title = cell_score, subtitle = sub_label, y ="Score", x =NULL) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(face="bold"),
          legend.position = "none") + #remove legend
    scale_fill_manual(values= fil.col) +
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 3, djust= 2)
  while (!is.null(dev.list()))  dev.off() #clear device
  pdf(paste0("result/barplot_",cell_score, prepost,".pdf"), width = 2, height = 3)
  print(p)
  dev.off()
  print(p)
}
plot_bar("Treg.vs.TILs","Pre")
plot_bar("Treg.vs.TILs","On")
plot_bar("Total.TILs","Pre")
plot_bar("Total.TILs","On")
plot_bar("Total.TILs","Pre")
plot_bar("CD8.vs.Exhausted.CD8","Pre")
plot_bar("CD8.vs.Exhausted.CD8","On")
plot_bar("CD8.vs.Treg","Pre")
plot_bar("CD8.vs.Treg","On")
