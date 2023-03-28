### celltypes.analysis.R  #########################################################################
#
# Connie Li, 	November 2022

### DESCRIPTION ###################################################################################
#
# Look at nanostring-predicted cell type differences between treatment status and response to 
# treatment
#

### LIBRARIES, ENV & INPUT ########################################################################
library(DESeq2)
library(ggplot2)
library(MASS)
library(EDASeq)
library(pheatmap)
library(ggplot2)
library(gridExtra)
setwd('./Projects/IpiNivo/RSession')

### MAIN ##########################################################################################
celltype.data <- read.table(
	'./Projects/IpiNivo/data/nanostring/RelativeCellTypeScore.txt',
	header = TRUE,
	sep = '\t',
	row.names = 1
	)
colnames(celltype.data) <- gsub('X', '', colnames(celltype.data))
nanostring.annotation <- read.table(
	'./Projects/IpiNivo/data/annotation/NanotringAnnotationAll_2022_renewed.txt',
	sep ='\t',
	header = TRUE,
	quote = '',
	row.names=1
	)

wanted.samples <- intersect(colnames(celltype.data), rownames(nanostring.annotation))
celltype.data <- celltype.data[,wanted.samples]
nanostring.annotation <- nanostring.annotation[wanted.samples,]

celltype.results <- data.frame(
	pre.pd.vs.sd.p = rep(NA, nrow(celltype.data)),
	pre.pd.vs.sd.q = rep(NA, nrow(celltype.data)),
	pre.pd.vs.pr.p = rep(NA, nrow(celltype.data)),
	pre.pd.vs.pr.q = rep(NA, nrow(celltype.data)),
	pre.sd.vs.pr.p = rep(NA, nrow(celltype.data)),
	pre.sd.vs.pr.q = rep(NA, nrow(celltype.data)),

	on.pd.vs.sd.p = rep(NA, nrow(celltype.data)),
	on.pd.vs.sd.q = rep(NA, nrow(celltype.data)),
	on.pd.vs.pr.p = rep(NA, nrow(celltype.data)),
	on.pd.vs.pr.q = rep(NA, nrow(celltype.data)),
	on.sd.vs.pr.p = rep(NA, nrow(celltype.data)),
	on.sd.vs.pr.q = rep(NA, nrow(celltype.data)),
	row.names = rownames(celltype.data)
	)
nanostring.annotation$SampleNCC <- gsub('NCC', '', nanostring.annotation$SampleNCC)

pd.samples <- nanostring.annotation[which(nanostring.annotation$Site_Response == 'PD'), 'SampleNCC']
sd.samples <- nanostring.annotation[which(nanostring.annotation$Site_Response == 'SD'), 'SampleNCC']
pr.samples <- nanostring.annotation[which(nanostring.annotation$Site_Response == 'PR'), 'SampleNCC']

for (each.celltype in rownames(celltype.data)) {
	pre.data <- celltype.data[each.celltype, grep('Pre', colnames(celltype.data))]
	post.data <- celltype.data[each.celltype, grep('Post', colnames(celltype.data))]

	pre.pd.vs.sd <- wilcox.test(
		as.numeric(pre.data[which(colnames(pre.data) %in% paste0(pd.samples, '_Pre'))]),
		as.numeric(pre.data[which(colnames(pre.data) %in% paste0(sd.samples, '_Pre'))])
		)
	pre.pd.vs.pr <- wilcox.test(
		as.numeric(pre.data[which(colnames(pre.data) %in% paste0(pd.samples, '_Pre'))]),
		as.numeric(pre.data[which(colnames(pre.data) %in% paste0(pr.samples, '_Pre'))])
		)
	pre.sd.vs.pr <- wilcox.test(
		as.numeric(pre.data[which(colnames(pre.data) %in% paste0(sd.samples, '_Pre'))]),
		as.numeric(pre.data[which(colnames(pre.data) %in% paste0(pr.samples, '_Pre'))])
		)
	post.pd.vs.sd <- wilcox.test(
		as.numeric(post.data[which(colnames(post.data) %in% paste0(pd.samples, '_Post'))]),
		as.numeric(post.data[which(colnames(post.data) %in% paste0(sd.samples, '_Post'))])
		)
	post.pd.vs.pr <- wilcox.test(
		as.numeric(post.data[which(colnames(post.data) %in% paste0(pd.samples, '_Post'))]),
		as.numeric(post.data[which(colnames(post.data) %in% paste0(pr.samples, '_Post'))])
		)
	post.sd.vs.pr <- wilcox.test(
		as.numeric(post.data[which(colnames(post.data) %in% paste0(sd.samples, '_Post'))]),
		as.numeric(post.data[which(colnames(post.data) %in% paste0(pr.samples, '_Post'))])
		)
	celltype.results[each.celltype,] <- c(
		pre.pd.vs.sd$p.value,
		NA,
		pre.pd.vs.pr$p.value,
		NA,
		pre.sd.vs.pr$p.value,
		NA,
		post.pd.vs.sd$p.value,
		NA,
		post.pd.vs.pr$p.value,
		NA,
		post.sd.vs.pr$p.value,
		NA
		)
}
celltype.results$pre.pd.vs.sd.q <- p.adjust(celltype.results$pre.pd.vs.sd.p, method='fdr')
celltype.results$pre.pd.vs.pr.q <- p.adjust(celltype.results$pre.pd.vs.pr.p, method='fdr')
celltype.results$pre.sd.vs.pr.q <- p.adjust(celltype.results$pre.sd.vs.pr.p, method='fdr')
celltype.results$on.pd.vs.sd.q <- p.adjust(celltype.results$on.pd.vs.sd.p, method='fdr')
celltype.results$on.pd.vs.pr.q <- p.adjust(celltype.results$on.pd.vs.pr.p, method='fdr')
celltype.results$on.sd.vs.pr.q <- p.adjust(celltype.results$on.sd.vs.pr.p, method='fdr')

write.table(
	celltype.results,
	file = 'ipinivo_nanostringcelltype.utest.results.txt',
	sep = '\t',
	quote = FALSE
	)
