### post.pre.signatures.analysis.R  ###############################################################
#
# Connie Li, 	November 2020 

### DESCRIPTION ###################################################################################
#
# Look at performance in post-signature in pre-samples and pre-signature in post-samples 

### LIBRARIES, ENV & INPUT ########################################################################
library(DESeq2)
library(ggplot2)
library(MASS)
library(EDASeq)
library(pheatmap)
library(ggplot2)
library(gridExtra)

### MAIN ##########################################################################################

eset <- readRDS('./Projects/IpiNivo/RSession/2022-10-19-ipinivo.nanostring.ruvgcorrected.esetobj.rds')
de.results <- read.table(
	'./Projects/IpiNivo/data/de_results/final_list_Post_Site_PDPR_with_fc_pvalue_Used_for_manuscript_first_submission.csv',
	fill = TRUE,
	sep = ',',
	header = TRUE
	)
nanostring.log2 <- read.table(
	'./Projects/IpiNivo/data/nanostring/nanostring.log2.norm.data.txt',
	sep = '\t',
	header = TRUE,
	row.names = 1
	)
colnames(nanostring.log2) <- gsub('X', 'NCC', colnames(nanostring.log2))

nanostring.annotation <- read.table(
	'./Projects/IpiNivo/data/annotation/NanotringAnnotationAll_2022_renewed.txt',
	sep ='\t',
	header = TRUE,
	quote = '',
	row.names=1
	)
rownames(nanostring.annotation) <- paste0('NCC', rownames(nanostring.annotation))

pretreatment.genes <- c('HLA-A', 'ICAM5', 'DKK1', 'COL6A3', 'ADAM12', 'SFRP1', 'PIK3R1')
overlapping.genes <- intersect(de.results$Gene, pretreatment.genes)
genes.of.interest <- unique(c(de.results$Gene, pretreatment.genes))

# U test results 
pre.on.utest.results <- data.frame(
	sample.type = NA,
	pre.median.pd = rep(NA, length(genes.of.interest)),
	pre.median.pr = rep(NA, length(genes.of.interest)),
	pre.median.diff = rep(NA, length(genes.of.interest)),
	pre.p.value = rep(NA, length(genes.of.interest)),
	pre.q.value = rep(NA, length(genes.of.interest)),	
	on.median.pd = rep(NA, length(genes.of.interest)),
	on.median.pr = rep(NA, length(genes.of.interest)),
	on.median.diff = rep(NA, length(genes.of.interest)),
	on.p.value = rep(NA, length(genes.of.interest)),
	on.q.value = rep(NA, length(genes.of.interest)),	
	row.names = genes.of.interest
	)
pre.on.utest.results[de.results$Gene, 'sample.type'] <- 'on-treatment'
pre.on.utest.results[pretreatment.genes, 'sample.type'] <- 'pre-treatment'
pre.on.utest.results[overlapping.genes, 'sample.type'] <- 'both'
pre.pd.samples <- rownames(nanostring.annotation)[which(nanostring.annotation$Treatment == 'Pre' & nanostring.annotation$Site_Response == 'PD')]
pre.pr.samples <- rownames(nanostring.annotation)[which(nanostring.annotation$Treatment == 'Pre' & nanostring.annotation$Site_Response == 'PR')]
post.pd.samples <- rownames(nanostring.annotation)[which(nanostring.annotation$Treatment == 'On' & nanostring.annotation$Site_Response == 'PD')]
post.pr.samples <- rownames(nanostring.annotation)[which(nanostring.annotation$Treatment == 'On' & nanostring.annotation$Site_Response == 'PR')]

wanted.genes <- overlapping.genes
wanted.genes <- rownames(pre.on.utest.results)

for (each.gene in wanted.genes) {
	pretreatment.utest <- wilcox.test(
		as.numeric(nanostring.log2[each.gene, pre.pd.samples]),
		as.numeric(nanostring.log2[each.gene, pre.pr.samples])
		)
	ontreatment.utest <- wilcox.test(
		as.numeric(nanostring.log2[each.gene, post.pd.samples]),
		as.numeric(nanostring.log2[each.gene, post.pr.samples])
		)
	pre.on.utest.results[each.gene,-1] <- c(
		median(as.numeric(nanostring.log2[each.gene, pre.pd.samples])),
		median(as.numeric(nanostring.log2[each.gene, pre.pr.samples])),
		NA,
		pretreatment.utest$p.value,
		NA,

		median(as.numeric(nanostring.log2[each.gene, post.pd.samples])),
		median(as.numeric(nanostring.log2[each.gene, post.pr.samples])),
		NA,
		ontreatment.utest$p.value,
		NA
		)
}

pre.on.utest.results$pre.median.diff <- pre.on.utest.results$pre.median.pr - pre.on.utest.results$pre.median.pd
pre.on.utest.results$on.median.diff <- pre.on.utest.results$on.median.pr - pre.on.utest.results$on.median.pd

pre.on.utest.results$pre.q.value <- p.adjust(pre.on.utest.results$pre.p.value)
pre.on.utest.results$on.q.value <- p.adjust(pre.on.utest.results$on.p.value)
pre.on.utest.results$gene <- rownames(pre.on.utest.results)

write.table(
	pre.on.utest.results,
	file = paste0(Sys.Date(), '-pre.in.post_post.in.pre_utests_lollipop.txt'),
	sep = '\t',
	quote = FALSE,
	row.names = FALSE
	)

top.twenty.sig <- de.results[order(de.results$pvalue.x), 'Gene'][1:20]
top.twenty.sig <- pre.on.utest.results[order(pre.on.utest.results$on.p.value),'gene'][1:20]
overlapping.genes <- c('HLA-A', 'ICAM5', 'DKK1', 'COL6A3', 'ADAM12', 'SFRP1', 'PIK3R1')
wanted.genes <- overlapping.genes

point.col <- rep('#5F1415', length(wanted.genes))
point.col[which(pre.on.utest.results[wanted.genes, 'on.median.diff'] < 0)] <- '#002F70'

point.size <- 2*-log10(pre.on.utest.results[wanted.genes, 'on.p.value'])
p1 <- ggplot(
	pre.on.utest.results[wanted.genes,], 
	aes(y=gene, x=on.median.diff)) +
	theme_classic() +
	xlab("On-treatment") +
	scale_x_continuous(limits = c(-2, 4.5)) +
	theme(axis.title.y=element_blank()) +
	geom_point(size = point.size, color = point.col, alpha = 0.9) +
	geom_segment(aes(y=gene,yend=gene, x=0, xend=on.median.diff)) +
	geom_vline(xintercept=0,size=1)
point.size <- 2*-log10(pre.on.utest.results[wanted.genes, 'pre.p.value'])
p2 <- ggplot(
	pre.on.utest.results[wanted.genes,], 
	aes(y=gene, x=pre.median.diff)
	) +
	theme_classic() +
	xlab("Pre-treatment") +
	scale_x_continuous(limits = c(-2, 4.5)) +
	theme(axis.title.y=element_blank()) +
	geom_point(size = point.size, color = point.col, alpha = 0.9) +
	geom_segment(aes(y=gene, yend=gene, x=0, xend=pre.median.diff)) +
	geom_vline(xintercept=0,size=1) +
	theme(legend.position="right") +
	scale_color_manual(values = c('orange', 'purple'))
pdf(paste0(Sys.Date(), "lollipop.presig.pdf"), width = 4.5, height = 1.7) #4.5 for on-treatment
grid.arrange(p1, p2, ncol=2)
dev.off()

legend.data <- data.frame(
	x = rep(1, 6),
	y = c(0.2, 0.4, 0.6, 0.8, 1, 1.2),
	size = c(-6, -4, -2, 2, 4, 6),
	colour = c('#002F70', '#002F70', '#002F70','#5F1415',  '#5F1415', '#5F1415')
	)
pdf(paste0(Sys.Date(), "lollipop.legend.pdf"), width = 4.5, height = 4.5)
ggplot(
	legend.data, 
	aes(y=y, x=x)
	) +
	theme_classic() +
	scale_y_continuous(limits = c(-0.5, 2)) +
	geom_point(size = abs(legend.data$size), color = legend.data$colour, alpha = 0.9)
dev.off()

patient.data <- pData(eset)
patient.data$time.site.response <- factor(paste0(patient.data$Treatment, '-', patient.data$Site_Response))
patient.data$Treatment <- factor(patient.data$Treatment, levels = c('Pre', 'On'))
patient.data$SampleNCC <- factor(patient.data$SampleNCC)
patient.data$Site_Response <- factor(patient.data$Site_Response, levels = c('PD', 'PR'))
patient.coldata <- patient.data[,c('Treatment', 'SampleNCC', 'Site_Response', 'W_1')]

keep.samples <- rownames(patient.data)[which(patient.data$Site_Response %in% c('PD', 'PR'))]

# Look at interaction
ns.dds <- DESeqDataSetFromMatrix(
	countData = counts(eset)[,keep.samples],
	colData = patient.coldata[keep.samples,],
	design = formula('~ Site_Response * Treatment + W_1')
	)
ns.dds <- DESeq(ns.dds)
ns.vsd <- varianceStabilizingTransformation(ns.dds)
normalised.ns <- assay(ns.vsd)

post.de.genes <- intersect(de.results$New, de.results$Padj)

# annotation
on.treatment.cov <- patient.data[which(patient.data$Treatment == 'On' & patient.data$Site_Response %in% c('PD', 'PR')), c('Site_Response', 'EBV', 'Smoking_Status', 'Sex')]
ebv_col30k = c("#37a124","#e5c064")
names(ebv_col30k) = levels(as.factor(on.treatment.cov $EBV))
annot_color = list(
	Site_Response = c(PR= "#f768a1", PD="#1e90ff", SD ="#e6e600"), 
	EBV = ebv_col30k,
	Smoking_Status = c(Smoker= "#282828", `Ex-smoker` ="#767676", `Non-smoker` = "#f6f6f6"), 
	Sex = c(Female ='plum', Male ='lightskyblue')
	)
on.treatment.data <- normalised.ns[post.de.genes,rownames(patient.data)[which(patient.data$Treatment == 'On' & patient.data$Site_Response %in% c('PD', 'PR'))]]

pheatmap(
	on.treatment.data ,
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = on.treatment.cov,
	annotation_colors = annot_color
)

pre.treatment.data <- normalised.ns[post.de.genes,rownames(patient.data)[which(patient.data$Treatment == 'Pre' & patient.data$Site_Response %in% c('PD', 'PR'))]]
pre.treatment.cov <- patient.data[which(patient.data$Treatment == 'Pre' & patient.data$Site_Response %in% c('PD', 'PR')), c('Site_Response', 'EBV', 'Smoking_Status', 'Sex')]
pdf(paste0(Sys.Date(), "on_Site_PRPD_postsig_heatmap_unlabelled.pdf"), width = 4.5, height = 6)
pheatmap(
	pre.treatment.data,
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = pre.treatment.cov,
	annotation_colors = annot_color
)
dev.off()

# Look at interaction
ns.dds <- DESeqDataSetFromMatrix(
	countData = counts(ruv.results$eset)[,keep.samples],
	colData = patient.coldata[keep.samples,],
	design = formula('~ Site_Response * Treatment + W_1')
	)

ns.dds <- DESeq(ns.dds)
# Look at relevant contrasts
resultsNames(ns.dds)
results(ns.dds, name = "Site_Response_PR_vs_PD")
results(ns.dds, name = "Treatment_On_vs_Pre")
results(ns.dds, name = "Site_ResponsePR.TreatmentOn")

# Question: Do genes whose expression changes between Pre->On have an association with Site response?
ns.dds <- DESeqDataSetFromMatrix(
	countData = counts(ruv.results$eset)[,keep.samples],
	colData = patient.coldata[keep.samples,],
	design = formula('~ Treatment + W_1')
	)
ns.dds <- DESeq(ns.dds)
# Look at relevant contrasts
resultsNames(ns.dds)
on.vs.pre.results <- results(ns.dds, name = "Treatment_On_vs_Pre")
on.vs.pre.sig.results <- on.vs.pre.results[which(on.vs.pre.results$padj < 0.05),]

# and are those genes associated with treatment?
ns.dds <- DESeqDataSetFromMatrix(
	countData = counts(ruv.results$eset)[,keep.samples],
	colData = patient.coldata[keep.samples,],
	design = formula('~ Site_Response + W_1')
	)
ns.dds <- DESeq(ns.dds)
# Look at relevant contrasts
pd.vs.pr.results <- results(ns.dds, name = "Site_Response_PR_vs_PD")
onvspre.in.pdvspr <- pd.vs.pr.results[rownames(on.vs.pre.sig.results),]

prevson.genechange.plot <- list()
for (each.gene in rownames(on.vs.pre.sig.results)) {
	change.graph.data <- data.frame(
		mrna = normalised.ns[each.gene, c(pre.samples, on.samples)],
		pre.or.on = factor(c(rep('Pre', length(pre.samples)), rep('Post', length(on.samples))), levels=c('Pre', 'Post')),
		same.pt = c(1:length(pre.samples), 1:length(on.samples)),
		response = patient.data[c(pre.samples, on.samples), 'Site_Response']
		)
	prevson.genechange.plot[[each.gene]] <- ggplot(
		change.graph.data, 
		aes(pre.or.on, mrna)
		) +
		xlab(each.gene) +
		geom_boxplot() +
		geom_jitter() +
		theme_classic() +
		guides(color = 'none', size = 'none')
}
for (each.gene in rownames(on.vs.pre.sig.results)) {
	change.graph.data <- data.frame(
		mrna = normalised.ns[each.gene, c(pre.samples, on.samples)],
		pre.or.on = factor(c(rep('Pre', length(pre.samples)), rep('Post', length(on.samples))), levels=c('Pre', 'Post')),
		same.pt = c(1:length(pre.samples), 1:length(on.samples)),
		response = patient.data[c(pre.samples, on.samples), 'Site_Response']
		)
	prevson.genechange.plot[[paste0(each.gene, '2')]] <- ggplot(
		change.graph.data, 
		aes(response, mrna)
		) +
		xlab(each.gene) +
		geom_boxplot() +
		geom_jitter() +
		theme_classic() +
		guides(color = 'none', size = 'none')
}
pdf("prevson_boxplot.pdf", width = 6, height = 3)
ggarrange(
	plots = prevson.genechange.plot,
	nrow = 2,
	ncol = 3,
	common.legend = TRUE,
	legend = "bottom"
	)
dev.off()

# Finally, a heatmap of all three
patient.cov <- pData(ruv.results$eset)[keep.samples, c('Site_Response', 'EBV', 'Smoking_Status', 'Sex', 'Treatment')]
pheatmap(
	normalised.ns[rownames(on.vs.pre.sig.results),keep.samples],
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = patient.cov,
	annotation_colors = annot_color
)
# CONCLUSION: Not sufficiently powered for this anslysis. But, have one gene to follow up on (the gene that is DE between prevson and seen in lisda's sig)

# Find DE genes for pre-treatment 
ns.dds <- DESeqDataSetFromMatrix(
	countData = counts(ruv.results$eset)[,keep.samples],
	colData = patient.coldata[keep.samples,],
	design = formula('~ Site_Response * Treatment + W_1')
	)

ns.dds <- DESeq(ns.dds)
# Look at relevant contrasts
resultsNames(ns.dds)
results(ns.dds, name = "Site_Response_PR_vs_PD")
results(ns.dds, name = "Treatment_On_vs_Pre")
results(ns.dds, name = "Site_ResponsePR.TreatmentOn")

# One gene is sig associated with response in pre-treatment mRNA
results(ns.dds, name = "Site_Response_PR_vs_PD")[which(results(ns.dds, name = "Site_Response_PR_vs_PD")$padj< 0.05),] 

# Mimic Lisda's leave one sample out strategy
patient.coldata <- patient.data[,c('Treatment', 'SampleNCC', 'Site_Response', 'W_1')]
pretreatment.de.loo.results <- list()
for (each.sample in keep.samples) {
	ns.dds <- DESeqDataSetFromMatrix(
		countData = counts(ruv.results$eset)[,keep.samples[-which(keep.samples == each.sample)]],
		colData = patient.coldata[keep.samples[-which(keep.samples == each.sample)],],
		design = formula('~ Site_Response * Treatment + W_1')
		)
	ns.dds <- DESeq(ns.dds)
	pretreatment.de.loo.results[[each.sample]] <- results(ns.dds, name = "Site_Response_PR_vs_PD")
}
tabulate.pretreatment.loo.results <- data.frame(
	ntimes.adjsig = rowSums(do.call(cbind, lapply(pretreatment.de.loo.results, function(x){x[,'padj'] < 0.05}))),
	ntimes.sig = rowSums(do.call(cbind, lapply(pretreatment.de.loo.results, function(x){x[,'pvalue'] < 0.01}))),
	mean.p = rowMeans(do.call(cbind, lapply(pretreatment.de.loo.results, function(x){x[,'pvalue']}))),
	mean.fc = rowMeans(do.call(cbind, lapply(pretreatment.de.loo.results, function(x){x[,'log2FoldChange']}))),
	row.names = row.names(pretreatment.de.loo.results[[1]])
	)
tabulate.pretreatment.loo.results <- tabulate.pretreatment.loo.results[order(tabulate.pretreatment.loo.results$ntimes.adjsig, decreasing=T),]
pre.de.genes <- rownames(tabulate.pretreatment.loo.results)[which(tabulate.pretreatment.loo.results$ntimes.sig ==29)]

pre.treatment.data <- normalised.ns[pre.de.genes,rownames(patient.data)[which(patient.data$Treatment == 'Pre' & patient.data$Site_Response %in% c('PD', 'PR'))]]
pre.treatment.cov <- patient.data[which(patient.data$Treatment == 'Pre' & patient.data$Site_Response %in% c('PD', 'PR')), c('Site_Response', 'EBV', 'Smoking_Status', 'Sex')]
pdf("pre_Site_PRPD_presig_heatmap_unlabelled.pdf", width = 4.5, height = 6)
pheatmap(
	pre.treatment.data,
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	#labels_col = gsub('NCC', '', colnames(pre.treatment.data)),
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = pre.treatment.cov,
	annotation_colors = annot_color
)
dev.off()

on.treatment.data <- normalised.ns[pre.de.genes,rownames(patient.data)[which(patient.data$Treatment == 'On' & patient.data$Site_Response %in% c('PD', 'PR'))]]
on.treatment.cov <- patient.data[which(patient.data$Treatment == 'On' & patient.data$Site_Response %in% c('PD', 'PR')), c('Site_Response', 'EBV', 'Smoking_Status', 'Sex')]
pdf("pre_Site_PRPD_presig_heatmap_unlabelled.pdf", width = 4.5, height = 6)
pheatmap(
	on.treatment.data,
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	#labels_col = gsub('NCC', '', colnames(pre.treatment.data)),
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = on.treatment.cov,
	annotation_colors = annot_color
)
dev.off()


lisda.post.de <- read.table(
	'./Projects/IpiNivo/data/fromLisda/de_results/full_list_Post_Site_PDPR_with_fc_pvalue.csv',
	sep = ',',
	header = TRUE
	)
post.de.genes <- lisda.post.de$Gene[which(lisda.post.de$freq == 14)]

# annotation
on.treatment.cov <- patient.data[which(patient.data$Treatment == 'On' & patient.data$Site_Response %in% c('PD', 'PR')), c('Site_Response', 'EBV', 'Smoking_Status', 'Sex')]
ebv_col30k = c("#37a124","#e5c064")
names(ebv_col30k) = levels(as.factor(on.treatment.cov $EBV))
annot_color = list(
	Site_Response = c(PR= "#f768a1", PD="#1e90ff", SD ="#e6e600"), 
	EBV = ebv_col30k,
	Smoking_Status = c(Smoker= "#282828", `Ex-smoker` ="#767676", `Non-smoker` = "#f6f6f6"), 
	Sex = c(Female ='plum', Male ='lightskyblue')
	)
on.treatment.data <- normalised.ns[post.de.genes,rownames(patient.data)[which(patient.data$Treatment == 'On' & patient.data$Site_Response %in% c('PD', 'PR'))]]
pheatmap(
	on.treatment.data ,
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = on.treatment.cov,
	annotation_colors = annot_color
)
post.de.genes <- lisda.post.de$Gene[which(lisda.post.de$freq == 14)]
pre.treatment.data <- normalised.ns[post.de.genes,rownames(patient.data)[which(patient.data$Treatment == 'Pre' & patient.data$Site_Response %in% c('PD', 'PR'))]]
pre.treatment.cov <- patient.data[which(patient.data$Treatment == 'Pre' & patient.data$Site_Response %in% c('PD', 'PR')), c('Site_Response', 'EBV', 'Smoking_Status', 'Sex')]

pdf("on_Site_PRPD_postsig_heatmap_unlabelled.pdf", width = 4.5, height = 6)
pheatmap(
	pre.treatment.data,
	scale = "row", 
	show_rownames = FALSE, 
	show_colnames = TRUE,
	#labels_col = gsub('NCC', '', colnames(pre.treatment.data)),
	border_color = NA,
	color = colorspace::diverging_hcl(n = 100, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), register = ),
	treeheight_row = 0, treeheight_col =1,#cellwidth = 10, cellheight = 8,
	fontsize_row = 9,
	annotation_col = pre.treatment.cov,
	annotation_colors = annot_color
)
dev.off()
