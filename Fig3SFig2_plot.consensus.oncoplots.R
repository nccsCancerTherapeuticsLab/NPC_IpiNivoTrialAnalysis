### plot.consensus.oncoplots.R  ##################################################################
#
# Connie Li, 	November 2020 

### DESCRIPTION ###################################################################################
#
# Make co-mutation and mutational signatures plots from consensus calls (Figure 3B-C)

### LIBRARIES, ENV & INPUT ########################################################################
library(BoutrosLab.plotting.general)
library(deconstructSigs)
library(maftools)
source('maftools.functions.R')
source('plot.consensus.sigs.R')

### FUNCTIONS  ####################################################################################
recode.mutations.numeric <- function(mutation.matrix, keep.silent = TRUE, keep.igr=TRUE) {
	unique.index <- 1	
	mutation.labels <- character()
	possible.types <- c('Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del|Frame_Shift_Ins', 'In_Frame_Del|In_Frame_Ins', 'Splice_Region|Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site|Translation_Start_Site In_Frame_Del|Translation_Start_Site In_Frame_Ins','3\'Flank|3\'UTR|5\'Flank|5\'UTR', 'RNA', 'Silent', 'IGR', 'Intron')
	if (!keep.silent) {
		possible.types <- possible.types[-which(possible.types == 'Silent')]
	}
	if (!keep.igr) {
		possible.types <- possible.types[-which(possible.types %in% c('IGR', 'Intron'))]
	}	
	recoded.matrix <- matrix(NA, nrow=nrow(mutation.matrix), ncol=ncol(mutation.matrix))
	for (each.type in possible.types) {
		if (any(grepl(each.type, mutation.matrix))) {
			for (i in 1:nrow(mutation.matrix)) {
				each.row <- mutation.matrix[i,]
				recoded.matrix[i, grep(each.type, each.row)] <- unique.index
			}
			mutation.labels[unique.index] <- each.type
			unique.index <- unique.index + 1
		}
	}
	recoded.matrix <- as.data.frame(apply(recoded.matrix, 2, as.numeric))
	rownames(recoded.matrix) <- rownames(mutation.matrix)
	colnames(recoded.matrix) <- colnames(mutation.matrix)
	# Then assign the colour 
	if (length(mutation.labels) <= 7) {
		avail.colours <- c("#336A90", "#65B4A2", "#B1D39A", "#F3E0A6", "#E0A890", "#FFE1EE", "#713E5A")
	} else {
		avail.colours <- default.colours(length(mutation.labels), palette.type='pastel')
	}
	mutation.colours.labels <- data.frame(
		colours = avail.colours[1:length(mutation.labels)],
		labels = mutation.labels
		)
	mutation.colours.labels$labels[grep('Frame_Shift', mutation.colours.labels$labels)] <- 'Frame Shift Indel'
	mutation.colours.labels$labels[grep('Translation', mutation.colours.labels$labels)] <- 'Translation Start Site'
	mutation.colours.labels$labels[grep('In_Frame_Del', mutation.colours.labels$labels)] <- 'In Frame Indel'
	mutation.colours.labels$labels[grep('Splice', mutation.colours.labels$labels)] <- 'Splice Site'
	mutation.colours.labels$labels[grep('Flank', mutation.colours.labels$labels)] <- "3\'/5\'UTR"
	mutation.colours.labels$labels <- gsub('_', ' ', mutation.colours.labels$labels)
	return(list(mutation.matrix = recoded.matrix, mutation.colour.labels = mutation.colours.labels))
}

### MAIN ####################################################################################

sample.list <- read.table("samplelist")
sample.list.no23 <- sample.list[-which(sample.list[,1] %in% c('NCC023')),]
sample.list <- sample.list[,1]
consensus.maf.list <- read.mafs(
	mafs.pwd = "./Projects/IpiNivo/data/mafFiles/consensus_mafs/",
	sample.list = sample.list.no23,
	tag = 'consensus'
	)
consensus.mafs <- merge_mafs(
	consensus.maf.list,
	)
gene.sample.matrix <- mutCountMatrix(
	consensus.mafs,
	includeSyn = FALSE,
	countOnly = NULL,
	removeNonMutated = TRUE
	)
clinical.data <- read.table(
	'./Projects/IpiNivo/data/annotation/SampleInformation.txt',
	sep ='\t',
	header = TRUE,
	quote = ''
	)
nanostring.data <- read.table(
	"./Projects/IpiNivo/data/nanostring_annotation.csv",
	sep =',',
	header = TRUE,
	row.names=1
	)
nanostring.data <- nanostring.data[grep('Pre', rownames(nanostring.data)),]
rownames(nanostring.data) <- gsub('Pre', '', nanostring.data$Alt_SampleID)
maf.list <- list.files('./Projects/IpiNivo/data/mafFiles/consensus_mafs/')

# Get gene lists
flag.genes <- read.table(
	"./Projects/IpiNivo/reference/FLAGS.txt",
	sep = '\t',
	header = TRUE
	)
cosmic.genes <- read.table(
	"./Projects/IpiNivo/reference/CosmicMutantExportCensus.tsv.gz",
	sep = '\t',
	header = TRUE
	)
nanostring.de.genes <- read.table(
	"./Projects/IpiNivo/data/fromLisda/DE_List_Ontreatment.csv",
	sep = ',',
	header = TRUE
	)
ade.npc.genes <- as.vector(t(read.table(
	'./RefFiles/NPC/npc.genelist.txt',
	sep = '\t',
	header = FALSE
	)))
dai.genes <- c('EGFR', 'PIK3CA', 'KRAS', 'HRAS', 'NRAS', 'BRAF', 'KIT', 'ABL1', 'AKT1', 'AKT2', 'CDK4', 'ERBB2', 'FGFR1', 'FGFR3', 'FLT3', 'JAK2', 'MET', 'RET', 'PDGFRA', 'ARID1A', 'BAP1', 'KMT2B', 'KMT2C', 'KMT2D', 'TSHZ3', 'HDAC4', 'PAXIP1', 'FAT1', 'FAT2', 'FAT3', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'ATG2A', 'ATG7', 'ATG13')
tcga.genes <- c('CDKN2A', 'FAT1', 'TP53','CASP8', 'AJUBA', 'PIK3CA', 'NOTCH1', 'KMT2D', 'NSD1', 'HLA-A', 'TGFBR2', 'HRAS', 'FBXW7', 'RB1', 'PIK3R1', 'TRAF3', 'NFE2L2', 'CUL3', 'PTEN')
genome.paper.genes <- c('TRAF3', 'NFKBIA', 'AEBP1', 'NLRC5', 'HLA-A', 'PLIN4', 'MUC21', 'SLC35G5', 'ERVW-1')
# http://cco.amegroups.com/article/view/9478/10771
de.genes <- nanostring.de.genes[,'Gene']
cosmic.genes <- unique(cosmic.genes$Gene.name)
flag.genes <- flag.genes[1:100,'FLAGS']

# Pick wanted genes: 
#1. Recurrently mutated in cohort
wanted.genes <- sort(rowSums(gene.sample.matrix>0),decreasing=TRUE)
wanted.genes <- names(wanted.genes[wanted.genes >= 2])
#2. Not FLAGs 
wanted.genes <- wanted.genes[which(!wanted.genes %in% flag.genes)]
#3. Are cancer/NPC genes 
wanted.genes <- intersect(wanted.genes, unique(c(cosmic.genes, dai.genes, ade.npc.genes, genome.paper.genes)))

# Fill in data frame with mutations information 
gene.sample.matrix.var <- data.frame(matrix(NA, nrow=length(wanted.genes), ncol=length(sample.list)))
rownames(gene.sample.matrix.var) <- wanted.genes
colnames(gene.sample.matrix.var) <- sample.list
for (each.sample in colnames(gene.sample.matrix.var)) {
	sample.maf <- read.table(
		paste0("./Projects/IpiNivo/data/mafFiles/consensus_mafs/", maf.list[grep(each.sample, maf.list)]),
		sep = '\t',
		header = TRUE,
		quote = ''
		)
	for (each.gene in rownames(gene.sample.matrix.var)) {
		gene.maf <- sample.maf[which(sample.maf$Hugo_Symbol == each.gene),]
		if (nrow(gene.maf) < 1) {next()}
		variant.class <- gene.maf$Variant_Classification
		if (length(variant.class) > 1) {
			print(paste0(each.sample, ' has ', length(variant.class), ' variants in ', each.gene))
			variant.class <- paste0(variant.class, collapse=',')
		} 
		gene.sample.matrix.var[each.gene, each.sample] <- variant.class
	} 
}
gene.sample.matrix.var2 <- gene.sample.matrix.var
# Prioritze mutation classes and convert to numeric 
recode.mutations <- recode.mutations.numeric(gene.sample.matrix.var, keep.igr = FALSE, keep.silent = FALSE)
gene.sample.matrix.var <- recode.mutations[[1]]
gene.sample.matrix.var <- gene.sample.matrix.var[nrow(gene.sample.matrix.var):1,]
mutation.colours.labels <- recode.mutations[[2]]
# Remove empty rows
gene.sample.matrix.var <- gene.sample.matrix.var[!(rowSums(is.na(gene.sample.matrix.var))==ncol(gene.sample.matrix.var)),]

# Put together annotation data
annotation.data <- clinical.data
rownames(annotation.data) <- annotation.data$Subject_no
annotation.data <- annotation.data[colnames(gene.sample.matrix.var),c(3,4,7,9,10,11,12,17,18)]
annotation.data$Best_overall_response <- factor(annotation.data$Best_overall_response, levels=c('PR', 'SD', 'PD'))
annotation.data$Site_Response <- nanostring.data[rownames(annotation.data), 'Site_Response']

# Recode Best_overall_response to Site Response
annotation.data$Best_overall_response <- factor(annotation.data$Site_Response,  levels=c('PR', 'SD', 'PD'))
annotation.data$Sex <- factor(annotation.data$Sex, levels=c('M', 'F'))
annotation.data$Smoking_Status <- factor(annotation.data$Smoking_Status, levels=c('Non-smoker', 'Ex-smoker', 'Smoker'))
# EV data
annotation.data$Baseline_EBV_DNA_copies.ml_ <- as.numeric(gsub('<|>', '', annotation.data$Baseline_EBV_DNA_copies.ml_))
annotation.data$Baseline_EBV_DNA_SI_ <- as.numeric(gsub('<|>', '', annotation.data$Baseline_EBV_DNA_SI_))
annotation.data$EBV <- factor(annotation.data$Baseline_EBV_DNA_SI_ > 3500)
annotation.data$EBV <- factor(annotation.data$Baseline_EBV_DNA_SI_ > 7800)

overlapping.samples <- intersect(colnames(gene.sample.matrix.var), rownames(annotation.data))
categorical.annotation <- annotation.data[overlapping.samples,c('Best_overall_response', 'EBV', 'Smoking_Status', 'Sex')]
levels(categorical.annotation$Sex) <- c(2,1)
levels(categorical.annotation$Best_overall_response) <- c(3,4,5)
levels(categorical.annotation$Smoking_Status) <- c(6,7,8)
levels(categorical.annotation$EBV) <- c(9,10)

# Calculated without silent mutations
tmb.maftools <- tmb(consensus.mafs, captureSize = 60)
tmb.data <- tmb.maftools$total_perMB
names(tmb.data) <- gsub('DNA_', '', tmb.maftools$Tumor_Sample_Barcode)
tmb.data['NCC023'] <- 0
compare.tmb.prevstudy(tmb.data)

#Set gene and sample order
# First, some ordering for genes by recurrence and pts 
ordering.matrix <- as.data.frame(cbind(categorical.annotation[colnames(gene.sample.matrix.var),], tmb = tmb.data[colnames(gene.sample.matrix.var)], t(gene.sample.matrix.var)))
gene.order <- c('TP53', 'FAM135B', 'CSMD3', 'MUC4', 'TRAF3', 'APC', 'PRKDC', 'CTNNA2', 'NOS3', 'NBEA', 'ATG2A', 'MYO5A', 'SIX1', 'PIK3CB', 'CDH10', 'MDC1', 'COL3A1', 'PIK3CA', 'EP300', 'CAMTA1', 'DMD', 'HLA-A' ,'HDAC6', 'CREBBP')
rest.order <- rownames(gene.sample.matrix.var)[which(!rownames(gene.sample.matrix.var) %in% gene.order)]
sample.order <- rownames(ordering.matrix)[order(
	ordering.matrix[,'TP53'], 
	-ordering.matrix[,'tmb'],
	ordering.matrix[,'FAM135B'],
	ordering.matrix[,'CSMD3'],
	ordering.matrix[,'MUC4'], 
	ordering.matrix[,'TRAF3'],
	ordering.matrix[,'APC'],
	ordering.matrix[,'PRKDC'],
	ordering.matrix[,'CTNNA2'],
	ordering.matrix[,'NOS3'],
	ordering.matrix[,'NBEA'],
	ordering.matrix[,'ATG2A'],
	ordering.matrix[,'MYO5A'],
	ordering.matrix[,'SIX1'],
	ordering.matrix[,'PIK3CB'],
	ordering.matrix[,'CDH10'],
	ordering.matrix[,'MDC1'],
	ordering.matrix[,'COL3A1'],
	ordering.matrix[,'PIK3CA'],
	ordering.matrix[,'EP300'],
	ordering.matrix[,'CAMTA1'],
	ordering.matrix[,'DMD'],
	ordering.matrix[,'HLA-A'],
	ordering.matrix[,'HDAC6'],
	ordering.matrix[,'CREBBP']
	)]

# Now plot
annotation.colours <- data.frame(
	colours = c('plum', 'lightskyblue', '#f768a1', '#e6e600', '#1e90ff', '#f6f6f6', '#767676', '#282828', '#37a124', '#e5c064'),
	labels = c('Female', 'Male', 'PR','SD', 'PD', 'Non-smoker', 'Ex-smoker', 'Smoker', 'Less than 30,000', 'More than 30,000')
	)
categorical.annotation <- categorical.annotation[sample.order,]
categorical.annotation <- apply(categorical.annotation, MARGIN=2, as.numeric)
covariate.heatmap <- create.heatmap(
	categorical.annotation,
	clustering.method = 'none',
	colour.scheme = annotation.colours$colours,
	yaxis.lab = c('Response', 'EBV', 'Smoker Status', 'Sex'),
	yaxis.cex = 0.7,
	fill.colour = 'transparent',
	total.colours = length(annotation.colours$colours) + 1,
	print.colour.key = FALSE
	)
all.legend <- legend.grob(
	list( 
		legend = list(
			colours = c('#37a124', '#e5c064'),
			labels = c( '<7,800 IU/mL', '>7,800 IU/mL'),
			title = 'EBV'
			),
		legend = list(
			colours = c('#f768a1', '#e6e600', '#1e90ff'),
			labels = c( 'PR','SD', 'PD'),
			title = 'Response'
			),
		legend = list(
			colours = c('#f6f6f6', '#767676', '#282828'),
			labels = c('Non-smoker', 'Ex-smoker', 'Smoker'),
			title = 'Smoker Status'
			),
		legend = list(
			colours = c('plum', 'lightskyblue'),
			labels = c('Female', 'Male'),
			title = 'Sex'
			),
		legend = list(
			colours = mutation.colours.labels$colours,
			labels = mutation.colours.labels$labels,
			title = 'Variant Class'
			)
		),
	label.cex = 0.6,
	title.cex = 0.7,
	size = 1,
	title.just = 'left'
	)

barplot.data <- data.frame(
	index = 1:length(sample.order),
	tmb = tmb.data[sample.order],
	row.names = sample.order
	)
tmb.barplot <- create.barplot(
	tmb ~ index,
	data = barplot.data,
	ylab.label = ' ',
	ylab.cex = 1,
	xlab.label = '',
	xat = -1,
	yat = seq(0,12,3),
	yaxis.cex = 0.7
	)
snv.heatmap <- create.heatmap(
	t(gene.sample.matrix.var[rev(gene.order),sample.order]),
	clustering.method = 'none',
	xlab.label = '', 
	ylab.label = '',
	xlab.cex = 1,
	ylab.cex = 1,
	#xaxis.lab = NA,
	yaxis.lab = NA,
	#yat = 1:50,
	xaxis.lab = gsub('NCC', '', sample.order),
	xaxis.cex = 0.7,
	yaxis.cex = 0.7,
	#colour.scheme = mutation.colours.labels$colours[sort(unique(unlist(gene.sample.matrix.var)))],
	colour.scheme = mutation.colours.labels$colours,
	fill.colour = 'transparent',
	total.colours = length(sort(unique(unlist(gene.sample.matrix.var))))+1,
	print.colour.key = FALSE
	)
create.multipanelplot(
	plot.objects = list(tmb.barplot, covariate.heatmap, snv.heatmap),
	plot.objects.heights = c(3,2,7),
	y.spacing = c(-12, -6.5),
 	layout.height = 3,
	legend = list(right = list(fun = all.legend)), 
	right.legend.padding = 0,
	resolution = 300,
	width = 5.7,
	height = 5.4,
	filename = paste0(Sys.Date(), '-IpiNivo-WES-consensus.mutations.heatmap.withtmb-medium.tiff')
	)

# And signatures 
consensus.sig.weights <- plot.consensus.sigs(consensus.mafs)
signature.order <- colnames(consensus.sig.weights)[which(colSums(consensus.sig.weights) != 0)]

colour.function <- function(x){
    colours <- rep('#bf360c', length(x));
    return(colours);
    }
spot.size.function <- function(x){1.5*abs(x)}
create.dotmap(
	x = t(consensus.sig.weights[sample.order,signature.order]),
	xlab.label = '',
	ylab.label = 'Signature weights',
	ylab.cex = 1,
	xlab.cex = 1,
	na.spot.size = 2,
	xaxis.cex = 0.7,
	yaxis.cex = 0.8,
	xaxis.rot = 90,
	xaxis.tck = 0.5,
	yaxis.tck = 0.5,
	spot.colour.function = colour.function,
	spot.size.function = spot.size.function,
	key = list(
		space = 'right',
		points = list(
			cex = spot.size.function(rev(c(0,0.25,0.5,0.75,1))),
			col = c(rep('#bf360c', 5)),
			pch = 19
			),
		text = list(
			lab = as.character(rev(c(0,0.25, 0.5, 0.75, 1))),
			cex = 0.8,
			adj = 1
			),
		padding.text = 3,
		background = 'white'
		),
	resolution = 300,
	width = 5,
	height = 3.4,
	filename = paste0(Sys.Date(), '-IpiNivo-WES-consensus.signatures.tiff')
	)

### Do some association tests
# APOBEC 
apobec.enrichment <- read.table(
	"./Projects/IpiNivo/analysis/MutSigs/2022-02-23-IpiNivo.consensus.APOBEC.trinucleotidematrix_APOBEC_enrichment.tsv",
	header = TRUE,
	sep = '\t',
	quote = ''
	)
rownames(apobec.enrichment) <- gsub('.*_', '', apobec.enrichment$Tumor_Sample_Barcode)
apobec.enrichment['NCC023',] <- NA
apobec.enrichment['NCC023', 'Tumor_Sample_Barcode'] <- 'DNA_NCC023'
apobec.enrichment <- apobec.enrichment[sample.order,]
model.data <- data.frame(
	apobec = apobec.enrichment[sample.order, 'APOBEC_Enriched'],
	sex = annotation.data[sample.order, 'Sex'],
	age = annotation.data[sample.order, 'Age'],
	response = annotation.data[sample.order, 'Site_Response'],
	is.pr = annotation.data[sample.order, 'Site_Response'] == 'PR',
	ebv = annotation.data[sample.order, 'EBV'],
	ebv.numeric = annotation.data[sample.order, 'Baseline_EBV_DNA_SI_'],
	smoking = annotation.data[sample.order, 'Smoking_Status'],
	tmb = tmb.data[sample.order],
	row.names = sample.order
	)
model.data['NCC023','apobec'] <- 'no'

apobec.response <- fisher.test(
	model.data$apobec,
	model.data$response
	)
apobec.ispr <- fisher.test(
	model.data$apobec,
	model.data$is.pr
	)
sex.response <- fisher.test(
	model.data$sex,
	model.data$response
	)
sex.ispr <- fisher.test(
	model.data$sex,
	model.data$is.pr
	)
smoking.response <- fisher.test(
	model.data$smoking,
	model.data$response
	)
smoking.ispr <- fisher.test(
	model.data$smoking,
	model.data$is.pr
	)
tmb.response <- fisher.test(
	model.data$tmb,
	model.data$response
	)
tmb.ispr <- wilcox.test(
	model.data$tmb[which(model.data$is.pr)],
	model.data$tmb[which(!model.data$is.pr)]
	)
# Some plots for supplementary 
model.data$tmb <- log10(model.data$tmb)
create.boxplot(
	tmb ~ response,
	data = model.data,
	xlab.label = 'Site-specific response',
	ylab.label = 'TMB (mut/Mb)',
	add.stripplot = TRUE,
	filename = 'tmb.vs.response.tiff',
	resolution = 200,
	ylimits = c(-2.05, 2.05),
	yat = c(-2:2),
	yaxis.lab = c(0.01, 0.1, 1, 10, 100),
	height = 3,
	width = 4,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	xlab.cex = 0.8,
	ylab.cex = 0.8
	)

barplot.data <- data.frame(
	proportion = c(0.375, 0.625, 0.5, 0.5, 0, 1),
	apobec = c('yes', 'no', 'yes', 'no', 'yes', 'no'),
	response = factor(c('PR', 'PR', 'SD', 'SD', 'PD', 'PD'), levels = c('PR', 'SD', 'PD'))
	)
bar.col = c('grey30', 'white')
create.barplot(
	proportion ~ response, 
	data = barplot.data,
	groups = barplot.data$apobec,
	col = bar.col,
	stack = TRUE,
	filename = 'apobec.vs.response.tiff',
	resolution = 200,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	xlab.cex = 0.8,
	ylab.cex = 0.8,
	xlab.label = 'Site-specific response',
	ylab.label = 'Proportion APOBEC-enriched',
	xaxis.lab = c('PR (n=8)', 'SD (n=6)', 'PD (n=6)'),
	height = 3,
	width = 3.8
	)
