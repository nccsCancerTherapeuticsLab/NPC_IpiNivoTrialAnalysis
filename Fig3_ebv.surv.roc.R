### ebv.surv.rocs.R  #############################################################################
#
# Connie Li, 	November 2020 

### DESCRIPTION ###################################################################################
#
# ROC analysis of EBV thresholds 
#
### LIBRARIES, ENV & INPUT ########################################################################
setwd("./Projects/IpiNivo/RSession")
library(survival)
library(survminer)
library(pROC)

clinical.data <- read.table(
	'./Projects/IpiNivo/data/annotation/SampleInformation.txt',
	sep ='\t',
	header = TRUE,
	quote = ''
	)
rownames(clinical.data) <- clinical.data$Subject_no
pfs.data <- read.table(
	'./Projects/IpiNivo/data/annotation/NPCtrial _PFS 20211001.csv',
	sep = ',',
	header = TRUE,
	quote = ''
	)
rownames(pfs.data) <- pfs.data$Patient.ID
rownames(pfs.data) <- gsub('^1', 'NCC', rownames(pfs.data))
rownames(pfs.data) <- gsub('^2', 'NUH', rownames(pfs.data))
rownames(pfs.data) <- gsub('^3', 'TW', rownames(pfs.data))

overlapping.pts <- intersect(rownames(pfs.data)[which(pfs.data$n.26.Phase.II.set == 1)], rownames(clinical.data))
surv.modeling.data <- cbind(
	clinical.data[overlapping.pts, c('Baseline_EBV_DNA_copies.ml_', 'Baseline_EBV_DNA_SI_', 'Current_status', 'Best_overall_response')],
	pfs.data[overlapping.pts, c('PFS.event.status', 'PFS.time.months')],
	clinical.data[overlapping.pts, c('Age', 'Sex', 'Ethnic_group', 'Smoking_Status', 'Staging_WHO_', 'Staging_AJCC_')]
	)
surv.modeling.data$Baseline_EBV_DNA_copies.ml_<- as.numeric(gsub('>|<', '', surv.modeling.data$Baseline_EBV_DNA_copies.ml_))
surv.modeling.data$Baseline_EBV_DNA_SI_ <- as.numeric(gsub('>|<', '', surv.modeling.data$Baseline_EBV_DNA_SI_))

# EBV cutoff is 7800
surv.modeling.data$ebv.high <- surv.modeling.data$Baseline_EBV_DNA_SI_ > 7800

pfs.survobj <- Surv(
	time = surv.modeling.data$PFS.time.months,
	event = surv.modeling.data$PFS.event.status
	)
pfs.ebv.survfit <- survfit(pfs.survobj ~ surv.modeling.data$ebv.high)
pfs.ebv.survfit <- survfit(
	Surv(PFS.time.months, PFS.event.status) ~ ebv.high,
	data = surv.modeling.data
	)
png(
	filename = paste0(Sys.Date(), '-pfs.ebv.kmplot.png'),
	width = 5.5,
	height = 5,
	units = 'in',
	res = 500
	)
ggsurvplot(
	fit =  pfs.ebv.survfit,
	palette = c('#37a124', '#e5c064'),
	risk.table = TRUE, 
	legend.labs = c(expression('EBV<7,800 iu/mL'), expression('EBV\u22657,800 iu/mL')),
	xlab = "Months", 
	ylab = "Progression-free survival"
	)
dev.off()
surv.modeling.data$pr.ornot <- surv.modeling.data$Best_overall_response == 'PR'
surv.modeling.data$pr.pd <- surv.modeling.data$Best_overall_response == 'PR'
surv.modeling.data$pr.pd[which(surv.modeling.data$Best_overall_response == 'SD')] <- NA
surv.modeling.data$pr.pd[which(surv.modeling.data$Best_overall_response == 'SD')] <- NA

# Then ROC curves for response 
png(
	filename = paste0(Sys.Date(), '-response.pr.ornot.ebv.roc.png'),
	width = 4,
	height = 4,
	units = 'in',
	res = 300
	)
plot.roc(
	surv.modeling.data$pr.ornot,
	surv.modeling.data$Baseline_EBV_DNA_SI_,
	print.auc = TRUE,
	print.thres = 'best',
	lwd=2,
	print.auc.y = 30,
	print.auc.x = 80,
	ci = TRUE,
	percent = TRUE
	)
dev.off()
# Get optimal cut point 
ebv.prornot.roc <- plot.roc(
	surv.modeling.data$pr.ornot,
	surv.modeling.data$Baseline_EBV_DNA_SI_,
	print.auc = TRUE,
	lwd=2,
	print.auc.y = 30,
	print.auc.x = 80,
	ci = TRUE,
	percent = TRUE
	)
optimal.coord <- coords(
	ebv.prornot.roc, 
	x = 'best'
	)

# Look at pr pd 
ebv.prpd.roc <- plot.roc(
	surv.modeling.data[which(!is.na(surv.modeling.data$pr.pd)),]$pr.pd,
	surv.modeling.data[which(!is.na(surv.modeling.data$pr.pd)),]$Baseline_EBV_DNA_SI_,
	add = TRUE,
	lty = 2,
	col = 'grey30',
	print.auc = TRUE,
	lwd=2,
	print.auc.y = 30,
	print.auc.x = 80,
	ci = TRUE,
	percent = TRUE
	)
prpd.optimal.coord <- coords(
	ebv.prpd.roc, 
	x = 'best'
	)

png(
	filename = paste0(Sys.Date(), '-response.pr.or.pd.ebv.roc.png'),
	width = 4,
	height = 4,
	units = 'in',
	res = 300
	)
plot.roc(
	surv.modeling.data$pr.pd,
	surv.modeling.data$Baseline_EBV_DNA_SI_,
	#add = TRUE,
	lty = 2,
	col = 'grey30',
	print.auc = TRUE,
	print.thres = 'best',
	lwd=2,
	print.auc.y = 30,
	print.auc.x = 80,
	ci = TRUE,
	percent = TRUE
	)
dev.off()

# Print combined ROC
png(
	filename = paste0(Sys.Date(), '-response.both.ebv.roc.png'),
	width = 4.2,
	height = 4.2,
	units = 'in',
	res = 300
	)
plot.roc(
	surv.modeling.data$pr.ornot,
	surv.modeling.data$Baseline_EBV_DNA_SI_,
	#print.auc = TRUE,
	print.thres = 'best',
	lwd=2,
	print.auc.y = 30,
	print.auc.x = 60,
	ci = TRUE,
	percent = TRUE
	)
plot.roc(
	surv.modeling.data[which(!is.na(surv.modeling.data$pr.pd)),]$pr.pd,
	surv.modeling.data[which(!is.na(surv.modeling.data$pr.pd)),]$Baseline_EBV_DNA_SI_,
	add = TRUE,
	lty = 2,
	col = 'grey30',
	lwd=2,
	print.auc.y = 30,
	print.auc.x = 60,
	ci = TRUE,
	percent = TRUE
	)
dev.off()

# View new thresholds as picked by ROC
surv.modeling.data <- cbind(
	clinical.data[overlapping.pts, c('Baseline_EBV_DNA_SI_', 'Current_status', 'Best_overall_response')],
	pfs.data[overlapping.pts, c('PFS.event.status', 'PFS.time.months')],
	clinical.data[overlapping.pts, c('Age', 'Sex', 'Ethnic_group', 'Smoking_Status', 'Staging_WHO_', 'Staging_AJCC_')]
	)
surv.modeling.data$Baseline_EBV_DNA_SI_ <- as.numeric(gsub('>|<', '', surv.modeling.data$Baseline_EBV_DNA_SI_))
surv.modeling.data$ebv.high <- surv.modeling.data$Baseline_EBV_DNA_SI_ > 7868
surv.modeling.data$ebv.high <- surv.modeling.data$Baseline_EBV_DNA_SI_ > 20680

pfs.survobj <- Surv(
	time = surv.modeling.data$PFS.time.months,
	event = surv.modeling.data$PFS.event.status
	)
pfs.ebv.survfit <- survfit(pfs.survobj ~ surv.modeling.data$ebv.high)
pfs.ebv.survfit <- survfit(
	Surv(PFS.time.months, PFS.event.status) ~ ebv.high,
	data = surv.modeling.data
	)
png(
	filename = paste0(Sys.Date(), '-pfs.ebv.7868.kmplot.png'),
	width = 5.5,
	height = 5,
	units = 'in',
	res = 500
	)
ggsurvplot(
	fit =  pfs.ebv.survfit,
	palette = c('#37a124', '#e5c064'),
	risk.table = TRUE, 
	legend.labs = c(expression('EBV<7,800 IU/mL'), expression('EBV\u22657,800 IU/mL')),
	xlab = "Months", 
	ylab = "Progression-free survival"
	)
dev.off()
surv_pvalue(pfs.ebv.survfit)
surv.modeling.data$pr.ornot <- surv.modeling.data$Best_overall_response == 'PR'
surv.modeling.data$pr.pd <- surv.modeling.data$Best_overall_response == 'PR'
surv.modeling.data$pr.pd[which(surv.modeling.data$Best_overall_response == 'SD')] <- NA

ebv.threshold.data <- cbind(
	clinical.data[overlapping.pts, c('Baseline_EBV_DNA_copies.ml_', 'Baseline_EBV_DNA_SI_', 'Best_overall_response')],
	pfs.data[overlapping.pts, c('PFS.event.status', 'PFS.time.months')]
	)
ebv.threshold.data$Baseline_EBV_DNA_copies.ml_ <- as.numeric(gsub('>|<', '', ebv.threshold.data$Baseline_EBV_DNA_copies.ml_))
ebv.threshold.data$ebv.30000 <- ebv.threshold.data$Baseline_EBV_DNA_copies.ml > 30000
ebv.threshold.data$Baseline_EBV_DNA_SI_ <- as.numeric(gsub('>|<', '', ebv.threshold.data$Baseline_EBV_DNA_SI_))
ebv.threshold.data$ebv.3500 <- ebv.threshold.data$Baseline_EBV_DNA_SI_ > 3500
ebv.threshold.data$ebv.7868 <- ebv.threshold.data$Baseline_EBV_DNA_SI_ > 7868
ebv.threshold.data$ebv.20680 <- ebv.threshold.data$Baseline_EBV_DNA_SI_ > 20680
ebv.threshold.data <- ebv.threshold.data[order(ebv.threshold.data$Baseline_EBV_DNA_SI_),]

# Extra supp plot
# For EBV, we want to look at all pts
model.data <- data.frame(
	response = surv.modeling.data[, 'Best_overall_response'],
	ebv = surv.modeling.data[,'Baseline_EBV_DNA_SI_'] > 7800,
	ebv.numeric = surv.modeling.data[, 'Baseline_EBV_DNA_SI_']
	)
model.data$is.pr <- model.data$response == 'PR'
model.data$response <- factor(model.data$response, levels = c('PR', 'SD', 'PD'))
points.col <- rep('#37a124', nrow(model.data))
points.col[which(model.data$ebv == TRUE)] <- '#e5c064'
create.boxplot(
	ebv.numeric ~ response,
	data = model.data,
	xlab.label = 'Site-specific response',
	ylab.label = 'EBV (IU/mL)',
	add.stripplot = TRUE,
	points.col = points.col, 
	filename = 'ebv.vs.response.tiff',
	resolution = 200,
	height = 3,
	width = 4,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	xlab.cex = 0.8,
	ylab.cex = 0.8
	)
ebv.response <- fisher.test(
	model.data$ebv,
	model.data$response
	)
ebv.ispr <- fisher.test(
	model.data$ebv,
	model.data$is.pr
	)
ebv.utest.ispr <- wilcox.test(
	model.data$ebv.numeric[which(model.data$is.pr)],
	model.data$ebv.numeric[which(!model.data$is.pr)]
	)
ebv.kw <- kruskal.test(
	ebv.numeric ~ response,
	data = model.data
	)