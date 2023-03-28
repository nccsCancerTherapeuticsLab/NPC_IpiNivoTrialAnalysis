compare.tmb.prevstudy <- function(tmb.data) {
	# Read in TMB from other studies 
	chhk.maf <- read.table(
		'TMBOtherStudy/Li-Lo_CHHK_NatComms.txt',
		sep = '\t',
		header = TRUE
		)
	# They formatted the table differently, so let's fill in the sample IDs 
	for (i in 2:nrow(chhk.maf)) {
		if (chhk.maf$Tumor.ID..N.111.[i] == '') {
			chhk.maf$Tumor.ID..N.111.[i] <- chhk.maf$Tumor.ID..N.111.[i-1]
		}
	}
	chhk.tmb <- as.numeric(summary(factor(chhk.maf$Tumor.ID..N.111.), maxsum = 200))
	names(chhk.tmb) <- unique(chhk.maf$Tumor.ID..N.111.)
	chhk.tmb <- chhk.tmb/35.1 #Mb 

	nus.tmb <- read.table(
		'Lin-Koeffler_NUS_NatGen.txt',
		sep = '\t', 
		header = TRUE
		)[1:56,]
	nus.tmb$tmb <- (as.numeric(nus.tmb$Missense) + as.numeric(nus.tmb$Stopgain..Framefshift..Splicing))/35.1 

	hunan.tmb <- read.table(
		'Tu-Xiong_Hunan_Carcinogenesis.txt',
		sep = '\t',
		header = TRUE
		)
	# Combine and plot 
	all.tmb.data <- rbind(
		cbind(tmb.data, rep('NCCS', length(tmb.data))),
		cbind(hunan.tmb$Mutation.Mb, rep('Tu et al', nrow(hunan.tmb))),
		cbind(nus.tmb$tmb, rep('Lin et al', nrow(nus.tmb))),
		cbind(chhk.tmb, rep('Li et al', length(chhk.tmb)))
		)
	colnames(all.tmb.data) <- c('tmb', 'study')
	all.tmb.data <- as.data.frame(all.tmb.data)
	all.tmb.data$tmb <- as.numeric(all.tmb.data$tmb)
	all.tmb.data$log.tmb <- log10(all.tmb.data$tmb)
	all.tmb.data$study <- factor(all.tmb.data$study, levels=c('NCCS', 'Li et al', 'Lin et al', 'Tu et al'))
	create.boxplot(
		log.tmb ~ study,
		data = all.tmb.data,
		xlab.label = '',
		ylab.label = 'TMB (mut/Mb)',
		xaxis.lab = c('Current\nstudy', 'Li et al', 'Lin et al', 'Tu et al'),
		ylimits = c(-2.05, 2.05),
		yat = c(-2:2),
		yaxis.lab = c(0.01, 0.1, 1, 10, 100),
		add.stripplot = TRUE,
		filename = 'compare.study.tmb.tiff',
		resolution = 200,
		height = 3,
		width = 4,
		xaxis.cex = 0.8,
		yaxis.cex = 0.8,
		xlab.cex = 0.8,
		ylab.cex = 0.8
		)
	# Do and return utests 
	utest.p <- rep(NA, 3)
	names(utest.p) <- c('Li', 'Lin', 'Tu')
	utest.p['Li'] <- wilcox.test(
		all.tmb.data$tmb[which(all.tmb.data$study == 'NCCS')],
		all.tmb.data$tmb[which(all.tmb.data$study == 'Li et al')],
		)$p.value
	utest.p['Lin'] <- wilcox.test(
		all.tmb.data$tmb[which(all.tmb.data$study == 'NCCS')],
		all.tmb.data$tmb[which(all.tmb.data$study == 'Lin et al')],
		)$p.value
	utest.p['Tu'] <- wilcox.test(
		all.tmb.data$tmb[which(all.tmb.data$study == 'NCCS')],
		all.tmb.data$tmb[which(all.tmb.data$study == 'Tu et al')],
		)$p.value
	return(utest.p)
}