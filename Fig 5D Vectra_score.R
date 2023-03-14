library(tidyverse)
library(broom)
library(readxl)

#Vectra Score is the Proportion of cells as counted by multiplex imaging
vectra_score_input <- read_excel("vectra_score_input.xlsx")
vectra_pdpr <- filter(vectra_score_input, Response != "SD")
vectra_pdpr$Response <- sub("PR", "1", vectra_pdpr$Response)
vectra_pdpr$Response <- sub("PD", "0", vectra_pdpr$Response)

#Creating a score/model
fit <- lm(Response~`PD1+CTLA4-CD8/CD8` + `PD1-CTLA4+CD8/CD8`, data = vectra_pdpr)
summary(fit)
a<-tidy(fit)

library(glmnet)
cv.fit1 <- cv.glmnet(x= data.matrix(vectra_pdpr[,c(4,5)]), y= as.factor(vectra_pdpr$Response), family = "binomial", alpha=0.5) #elastic net
plot(cv.fit1)
tmp_coeffs1 <- coef(cv.fit1, s = "lambda.min") #or lambda.1se
tmp_coeffs1

c <- tmp_coeffs1@x[1]
pd1 <- tmp_coeffs1@x[2]
ctla4 <- tmp_coeffs1@x[3]

vectra_score_input <- vectra_score_input %>% mutate(Score = c + pd1*`PD1+CTLA4-CD8/CD8`+ ctla4*`PD1-CTLA4+CD8/CD8`)
vectra_score_input$Response <- factor(vectra_score_input$Response, levels = c("PR", "SD", "PD"))

#Function to plot
fil.col <- c('#f768a1','#e6e600', '#1e90ff')
ggplotRegression <- function (db,dep_val,gene) {
  require(ggplot2)
  db <- na.omit(db) #REMOVE THIS FOR PFS
  eqn <- formula(paste0(dep_val," ~ ",gene))
  t <- cor(db[,dep_val], db[,gene], method = "pearson") #calculate r value
  t<- round(t, digits = 3) #3 decimal place
  fit <- lm(eqn, data = db)
  ggplot(db, aes(x = !!ensym(dep_val), y = !!ensym(gene))) + 
    geom_point(aes(fill = Response), colour="black",pch=21, size = 5) +  guides(fill="legend") +
    stat_smooth(method = "lm", col = "red") + theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(hjust = 0.5, size = 20), aspect.ratio = 1,
          text = element_text(size = 17),
          legend.position = "none") + #remove legend
    scale_fill_manual(values= fil.col)+ labs(title = paste(gene, "     r = ", t ,",  P =",signif(summary(fit)$coef[2,4], 4))) +
    ylab("PD1 CTLA4 Score") 
}


a<- ggplotRegression(vectra_score_input, "TTP", "Score")
while (!is.null(dev.list()))  dev.off() #prevent pdf write problem
pdf("vectra_score.pdf", width = 5.5, height = 4)
a
dev.off()

