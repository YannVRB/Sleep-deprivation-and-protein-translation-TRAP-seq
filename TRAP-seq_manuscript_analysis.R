######################################################################################################################################################
########## Title
########## Yann Vanrobaeys
########## August 10, 2020

######################################################################################################################################################
#################################################################################################### LOADING PACKAGES.
library(edgeR) # edgeR_3.28.1
library(EDASeq) # EDASeq_2.20.0
library(ggplot2) # ggplot2_3.3.2
library(ggrepel) # ggrepel_0.8.2
library(ggpmisc) #ggpmisc_0.3.5

######################################################################################################################################################
#################################################################################################### Upper quartile normalization for sequencing depth.
# setting up groups based on the two conditions (5 non-sleep deprived VS 6 sleep deprived).
group <- factor(c(1,1,1,1,1,2,2,2,2,2,2))
# create DGEList for edgeR.
y <- DGEList(counts=raw_counts_matrix,group=group)

# filter undetected / lowly expressed genes.
keep <- filterByExpr(y) 
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
# create design matrix.
design <- model.matrix(~group)
# edgeR glm workflow.
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
topRsFC <- topTags(lrt, n=Inf)$table
write.table(topRsFC, file=" DGEanalysis_UQ.csv", sep = ",")

######################################################################################################################################################
##################################################################################################### Quadrant diagrams 
# Significant DEGs in common
colnames(`49DEGsoverlap_FC`) <- c("Vecsey", "TRAP")
ggplot(`49DEGsoverlap_FC`, aes(x= Vecsey, y= TRAP)) + 
  geom_point() +
  geom_text_repel(label=rownames(`49DEGsoverlap_FC`)) +
  xlab("Fold Change from Vecsey paper") +
  ylab("Fold Change from TRAP-seq") +
  ggtitle("Comparison of 49 DEGs overlapping Vecsey and TRAP data") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)

# significant DEGs from Vecsey et al. (2019) with fold change from TRAP (this study). Genes are labeled if they have an absolute delta fold change above 0.4.
colnames(X280VecseyDEGS_withTRAPFC_absDeltaFC) <- c("gene", "Vecsey", "TRAP", "Delta", "AbsDelta")
dat2 <- subset(X280VecseyDEGS_withTRAPFC_absDeltaFC, AbsDelta > 0.4)
ggplot(X280VecseyDEGS_withTRAPFC_absDeltaFC, aes(x= Vecsey, y= TRAP, label=gene)) + 
  geom_point() +
  geom_text_repel(data = dat2, force = 9) +
  xlab("Fold Change from Vecsey paper") +
  ylab("Fold Change from TRAP-seq") +
  ggtitle("Comparison of 280 DEGs from Vecsey with respective FC TRAP") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  scale_x_continuous(limit = c(-2, 2)) +
  scale_y_continuous(limit = c(-2, 2))

# Significant DEGs from TRAP (this study) with fold change from Vecsey et al. (2019). Genes are labeled if they have an absolute delta fold change above 0.4.
colnames(X149TRAPDEGs_withVecseyFC_absDeltaFC) <- c("gene", "TRAP", "Vecsey", "Delta", "AbsDelta")
dat2 <- subset(X149TRAPDEGs_withVecseyFC_absDeltaFC, AbsDelta > 0.4)
ggplot(X149TRAPDEGs_withVecseyFC_absDeltaFC, aes(x= Vecsey, y= TRAP, label=gene)) + 
  geom_point() +
  geom_text_repel(data = dat2, force = 9) +
  xlab("Fold Change from Vecsey paper") +
  ylab("Fold Change from TRAP-seq") +
  ggtitle("Comparison of 149 DEGs from TRAP with respective FC Vecsey") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  scale_x_continuous(limit = c(-2, 2)) +
  scale_y_continuous(limit = c(-2, 2))
