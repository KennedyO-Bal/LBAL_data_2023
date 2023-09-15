## Load libraries----
library(DESeq2)
library(ggplot2)
library(ggside)
library("pheatmap")
library(RColorBrewer)
library(tidyverse)
library(sva)
library(PCAtools)
library(viridis)
library(ggrepel)
library(EnhancedVolcano)
library(VennDiagram)
library(beepr)
library(vegan) 
library(gplots) 
library(ggplot2) 
library(pheatmap) 
library(UpSetR)  

getwd() 
setwd("~/Desktop/Paper Witing, Genomic & RNA seqs/Miseq and novaseq data_2021_NEB/Japan data_2021 collections/novaseq/Round 2_counts/subsets")


#load count files----
Counts<-read.delim("adult.brain_feature.counts new.txt", header = T, row.names = 1, stringsAsFactors = F, check.names=FALSE)
colSums(Counts)
str(Counts)

#Load unfiltered trait data
sample.info.all<- read.csv("sample.info.metadata.csv", header = TRUE )

sample.info <- sample.info.all[which(sample.info.all$sex == 'female' & 
                                       sample.info.all$caste == 'adult' & 
                                       sample.info.all$tissue == 'brain' &
                                       sample.info.all$total_gene.counts > 0),] 

# Graphics to observe before normalization

barplot(colSums(Counts)*1e-6, col="light blue",names=colnames(Counts), ylab="Library size (millions)")

## Drop samples with less than x million reads per sample
filt.counts <- Counts[ ,(colSums(Counts) >= 2*1e6) ]
sample.info.filt  <- sample.info [which(sample.info $sample.name %in% colnames(filt.counts)),]

# reorder to make sure colnames and rownames are in same order
reorder_idx <- match(sample.info.filt$sample.name, colnames(filt.counts))
filt.counts <- filt.counts[reorder_idx]
nrow(sample.info.filt) == ncol(filt.counts)

# check that sample names match
sample.info.filt$sample.name == colnames(filt.counts)
#filt.counts

#set groupings
behavior<- factor(sample.info.filt$behavior)
sociality<- factor(sample.info.filt$sociality)
sample.name<-factor(sample.info.filt$sample.name)
nest<-factor(sample.info.filt$nest)
probe<-factor(sample.info.filt$temp_probe)
Location<-factor(sample.info.filt$Location)
Group<-factor(sample.info.filt$Group)

all((sample.info.filt$sample.name) == colnames(filt.counts))



##Buid DESeq2 object----

dds <- DESeqDataSetFromMatrix(countData = filt.counts,
                              colData = sample.info.filt, 
                              design= ~behavior + nestID) #with unfiltered counts

#check samples match
colnames(dds) == sample.info.filt$sample.name
length(colnames(dds))

## Normalizing using the method for an object of class 'count data' 
dds <-  estimateSizeFactors(dds)
counts(dds, normalized=TRUE) #divides each column by its size factor.

keep <- rowSums(counts(dds) >= 2) >= 10# rows with sum of all counts bigger than 10
dim(dds)
dds<-dds[keep,]
dim(dds) #11902 genes

dds <- DESeq(dds)

ddsClean <- dds[which(mcols(dds)$betaConv),] #this removes remaining rows that that don't converge in beta
dim(ddsClean) #11854

dds<-ddsClean


dds.sva=dds #use this for svaseq in PART 2.normalized counts

# Including sva cariables

# SVAseq for latent variable characterization
dat  <- counts(dds.sva, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ behavior , colData(dds.sva))
mod0 <- model.matrix(~   1, colData(dds.sva))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

# use latent variables in downstream DE analyses
ddssva <- dds.sva
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]

# plot latent variables by temp probe
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$temp_probe, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

# plot latent variables by nest
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$nestID, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

# plot latent variables by behavior
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$behavior, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}


# does not correlate with behavior
dev.off()

plot(svseq$sv[,1:2], col=ddssva$behavior, pch=16, cex=2,
     xlab="SV1", ylab="SV2")
legend("bottomright", levels(ddssva$behavior), pch=16,
       col=1:4, cex=.8, ncol=4, title="behavior")


dds.2 <- DESeqDataSet(ddssva, design = ~ SV1 + SV2 + behavior + nestID)

# filter low count genes
dds.2 <- dds.2[ rowSums(counts(dds.2) >= 2) >= 10, ]
dim(dds.2) #11858

## Make DESeq table with svas included
dds.2 <- DESeq(dds.2)
head(results(dds.2, tidy=TRUE))
resultsNames(dds.2)

ddsClean2 <- dds.2[which(mcols(dds.2)$betaConv),] #this removes remaining rows that that don't converge in beta

dim(ddsClean2) #11854

dds.2<-ddsClean2

# Graphics to observe after normalization
dds.2_after<-counts(dds.2, normalized = TRUE)
barplot(colSums(dds.2_after), col = "light blue", ylab = "Total number of read counts")

# Estimate data reproducibility
dds.rep = estimateDispersions(dds.2)
plotDispEsts(dds.rep)

## To be used for comparing Workers to all reproductives(Queens and sol_reproductive)
### DEGs for combinatorial contrasts ###
dds.3 <- DESeqDataSet(ddssva, design = ~ 0 + behavior + SV1 + SV2 + nestID)
# filter low count genes
dds.3 <- dds.3[ rowSums(counts(dds.3) >= 2) >= 10, ]
dim(dds.3) #11854
## Make DESeq table with svas included
dds.3 <- DESeq(dds.3)
head(results(dds.3, tidy=TRUE))
resultsNames(dds.3)
ddsClean3 <- dds.3[which(mcols(dds.3)$betaConv),] #this removes remaining rows that that don’t converge in beta
dim(ddsClean3) #11854
dds.3<-ddsClean3

#PCA, Heatmatrix with SVA variables included. Same code as in part 1----

#Extract transformed values

#Heatmaps of count matrix----

vsd.sva <- vst(dds.2, blind=FALSE) 
head(assay(vsd.sva), 3)


hmcol <- viridis(150)

anot_col<- brewer.pal(12,"Paired") 

annotation_colors <- list(
  sociality = c(soc=anot_col[1], sol=anot_col[11],unknown=anot_col[10]),
  behavior=c(foundress = anot_col[12], queen = anot_col[6], sol_reprod=anot_col[2], worker =anot_col[7])
)

df.sva <- as.data.frame(colData(dds.2)[,c("behavior","sociality")])

#selection of 250  genes
select.sva <- order(rowMeans(counts(dds.2,normalized=TRUE)), 
                    decreasing=TRUE)[1:250]

#all genes.!!Running this takes a lot of time

allGenes.sva<-order(rowMeans(counts(dds.2,normalized=TRUE)),
                    decreasing=TRUE)

#vst transformed heatmap(variance stabilized gene expression values). 
pheatmap(assay(vsd.sva)[select.sva,], #or allGenes.sva
         cluster_rows=TRUE, 
         show_rownames=FALSE, 
         show_colnames = F, 
         cluster_cols=TRUE, 
         annotation=df.sva,
         annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         clustering_distance_cols = "manhattan",
         cutree_rows = 1,
         cutree_cols = 2,
         color=hmcol,
         width = 7,
         height=6)
dev.off()

#Heatmap of sample distances. Overview over similarities and dissimilarities between samplessampleDists.sva <- dist(t(assay(vsd.sva)))
sampleDists.sva <- dist(t(assay(vsd.sva)))
sampleDistMatrix.sva <- as.matrix(sampleDists.sva)
rownames(sampleDistMatrix.sva) <- paste(vsd.sva$behavior, sep="-")
colnames(sampleDistMatrix.sva) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) ) (54)  #54 samples
pheatmap(sampleDistMatrix.sva,
         clustering_distance_rows=sampleDists.sva,
         clustering_distance_cols=sampleDists.sva,
         col=colors)

pcaData<-plotPCA(vsd.sva, intgroup=c( "behavior"),returnData = T)# returnData=T returns the data frame useful for plotting in ggplot


percentVar <- round(100 * attr(pcaData, "percentVar"))# percent variance explained. 


###PCA type 1----all genes
pca1<-ggplot(pcaData,aes(PC1,PC2,shape= behavior, color= behavior))+
  geom_point(size=3.5,color="black")+
  geom_point(size=3,aes(color=behavior))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  #stat_ellipse() +
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_hline(yintercept = 0, linetype="dotted")+
  # scale_color_manual(values =anot_col[c(1,11)])+
  theme_bw() +
  scale_x_continuous(limits=c(-40, 60)) +
  scale_y_continuous(limits=c(-40, 40)) 


pca1

gg<-pca1 + scale_color_manual(values=c("#f1d37c", "#511314", "#345c67", "#9e2b2b")) + theme_bw() + scale_shape_manual(values=c(17,16,15,18))  + theme_bw() + theme(legend.title=element_blank())

gg

#ggsave("adult_brains.pdf",gg, width=3, height=3, units="in", scale=3)


#These plots blank themes. Great for poster ppt.
pca1  +    
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )



##########################################
#Regression plots

lm(pcaData$PC1 ~ vsd.sva$longest.oocyte, data = pcaData) |> summary()
plot(-pcaData$PC1 ~ vsd.sva$longest.oocyte)
abline(lm(-pcaData$PC1 ~ vsd.sva$longest.oocyte))

#OR 

fit <- lm(-pcaData$PC1 ~ vsd.sva$longest.oocyte, data = pcaData)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  theme_bw() +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))

#Anova
#PC1 and behavior association

brain_anov<-aov(pcaData$PC1 ~ behavior,
                data=pcaData)
summary(brain_anov)

########################################

# exporting the vst values for the WGCNA analysis. Normalized transformed
dds_vst.df.sva <- as.data.frame(assay(vsd.sva))
dds_vst_transposed.sva <- t(dds_vst.df.sva)
#write.csv(dds_vst_transposed.sva, file = "vst_counts_core_for WGCNA_adult.brains_sva.csv") #can be saved as untransposed too

######PCAs
#PCA 1
z.sva <- plotPCA(vsd.sva, intgroup = "behavior" )
z.sva +theme_bw() + geom_text_repel(aes(label = name)) 

z.sva + geom_label(aes(label=row.names(sample.info.filt))) + theme_bw()
z.sva$data

z.df.sva<-as.tibble(z.sva$data)


#PCA 2. Detailed PCA
pcaData.sva <- plotPCA(vsd.sva, intgroup="behavior", returnData=TRUE)
percentVar.sva <- round(100 * attr(pcaData.sva, "percentVar"))
pca2 <- ggplot(pcaData.sva, aes(PC1, PC2, shape=behavior, group=behavior, color= behavior)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.sva[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.sva[2],"% variance")) + 
  #geom_text(aes(label=rownames(sample.info.filt)),hjust=.5, vjust=-.8, size=3) +
  # geom_density2d(alpha=.5) +
  stat_ellipse() +
  theme_bw() +
  coord_fixed()
# pca2 + scale_color_manual(values=c("orange2", "lightsalmon4", "darkslategray4", "tomato3")) + theme_bw()
pca2 + scale_shape_manual(values=c(17,16,15,2,1,0))  + theme_bw() + theme(legend.title=element_blank())



# Detailed PCA of small multiple plots

new.z.sva.df<- z.sva$data [,1:2] %>% # PC's 1$2
  as_tibble() %>%
  add_column(behavior) %>%
  add_column(sample.name) %>%
  add_column(sociality)


z.pivot<- pivot_longer (new.z.sva.df,
                        cols = PC1:PC2,
                        names_to = "PC",
                        values_to = "loadings")
p3<-ggplot(z.pivot) +
  aes(x=behavior, y=loadings, fill=sociality) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC) +
  labs(title = "PCA of 'small multiples' plot",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
#beep(sound = 6)

#PCA type 3_ggside
ggplot(z.df.sva, aes(PC1, PC2, color=behavior)) + 
  geom_point(size = 2) +
  geom_xsideboxplot(aes(y =behavior), orientation = "y") +
  scale_xsidey_discrete() + #In order to use xsideboxplot with a main panel that uses
  #continuous data, force y axis in xside panel to be discrete
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL) +
  theme(ggside.panel.scale = .5) +
  theme_bw() 

pca3

#Pairwise comparisons_ with SVAs----

## comparing Queens and Workers
Q_W_results.sva <- results(dds.2, contrast = c("behavior","queen","worker"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(Q_W_results.sva, na.rm=TRUE)
Q_W_results.sva <- as.data.frame(Q_W_results.sva)
Q_W_results.sva_sig <- subset(Q_W_results.sva, padj < 0.001)
Q_W_results.sva_sig <- as.data.frame(Q_W_results.sva_sig)
length(Q_W_results.sva_sig$padj)#166

## comparing Queens and Solitary reproductives
Q_S_results.sva <- results(dds.2, contrast = c("behavior","queen","sol_reprod"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(Q_S_results.sva, na.rm=TRUE)
Q_S_results.sva <- as.data.frame(Q_S_results.sva)
Q_S_results.sva_sig <- subset(Q_S_results.sva, padj < 0.001)
Q_S_results.sva_sig <- as.data.frame(Q_S_results.sva_sig)
length(Q_S_results.sva_sig$padj)#0

## comparing workers and solitary reproductive
W_S_results.sva <- results(dds.2, contrast = c("behavior","worker","sol_reprod"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(W_S_results.sva, na.rm=TRUE)
W_S_results.sva <- as.data.frame(W_S_results.sva)
W_S_results.sva_sig <- subset(W_S_results.sva, padj < 0.001)
W_S_results.sva_sig <- as.data.frame(W_S_results.sva_sig)
length(W_S_results.sva_sig$padj)#52

## comparing Queens and Foundress
Q_F_results.sva <- results(dds.2, contrast = c("behavior","queen","foundress"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(Q_F_results.sva, na.rm=TRUE)
Q_F_results.sva <- as.data.frame(Q_F_results.sva)
Q_F_results.sva_sig <- subset(Q_F_results.sva, padj < 0.001)
Q_F_results.sva_sig <- as.data.frame(Q_F_results.sva_sig)
length(Q_F_results.sva_sig$padj)#8

## comparing Foundress and Solitary reproductives
F_S_results.sva <- results(dds.2, contrast = c("behavior","foundress","sol_reprod"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(F_S_results.sva, na.rm=TRUE)
F_S_results.sva <- as.data.frame(F_S_results.sva)
F_S_results.sva_sig <- subset(F_S_results.sva, padj < 0.001)
F_S_results.sva_sig <- as.data.fram0e(F_S_results.sva_sig)
length(F_S_results.sva_sig$padj)#1

## comparing Foundress and Workers
F_W_results.sva <- results(dds.2, contrast = c("behavior","foundress","worker"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(F_W_results.sva, na.rm=TRUE)
F_W_results.sva <- as.data.frame(F_W_results.sva)
F_W_results.sva_sig <- subset(F_W_results.sva, padj < 0.001)
F_W_results.sva_sig <- as.data.frame(F_W_results.sva_sig)
length(F_W_results.sva_sig$padj)#265

## saving the DEG files for each comparison

setwd()

# write.csv(Q_W_results.sva_sig, file = "Q_W_DEGs_adult_brains_sva.csv")
#write.csv(Q_S_results.sva_sig, file = "Q_S_DEGs_adult_brains_sva.csv")
#write.csv(W_S_results.sva_sig, file = "W_S_DEGs_adult_brains_sva.csv")
#write.csv(Q_F_results.sva_sig, file = "Q_F_DEGs_adult_brains_sva.csv")
#write.csv(F_S_results.sva_sig, file = "F_S_DEGs_adult_brains_sva.csv")
#write.csv(F_W_results.sva_sig, file = "F_W_DEGs_adult_brains_sva.csv")



### comparing Workers to all reprods
### DEGs for combinatorial contrasts ###
W_vs_QS_results.sva <- results(dds.3,
                               contrast = list( c("behaviorworker"),c("behaviorqueen", "behaviorsol_reprod") ),
                               listValues = c(1,-1/2),
                               alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(W_vs_QS_results.sva, na.rm=TRUE)
W_vs_QS_results.sva <- as.data.frame(W_vs_QS_results.sva)
W_vs_QS_results_sig <- subset(W_vs_QS_results.sva, padj < 0.05)
W_vs_QS_results_sig <- as.data.frame(W_vs_QS_results_sig)
length(W_vs_QS_results_sig$padj)#1585
#write.csv(W_vs_QS_results_sig, file = "W_vs_QS_results_sig.csv")

#############################################################################################################################################################
#---------------- Volcano plots
#---------------- Figure 
#############################################################################################################################################################
0 #workers v Reproductive(QueensVsSolRep)

W_vs_QS <- results(dds.3, contrast = list( c("behaviorworker"),c("behaviorqueen", "behaviorsol_reprod") ),listValues = c(1,-1/2),)
summary(W_vs_QS, na.rm=TRUE)

W_vs_QS.df <- as.data.frame(W_vs_QS)

W_vs_QS.df<- W_vs_QS.df %>%         #label the rows
  as_tibble(rownames= "geneID")


#The default cut-off for log2FC is 2; the default cut-off for P value is 10e-6.
lab_italics <- paste0("italic('", rownames(W_vs_QS), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_02909','LBAL_11453', 'LBAL_08120','LBAL_03296','LBAL_09148','LBAL_10596'),
  "')")

EnhancedVolcano(W_vs_QS,
                lab = lab_italics,
                xlim = c(-22, 22), ylim = c(0,35),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = NA,#NULL for all labels
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10] ~ italic(P)),
                title = 'Volcano Plot',
                subtitle = 'Workers vs Reproductive_brain',
                pCutoff = 0.001,
                FCcutoff = 0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('grey26', 'grey26', 'grey26', 'purple'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()




1 #Queen vs workers

Q_W_vol1 <- results(dds.2, contrast = c("behavior","queen","worker"))


summary(Q_W_vol1)

Q_W.df <- as.data.frame(Q_W_vol1)

Q_W.df<- Q_W.df %>%         #label the rows
  as_tibble(rownames= "geneID")


#The default cut-off for log2FC is 2; the default cut-off for P value is 10e-6.
lab_italics <- paste0("italic('", rownames(Q_W_vol1), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_02909','LBAL_11453', 'LBAL_08120','LBAL_03296','LBAL_09148','LBAL_10596'),
  "')")
EnhancedVolcano(Q_W_vol1,
                lab = lab_italics,
                xlim = c(-10, 10), ylim = c(0, 70),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = NA,#NULL for all labels
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10] ~ italic(P)),
                title = 'Volcano Plot',
                subtitle = 'queens vs workers',
                pCutoff = 10e-10,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'green4', 'blue3', 'goldenrod1'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()


2 #Queen vs solitary reproductives

Q_SR_vol2 <- results(dds.2, contrast = c("behavior","queen","foundress"))
summary(Q_SR_vol2)

QSR.df <- as.data.frame(Q_SR_vol2)

QSR.df<- QSR.df %>%         #label the rows
  as_tibble(rownames= "geneID")

lab_italics <- paste0("italic('", rownames(Q_SR_vol2), "')")
selectLab_italics = paste0(
  "italic('",
  c('LLBAL_11437'),"')")

EnhancedVolcano(Q_SR_vol2,
                lab = lab_italics,
                xlim = c(-16, 16), ylim = c(0, 70),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Volcano Plot',
                subtitle = 'queens vs solitary reproductives',
                pCutoff = 10e-10,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'green4', 'blue3', 'darkorchid1'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()



3 #Workers vs solitary reproductives

W_SR_vol3 <- results(dds.2, contrast = c("behavior","worker","sol_reprod"))

summary(W_SR_vol3)

WSR.df <- as.data.frame(W_SR_vol3)

WSR.df<- WSR.df %>%         #label the rows
  as_tibble(rownames= "geneID") 

lab_italics <- paste0("italic('", rownames(W_SR_vol3), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_03018','LBAL_07230','LBAL_08133','LBAL_03921','LBAL_11514','LBAL_11323'),
  "')")

EnhancedVolcano(W_SR_vol3,
                lab = lab_italics,
                xlim = c(-16, 16), ylim = c(0, 70),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Volcano Plot',
                subtitle = 'workers vs Solitary reproductives',
                pCutoff = 10e-10,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'green4', 'blue3', 'darkorange2'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') #+ coord_flip()


4 #Queens vs Foundress

Q_F_vol4<- results(dds.2, contrast = c("behavior","queen","foundress"))

summary(Q_F_vol4)

QF.df <- as.data.frame(Q_F_vol4)

QF.df<- QF.df %>%         #label the rows
  as_tibble(rownames= "geneID") 

lab_italics <- paste0("italic('", rownames(Q_F_vol4), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_04609',''),
  "')")

EnhancedVolcano(Q_F_vol4,
                lab = lab_italics,
                xlim = c(-15, 15), ylim = c(0, 70),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Volcano Plot',
                subtitle = 'queens vs foundress',
                pCutoff = 10e-10,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'green4', 'blue3', 'firebrick'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()

5 #worker vs foundress

F_W_vol5 <- results(dds.2, contrast = c('behavior','foundress','worker'))

summary(F_W_vol5)
FW.df <- as.data.frame(F_W_vol5)

FW.df<- FW.df %>%         #label the rows
  as_tibble(rownames= "geneID")


#The default cut-off for log2FC is 2; the default cut-off for P value is 10e-6.
lab_italics <- paste0("italic('", rownames(F_W_vol5), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_11453'),
  "')")

EnhancedVolcano(F_W_vol5,
                lab = lab_italics,
                xlim = c(-30, 30), ylim = c(0, 70),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Volcano Plot',
                subtitle = 'foundress vs worker',
                pCutoff = 10e-10,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'green4', 'blue3', 'dodgerblue2'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()



6 #foundress vs sol_reproductive

#NO DEGS

F_S_vol6 <- results(dds.2, contrast = c('behavior','foundress', 'sol_reprod'))

summary(F_S_vol6)

FS.df <- as.data.frame(F_S_vol6)

FS.df<- FS.df %>%         #label the rows
  as_tibble(rownames= "geneID")


#The default cut-off for log2FC is 2; the default cut-off for P value is 10e-6.
lab_italics <- paste0("italic('", rownames(F_S_vol6), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_09926'),
  "')")

EnhancedVolcano(F_S_vol6,
                lab = lab_italics,
                xlim = c(-13, 13), ylim = c(0, 70),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Volcano Plot',
                subtitle = 'foundress vs solitary reproductive',
                pCutoff = 10e-10,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'green4', 'blue3', 'turquoise1'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()




#############################################################################################################################################################
#-------------- UPSETR plot 
#-------------- Figure 
#############################################################################################################################################################





Q_W.Genes <-  rownames(Q_W_results.sva_sig)[
  which(p.adjust(Q_W_results.sva_sig$padj, "fdr") <= 0.05)
]

Q_S.Genes <-  rownames(Q_S_results.sva_sig)[
  which(p.adjust(Q_S_results.sva_sig$padj, "fdr") <= 0.05)
]

W_S.Genes <-  rownames(W_S_results.sva_sig)[
  which(p.adjust(W_S_results.sva_sig$padj, "fdr") <= 0.05)
]

Q_F.Genes <-  rownames(Q_F_results.sva_sig)[
  which(p.adjust(Q_F_results.sva_sig$padj, "fdr") <= 0.05)
]

F_W.Genes <-  rownames(F_W_results.sva_sig)[
  which(p.adjust(F_W_results.sva_sig$padj, "fdr") <= 0.05)
]

F_S.Genes <-  rownames(F_S_results.sva_sig)[
  which(p.adjust(F_S_results.sva_sig$padj, "fdr") <= 0.05)
]


length(Q_W.Genes)
length(Q_S.Genes) # 0 DEGs
length(W_S.Genes)
length(Q_F.Genes)
length(F_W.Genes)
length(F_S.Genes) # 1 DEGs


#convert these from a list of character vectors to a presence/absence matrix

UpSetR::upset(fromList(list(QvW = Q_W.Genes,  WvS = W_S.Genes, QvF = Q_F.Genes, FvW = F_W.Genes, FvS = F_S.Genes)),
              nsets = 6,
              point.size = 4, 
              sets.bar.color =c("dodgerblue2","goldenrod1", "darkorange2", "firebrick", "darkorchid1"),
              line.size = 1,
              text.scale=1.5,
              mainbar.y.label = "DEGs",
              order.by = c("freq", "degree"))




#EXTRA


get_upregulated <- function(Q_W_results.sva_sig){
  
  key <- intersect(rownames(Q_W_results.sva_sig)[which(Q_W_results.sva_sig$log2FoldChange>= 1)], rownames(Q_W_results.sva_sig)[which(Q_W_results.sva_sig$padj<=0.05)])
  
  results <- as.data.frame((Q_W_results.sva_sig)[which(rownames(Q_W_results.sva_sig) %in% key),])
  return(results)
}

get_downregulated <- function(Q_W_results.sva_sig){
  
  key <- intersect(rownames(Q_W_results.sva_sig)[which(Q_W_results.sva_sig$log2FoldChange <= -1)], rownames(Q_W_results.sva_sig)[which(Q_W_results.sva_sig$padj<=0.05)])
  
  results <- as.data.frame((Q_W_results.sva_sig)[which(rownames(Q_W_results.sva_sig) %in% key),])
  return(results)
}

up_QW <- get_upregulated (Q_W_results.sva_sig) 
down_QW <- get_downregulated(Q_W_results.sva_sig)
