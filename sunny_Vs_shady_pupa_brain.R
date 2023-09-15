#PART 1(WITHOUT SVA's)
## Load libraries----
library(DESeq2)
library(ggplot2)
library(ggside)
library("pheatmap")
library(RColorBrewer)
library(tidyverse) 
library(sva)
library(ggrepel)
library(viridis)
library(EnhancedVolcano)


#load count files----
setwd("~/Desktop/Paper Witing, Genomic & RNA seqs/Miseq and novaseq data_2021_NEB/Japan data_2021 collections/novaseq/Round 2_counts/subsets")
Counts<-read.delim("sunny&shady_brain new.txt", header = T, row.names = 1, stringsAsFactors = F, check.names=FALSE)
colSums(Counts)
str(Counts)

#Load unfiltered trait data
sample.info.all<- read.csv("sample.info.metadata.csv", header = TRUE )

sample.info <- sample.info.all[which(sample.info.all$sex == 'female' & 
                                       sample.info.all$caste == 'pupa' & 
                                       sample.info.all$tissue == 'brain' &
                                       sample.info.all$total_gene.counts > 0),] 

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
nestID<-factor(sample.info.filt$nestID)
temp_probe<-factor(sample.info.filt$temp_probe)
Location<-factor(sample.info.filt$Location)
Group<-factor(sample.info.filt$Group)

all((sample.info.filt$sample.name) == colnames(filt.counts))


##Buid DESeq2 object----

dds <- DESeqDataSetFromMatrix(countData = filt.counts,
                              colData = sample.info.filt, 
                              design= ~ Location ) #with unfiltered counts

keep <- rowSums(counts(dds) >= 2) >= 10 # 10 samples (smallest sample size for behavior) have counts >=3. ##Go  through each row and determine every value > or= 3.
dim(dds)
dds<-dds[keep,] 
dim(dds) # 11248 genes

#check samples match
colnames(dds) == sample.info.filt$sample.name
length(colnames(dds))

## Normalizing using the method for an object of class 'count data' 
dds <-  estimateSizeFactors(dds)
counts(dds, normalized=TRUE) 

dds <- DESeq(dds)
dim(dds) #11248

# Estimate data reproducibility
dds.rep = estimateDispersions(dds)
plotDispEsts(dds.rep)

#Extract transformed values

#Heatmaps of count matrix----


vsd <- vst(dds, blind=FALSE) 
head(assay(vsd), 3)

# exporting the vst values
vst.df <- as.data.frame(assay(vsd))
vst.df <- t(dds_vst.df.sva)
#write.csv(vst.df, file = "counts_core_for WGCNA_pupa.sunny v shady brn_sva.csv")


hmcol <- viridis(150)
anot_col<- brewer.pal(12,"Paired") 



df<- as.data.frame(colData(dds)["Location"])


annotation_colors <- list(
  Location = c(sunny=anot_col[1], shady=anot_col[6])
)


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:250]

allGenes<-order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)


#!!!!Running this heatmap takes a lot of time
pheatmap(assay(vsd)[select,], #or allGenes.sva
         cluster_rows=TRUE, 
         show_rownames=FALSE, 
         show_colnames = F, 
         cluster_cols=TRUE, 
         annotation=df,
         annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         clustering_distance_cols = "manhattan",
         cutree_rows = 1, #change these
         cutree_cols = 2,
         color=hmcol,
         width = 7,
         height=6)

#Heatmap of sample distances. Overview over similarities and dissimilarities between samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Location, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues"))) (22)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData1<-plotPCA(vsd, intgroup=c( "Location"),returnData = T)##including returnData=T returns the data frame useful for plotting in ggplot


percentVar1 <- round(100 * attr(pcaData1, "percentVar"))# percent variance explained

###PCA type 1----all genes
pca1<-ggplot(pcaData1,aes(PC1,PC2,shape=Location, color= Location))+
  geom_point(size=3.5,color="black")+ 
  # geom_text_repel(aes(label = sample.name)) +
  geom_point(size=3,aes(color=Location))+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  coord_fixed()+
  # stat_ellipse() +
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_hline(yintercept = 0, linetype="dotted")+
  # scale_color_manual(values =anot_col[c(1,11)])+
  theme_bw() +
  scale_x_continuous(limits=c(-40, 60)) + #---->has 4 missing values
  scale_y_continuous(limits=c(-40, 40)) 

gg<-pca1 + scale_color_manual(values=c( "#345c67", "#9e2b2b")) + theme_bw()

gg
#ggsave("Sunny v Shady broods.pdf",pca1, width=3, height=3, units="in", scale=3)

#################
#Regression plot

lm(pcaData1$PC1 ~ vsd$Location, data = pcaData1) |> summary()
plot(pcaData1$PC1 ~ vsd$Location)
lm(pcaData1$PC1 ~ vsd$Location)

# making PCA plot_type 2----

# exporting the vst values for the WGCNA analysis
dds_vst.df <- as.data.frame(assay(vsd))
dds_vst_transposed <- t(dds_vst.df)
#write.csv(dds_vst.df, file = "counts_core_for WGCNA_pupa.brains_brood 1.csv")


######more PCAs
z <- plotPCA(vsd, intgroup = "Location" )
z +theme_bw()   + geom_text_repel(aes(label = name)) 

z + geom_label(aes(label=row.names(sample.info.filt))) + theme_bw()
z$data

z.df<-as.tibble(z$data)


#PCA type 3_ggside
ggplot(z.df, aes(PC1, PC2, color=Location)) + 
  geom_point(size = 2) +
  geom_xsideboxplot(aes(y = Location), orientation = "y") +
  scale_xsidey_discrete() + #In order to use xsideboxplot with a main panel that uses
  #continuous data, force y axis in xside panel to be discrete
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL) +
  theme(ggside.panel.scale = .5) +
  theme_bw() 
# Detailed PCA of small multiple plots

new.z.df<- z$data [,1:2] %>% # PC's 1$2
  as_tibble() %>%
  add_column(Location) %>%
  add_column(sample.name) %>%
  add_column(sociality)


z.pivot<- pivot_longer (new.z.df,
                        cols = PC1:PC2,
                        names_to = "PC",
                        values_to = "loadings")
pca3<-ggplot(z.pivot) +
  aes(x=Location, y=loadings, fill=sociality) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC) +
  labs(title = "PCA of 'small multiples' plot",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
beep(sound = 6)

###Pairwise comparisons Deseq----

ddsMF <- dds
resultsNames(ddsMF)


## comparing sunny and shady
SunnyVshady_results <- results(ddsMF, contrast = c("Location","sunny","shady"),alpha=0.001, lfc=0, pAdjustMethod="BH")
summary(SunnyVshady_results, na.rm=TRUE)
SunnyVshady_results <- as.data.frame(SunnyVshady_results)
SunnyVshady_results_sig <- subset(SunnyVshady_results, padj < 0.001)
SunnyVshady_results_sig <- as.data.frame(SunnyVshady_results_sig)
length(SunnyVshady_results_sig$padj)#1121


## saving the DEG files for each comparison

setwd()

#write.csv(SunnyVshady_results_sig, file = "SunnyVshady_pupa_brains_brood 1.csv")


#PLOT VOLCANO PLOTS----

#sunny vs shady

SunnyVshady_results <- results(ddsMF, contrast = c("Location","sunny","shady"))

summary(SunnyVshady_results)

SunnyVshady.df <- as.data.frame(SunnyVshady_results)

SunnyVshady.df<- SunnyVshady.df %>%         #label the rows
  as_tibble(rownames= "geneID")


lab_italics <- paste0("italic('", rownames(SunnyVshady_results), "')")
selectLab_italics = paste0(
  "italic('",
  c('LBAL_00243','LBAL_03711','LBAL_12291','LBAL_04226'),
  "')")

EnhancedVolcano(SunnyVshady_results,
                lab = lab_italics,
                xlim = c(-9, 9), ylim = c(0, 35),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Volcano Plot',
                subtitle = 'sunny v shady',
                pCutoff = 0.001,
                FCcutoff = 0,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('grey26', 'grey26', 'grey26', 'mediumvioletred'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')# + coord_flip()



