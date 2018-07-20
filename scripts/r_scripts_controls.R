# load DESeq2 library and functions:
library("DESeq2")

##########################################################################
###### all >12hr Oocerea controls (genes) ###############################
##########################################################################

# set all output to be written to file:
zz <- file("deseq2_analysis_meth_bs_controls_no_outliers.R.out", open = "wt")
sink(zz)
sink(zz, type = "output")

# set location of HTSeq count files:
directory <- "/Volumes/antqueen/genomics/experiments/analyses/PRO20150424_star_all_cerapachys/htseq"

# extract sample files:
sampleFiles <- grep("[PM].*(_SP|_FL).*htseq.gene.out", list.files(directory), value=TRUE)
sampleFiles

# assign phase and batch to each sample, based on filename:
sampleCondition <- sub(".*(_SP|_FL).*", "\\1", sampleFiles)
sampleCondition

sampleBatch <- sub("[PM]_(B[0-9]).*", "\\1", sampleFiles)
sampleBatch


# construct table, including design (incorporates batch effect in model)
sampleTable <- data.frame(sampleName = sampleFiles,
                           fileName = sampleFiles,
                           condition = sampleCondition,
                           batch = sampleBatch)

### PERFORM TESTS ###
# standard
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        design= ~ batch + condition)
dds
dds_std <- DESeq(dds, minReplicatesForReplace=5)

# look at results of tests:
resQF <- results(dds_std, alpha=0.05, contrast=c("condition", "_SP", "_FL"))


# produce and save the variance adjusted count values for each gene and save to file:
vsd <- varianceStabilizingTransformation(dds_std)
vstMat <- assay(vsd)
write.table(as.data.frame(vstMat), file="Cerapachys.genes.all_samples.vst.tbl", quote=FALSE)

# order results on adjusted pvalue:
resQFOrdered <- resQF[order(resQF$padj),]
head(resQFOrdered)
summary(resQF)



# count number of genes with adjusted pvalue less than 0.05
SigQF <- subset(resQFOrdered, padj < 0.05)
length(SigQF$padj)

# graph all variance adjusted genes and save to file:
jpeg("gene_variance_controls_meth_bs.jpeg")
plotMA(resQF, main="statary and foraging Cerapachys Samples (Genes)", ylim=c(-2,2))
dev.off()

# save pvalue and log fold change results for all genes to file:
write.table(as.data.frame(resQFOrdered), file="Cerapachys.genes.all_samples.pvalues.out", quote=FALSE)


# stop reporting to file. If reporting is restarted, don't forget to add append = "TRUE" !!!
sink(zz, type = "output")
sink()


##########################################################################
###### all Methylation Oocerea controls (genes) #######################
##########################################################################

# set all output to be written to file:
zz <- file("deseq2_analysis_meth_controls.R.out", open = "wt")
sink(zz)
sink(zz, type = "output")

# set location of HTSeq count files:
directory <- "/Volumes/antqueen/genomics/experiments/analyses/PRO20150424_star_all_cerapachys/htseq"

# extract sample files:
sampleFiles <- grep("[M].*(_SP|_FL).*htseq.gene.out", list.files(directory), value=TRUE)
sampleFiles

# assign phase and batch to each sample, based on filename:
sampleCondition <- sub(".*(_SP|_FL).*", "\\1", sampleFiles)
sampleCondition

sampleBatch <- sub("[M]_(B[0-9]).*", "\\1", sampleFiles)
sampleBatch


# construct table, including design (incorporates batch effect in model)
sampleTable <- data.frame(sampleName = sampleFiles,
                           fileName = sampleFiles,
                           condition = sampleCondition,
                           batch = sampleBatch)

### PERFORM TESTS ###
# standard
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        design= ~ batch + condition)
dds
dds_std <- DESeq(dds, minReplicatesForReplace=5)

# look at results of tests:
resQF <- results(dds_std, alpha=0.05, contrast=c("condition", "_SP", "_FL"))


# produce and save the variance adjusted count values for each gene and save to file:
vsd <- varianceStabilizingTransformation(dds_std)
vstMat <- assay(vsd)
write.table(as.data.frame(vstMat), file="Cerapachys.genes.all_meth_samples.vst.tbl", quote=FALSE)

# order results on adjusted pvalue:
resQFOrdered <- resQF[order(resQF$padj),]
head(resQFOrdered)
summary(resQF)



# count number of genes with adjusted pvalue less than 0.05
SigQF <- subset(resQFOrdered, padj < 0.05)
length(SigQF$padj)

# graph all variance adjusted genes and save to file:
jpeg("gene_variance_controls_meth.jpeg")
plotMA(resQF, main="statary and foraging Cerapachys Samples (Genes)", ylim=c(-2,2))
dev.off()

# save pvalue and log fold change results for all genes to file:
write.table(as.data.frame(resQFOrdered), file="Cerapachys.genes.all_meth_samples.pvalues.out", quote=FALSE)


# stop reporting to file. If reporting is restarted, don't forget to add append = "TRUE" !!!
sink(zz, type = "output")
sink()



##########################################################################
###### all broodswap >12 hr Oocerea controls (genes) ####################
##########################################################################

# set all output to be written to file:
zz <- file("deseq2_analysis_12hr_bs_controls.R.out", open = "wt")
sink(zz)
sink(zz, type = "output")

# set location of HTSeq count files:
directory <- "/Volumes/antqueen/genomics/experiments/analyses/PRO20150424_star_all_cerapachys/htseq"

# extract sample files:
sampleFiles <- grep("[P].*(_SP|_FL).*htseq.gene.out", list.files(directory), value=TRUE)
sampleFiles

# assign phase and batch to each sample, based on filename:
sampleCondition <- sub(".*(_SP|_FL).*", "\\1", sampleFiles)
sampleCondition

sampleBatch <- sub("[P]_(B[0-9]).*", "\\1", sampleFiles)
sampleBatch


# construct table, including design (incorporates batch effect in model)
sampleTable <- data.frame(sampleName = sampleFiles,
                           fileName = sampleFiles,
                           condition = sampleCondition,
                           batch = sampleBatch)

### PERFORM TESTS ###
# standard
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        design= ~ batch + condition)
dds
dds_std <- DESeq(dds, minReplicatesForReplace=5)

# look at results of tests:
resQF <- results(dds_std, alpha=0.05, contrast=c("condition", "_SP", "_FL"))


# produce and save the variance adjusted count values for each gene and save to file:
vsd <- varianceStabilizingTransformation(dds_std)
vstMat <- assay(vsd)
write.table(as.data.frame(vstMat), file="Cerapachys.genes.all_12hr_bs_samples.vst.tbl", quote=FALSE)

# order results on adjusted pvalue:
resQFOrdered <- resQF[order(resQF$padj),]
head(resQFOrdered)
summary(resQF)



# count number of genes with adjusted pvalue less than 0.05
SigQF <- subset(resQFOrdered, padj < 0.05)
length(SigQF$padj)

# graph all variance adjusted genes and save to file:
jpeg("gene_variance_controls_12hr_bs.jpeg")
plotMA(resQF, main="statary and foraging Cerapachys Samples (Genes)", ylim=c(-2,2))
dev.off()

# save pvalue and log fold change results for all genes to file:
write.table(as.data.frame(resQFOrdered), file="Cerapachys.genes.all_12hr_bs_samples.pvalues.out", quote=FALSE)


# stop reporting to file. If reporting is restarted, don't forget to add append = "TRUE" !!!
sink(zz, type = "output")
sink()

