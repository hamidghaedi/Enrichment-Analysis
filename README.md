# Enrichment-Analysis
R codes I am using for getting from RNA-seq raw count to Pathways


## Introduction
Assume we performed an RNA-seq (or microarray gene expression) experiment and now want to know what pathway/biological process shows enrichment for our [differentially expressed] genes. There are terminologies someone may need to be familiar with before diving into the topic thoroughly: 

**Gene set**

A gene set is an unordered collection of genes that are functional related. 

**Pathway**

A pathway can be interpreted as a gene set by ignoring functional relationships among genes.

**Gene Ontology (GO)**

GO describe gene function. A gene role/function could be attributed to the three main classes: 

*•	Molecular Function (MF)* : which define molecular activities of gene products.

*•	Cellular Component (CC)* : which describe where gene products are active/localized.

*•	Biological Process (BP)* : which describe pathways and processes that the gene product is active in.

**Kyoto Encyclopedia of Genes and Genomes (KEGG)**

KEGG is a collection of manually curated pathway maps representing molecular interaction and reaction networks. 



There are different methods widely used for functional enrichment analysis: 

**1-	Over Representation Analysis (ORA):**

This is the simplest version of enrichment analysis and at the same time the most widely used approach. The concept in this approach is based on a Fisher exact test p value in a contingency table. For example, supposed we come up with 160 differentially expressed genes from a microarray expression experiment which was able to investigate 17,000 gene expression levels. We found that 30 genes from our findings are somehow member of a pathway which has 50 gene members (call it pathway X). To find enrichment we can perform a Fisher exact test on the following contingency table:

|                 | not.intrested.genes   |intrested.genes |
| :-------------- |:---------------------:|---------------:|
| in pathwayX     | 20                    | 30             |
| not.in pathwayX | 16820                 | 130            |

There are relatively large number of web-tools R package for ORA. Personally I am a fan of [DAVID](https://david.ncifcrf.gov/home.jsp) webtools however its lat update was in 2016 (DAVID 6.8 Oct. 2016).  

**2-	Gene Set Enrichment Analysis (GSEA):**

It was developed by Broad Institute. This is the preferred method when genes are coming from an expression experiment like microarray and RNA-seq. However, the original methodology was designed to work on microarray but later modification made it suitable for RNA-seq also. In this approach, you need to rank your genes based on a statistic (like what DESeq2 provide), and then perform enrichment analysis against different pathways (= gene set). You have to download the gene set files into you local system. The point is that here the algorithm will use all genes in the ranked list for enrichment analysis. [in contrast to ORA where only genes passed a specific threshold (like DE ones) would be used for enrichment analysis]. You can find more details about the methodology on the original [PNAS paper](https://www.pnas.org/content/102/43/15545.abstract), here is a summary of why one should use this approach instead of ORA:

1- After correcting for multiple hypotheses testing, no individual gene may meet the threshold for statistical significance.

2- On the other hand, one may be left with a long list of statistically significant genes without any unifying biological theme.

3- Cellular processes often affect sets of genes acting in concert, using ORA may lead to miss important effects on pathways.

GSEA software maybe finds on its [homepage](https://www.gsea-msigdb.org/gsea/index.jsp). However, there are some Bioconductor packages which use a similar approach to do GSEA, I like to use this one : [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html). Also there are some R package which can do ROA and GSEA for you like [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html). 

In the analysis after getting a ranked gene from diffrential expression analysis, we need to have gene lists for GSEA. The Molecular Signatures Database (MSigDB) is a collection of annotated gene sets for use with GSEA software and possibly those works like GSEA. The MSigDB gene sets are divided into 9 major collections:

H: hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.

C1: positional gene sets  for each human chromosome and cytogenetic band.

C2: curated gene sets  from online pathway databases, publications in PubMed, and knowledge of domain experts.

C3: regulatory target gene sets  based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.

C4: computational gene sets  defined by mining large collections of cancer-oriented microarray data.

C5: ontology gene sets  consist of genes annotated by the same ontology term.

C6: oncogenic signature gene sets  defined directly from microarray gene expression data from cancer gene perturbations.

C7: immunologic signature gene sets  defined directly from microarray gene expression data from immunologic studies.

C8: cell type signature gene sets  curated from cluster markers identified in single-cell sequencing studies of human tissue.

To download these gene sets in a folder go the MSigDB [website](https://www.gsea-msigdb.org/gsea/login.jsp#msigdb), registe with your email and download the data. 

```R
#_______________Loading packages______________________________#
library(DESeq2)
library(org.Hs.eg.db)
library(tibble)
library(dplyr)
library(tidyr)
library(fgsea)
library(ggplot2)
library(reshape2)
 
#________________diffrential expression analysis______________#
# reading expression data matrix and getting rid of what we dont like

# my data is raw count data comming from RNA-seq on 476 bladder cancer sample
rna <- read.table("Uromol1_CountData.v1.csv", header = T, sep = ",")

dim(rna)
#  [1] 60488   477
head(rna, 5)

#  genes               U0001 U0002 U0006 U0007 U0010 U0011 U0012 U0013 U0015 U0018 U0023 U0024 U0026 U0027 U0028
#1 ENSG00000000003.13  1458   228  1800  3945   293  1737  6362  8869  1728  4822  5913   103  2138  4196  4949
#2 ENSG00000000005.5     0     0     9     0     0     0     0     5     2     2     2     0     2     4     2

# setting row name for expression matrix
rownamefordata <- substr(rna$genes,1,15)
rownames(rna) <- rownamefordata
rm(rownamefordata)
# removing column with gene name
rna <- rna[, -1]

# reading clinical data 
# read uromol clinical
uromol_clin <- read.table("uromol_clinic.csv", sep = ",", header = T)
names(uromol_clin)

#[1] "UniqueID"                           "CLASS"                              "BASE47"                            
#[4] "CIS"                                "X12.gene.signature"                 "Lund"                              
#[7] "Stage"                              "Grade"                              "EORTC.risk.score..NMIBC."          
#[10] "sex"                                "age"                                "TumorSize"                         
#[13] "GrowthPattern"                      "BCGTreatment"                       "CIS.in.disease.course"             
#[16] "Cystectomy"                         "Progression.to.T2."                 "Progression.free.survival..months."

#Setting row name for clinical data
ids <- uromol_clin$UniqueID
rownames(uromol_clin) <- ids
rm(ids)

# I want to compare samples with Ta stage against those in higher stage (non_Ta)
# making column for comparison

table(uromol_clin$Stage)

# CIS   T1 T2-4   Ta 
# 3     112   16  345 

# selecting only Ta samples
uromol_clin$isTa <- as.factor(ifelse(uromol_clin$Stage == "Ta", "Ta", "non_Ta"))
levels(uromol_clin$isTa)
# [1] "non-Ta" "Ta" 

# making sure dataset are compatible regrding to sample order
rna <- rna[, row.names(uromol_clin)]
all(rownames(uromol_clin) %in% colnames(rna))

#_______Making_Expression_Object__________#
#Making DESeqDataSet object which stores all experiment data
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = uromol_clin,
                              design = ~ isTa)

# prefilteration: it is not necessary but recommended to filter out low expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rm(keep)
# data tranfromation
vsd <- vst(dds, blind=FALSE)

#________________DE_analysis_____________#
dds <- DESeq(dds) #This would take some time
res <- results(dds, alpha=0.05)
summary(res)

#_________________GSEA___________________#
# Steps toward doing gene set emrichment analysis (GSEA):

# 1- obtaining a stats for ranking genes in your experiment,
# 2- creating a named vector out of the DESeq2 result
# 3- Obtaining a gene set from mysigbd
# 4- doing analysis


# already we performed DESeq2 analysis and have statitics for workig on it
res$row <- rownames(res)

# important notice: if you have not such stats in your result (say comming from edgeR),
# you may need to create a rank metric for your genes. To do this:
# metric = -log10(pvalue)/sign(log2FC)



# Map Ensembl gene IDs to symbol. First create a mapping table.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
names(ens2symbol)[1] <- "row"

ens2symbol <- as_tibble(ens2symbol)
ens2symbol
# joining
res <- merge(data.frame(res), ens2symbol, by=c("row"))

# remove the NAs, averaging statitics for a multi-hit symbol
res2 <- res %>% 
  select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

# creating  a named vector [ranked genes]
ranks <- res2$stat
names(ranks) <- res2$SYMBOL


# Load the pathway (gene set) into a named list
# downloaded mysigdb were located in my "~" directory:
pathways.hallmark <- gmtPathways("~/mysigdb/h.all.v7.2.symbols.gmt")

# show few lines from the pathways file
head(pathways.hallmark)

#Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

# To see what genes are in each of these pathways:
gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  inner_join(res, by="SYMBOL")

#______________________VISUALIZATION______________________________#

#__________bar plot _______________#
# Plot the normalized enrichment scores. 
#Color the bar indicating whether or not the pathway was significant:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Hallmark pathways Enrichment Score from GSEA")
  ```
  ![alt text](https://github.com/hamidghaedi/Enrichment-Analysis/blob/main/gsea.png)

```R
#__________Enrichment  Plot_______#
# Enrichment plot for E2F target gene set
plotEnrichment(pathway = pathways.hallmark[["HALLMARK_E2F_TARGETS"]], ranks)
```
![alt text](https://github.com/hamidghaedi/Enrichment-Analysis/blob/main/gse2.png)
```R
#
plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes, 
                gseaParam=0.5)
 ```
 ![alt text](https://github.com/hamidghaedi/Enrichment-Analysis/blob/main/gse3.png)
 
 ```R
#________ Heatmap Plot_____________#
# pathways with significant enrichment score
sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
sig.gen <- unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))

### create a new data-frame that has '1' for when a gene is part of a term, and '0' when not
h.dat <- dcast(gene.in.pathway[, c(1,2)], SYMBOL~pathway)
rownames(h.dat) <- h.dat$SYMBOL
h.dat <- h.dat[, -1]

h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
h.dat <- h.dat[, colnames(h.dat) %in% sig.path]

# keep those genes with 3  or more occurnes
table(data.frame(rowSums(h.dat)))

# 1       2    3    4    5    6 
# 1604  282   65   11    1    1 
h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]


```

## Refrences
1- [clusterProfiler: universal enrichment tool for functional and comparative study](http://yulab-smu.top/clusterProfiler-book/)

2- [Fast Gene Set Enrichment Analysis](https://bioconductor.org/packages/release/bioc/html/fgsea.html)

3- [Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles](https://www.pnas.org/content/102/43/15545.abstract)

4- [DAVID Bioinformatics Resources 6.8](https://david.ncifcrf.gov/home.jsp)

5- [DESeq results to pathways in 60 Seconds with the fgsea package](https://stephenturner.github.io/deseq-to-fgsea/)

6- [Rank-rank hypergeometric overlap: identification of statistically significant overlap between gene-expression signatures](https://pubmed.ncbi.nlm.nih.gov/20660011/)


