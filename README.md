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



## Refrences
1- [clusterProfiler: universal enrichment tool for functional and comparative study](http://yulab-smu.top/clusterProfiler-book/)

2- [Fast Gene Set Enrichment Analysis](https://bioconductor.org/packages/release/bioc/html/fgsea.html)

3- [Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles](https://www.pnas.org/content/102/43/15545.abstract)

4- [DAVID Bioinformatics Resources 6.8](https://david.ncifcrf.gov/home.jsp)

5- [DESeq results to pathways in 60 Seconds with the fgsea package](https://stephenturner.github.io/deseq-to-fgsea/)

6- [Rank-rank hypergeometric overlap: identification of statistically significant overlap between gene-expression signatures](https://pubmed.ncbi.nlm.nih.gov/20660011/)


