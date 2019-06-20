Gene Set Enrichment Analysis with ClusterProfiler
================
Mohammed Khalfan
5/19/2019

This R Notebook describes the implementation of gene set enrichment analysis (GSEA) using the clusterProfiler package. For more information please see the full documentation here: <https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html>

Install and load required packages
==================================

``` r
#BiocManager::install("clusterProfiler", version = "3.8")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
```

Annotations
===========

I'm using *D melanogaster* data, so I install and load the annotation "org.Dm.eg.db" below. See all annotations available here: <http://bioconductor.org/packages/release/BiocViews.html#___OrgDb> (there are 19 presently available).

``` r
# SET THE DESIRED ORGANISM HERE
organism = "org.Dm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
```

Prepare Input
=============

``` r
# reading in data from deseq2
df = read.csv("drosphila_example_de.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
```

Gene Set Enrichment
-------------------

Params:

**keyType** This is the source of the annotation (gene ids). The options vary for each annotation. In the example of *org.Dm.eg.db*, the options are:

"ACCNUM" "ALIAS" "ENSEMBL" "ENSEMBLPROT" "ENSEMBLTRANS" "ENTREZID"
"ENZYME" "EVIDENCE" "EVIDENCEALL" "FLYBASE" "FLYBASECG" "FLYBASEPROT"
"GENENAME" "GO" "GOALL" "MAP" "ONTOLOGY" "ONTOLOGYALL"
"PATH" "PMID" "REFSEQ" "SYMBOL" "UNIGENE" "UNIPROT"

Check which options are available with the `keytypes` command, for example `keytypes(org.Dm.eg.db)`.

**ont** one of "BP", "MF", "CC" or "ALL"
**nPerm** permutation numbers, the higher the number of permutations you set, the more accurate your results is, but it will also cost longer time for running permutation.
**minGSSize** minimal size of each geneSet for analyzing.
**maxGSSize** maximal size of genes annotated for testing.
**pvalueCutoff** pvalue Cutoff.
**pAdjustMethod** one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

``` r
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
```

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in fgsea(pathways = geneSets, stats = geneList, nperm = nPerm, minSize = minGSSize, : There are ties in the preranked stats (10.64% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...

    ## done...

Output
======

Table of results
----------------

``` r
head(gse)
```

    ##            ONTOLOGY         ID
    ## GO:0031226       CC GO:0031226
    ## GO:0005887       CC GO:0005887
    ## GO:0060089       MF GO:0060089
    ## GO:0004888       MF GO:0004888
    ## GO:0007186       BP GO:0007186
    ## GO:0004930       MF GO:0004930
    ##                                             Description setSize
    ## GO:0031226       intrinsic component of plasma membrane     448
    ## GO:0005887        integral component of plasma membrane     435
    ## GO:0060089                molecular transducer activity     411
    ## GO:0004888    transmembrane signaling receptor activity     321
    ## GO:0007186 G protein-coupled receptor signaling pathway     211
    ## GO:0004930          G protein-coupled receptor activity     113
    ##            enrichmentScore       NES       pvalue     p.adjust    qvalues
    ## GO:0031226      -0.4155733 -1.684353 0.0001297185 0.0001297185 0.09959456
    ## GO:0005887      -0.4253170 -1.719741 0.0001302253 0.0001302253 0.09959456
    ## GO:0060089      -0.3840571 -1.544868 0.0001317349 0.0001317349 0.09959456
    ## GO:0004888      -0.4112261 -1.621466 0.0001368363 0.0001368363 0.09959456
    ## GO:0007186      -0.5276227 -1.990696 0.0001460707 0.0001460707 0.09959456
    ## GO:0004930      -0.5803457 -2.016342 0.0001558118 0.0001558118 0.09959456
    ##            rank                   leading_edge
    ## GO:0031226 1351  tags=22%, list=9%, signal=20%
    ## GO:0005887 1351  tags=22%, list=9%, signal=21%
    ## GO:0060089 1458 tags=18%, list=10%, signal=17%
    ## GO:0004888 1458 tags=21%, list=10%, signal=19%
    ## GO:0007186 1344  tags=28%, list=9%, signal=26%
    ## GO:0004930 1016  tags=32%, list=7%, signal=30%
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            core_enrichment
    ## GO:0031226 FBgn0040507/FBgn0036278/FBgn0000037/FBgn0027843/FBgn0032006/FBgn0263916/FBgn0037546/FBgn0051146/FBgn0003861/FBgn0261794/FBgn0015380/FBgn0037167/FBgn0052600/FBgn0259245/FBgn0010395/FBgn0037386/FBgn0036125/FBgn0263116/FBgn0260971/FBgn0035331/FBgn0040726/FBgn0005614/FBgn0028875/FBgn0050361/FBgn0038880/FBgn0036043/FBgn0038309/FBgn0040506/FBgn0024963/FBgn0040238/FBgn0260446/FBgn0022710/FBgn0263289/FBgn0265575/FBgn0028704/FBgn0039294/FBgn0031275/FBgn0265413/FBgn0034136/FBgn0004573/FBgn0033137/FBgn0004622/FBgn0004620/FBgn0085409/FBgn0053531/FBgn0053310/FBgn0261046/FBgn0000535/FBgn0033932/FBgn0265416/FBgn0015399/FBgn0038140/FBgn0050106/FBgn0032593/FBgn0024944/FBgn0032456/FBgn0085384/FBgn0035170/FBgn0266137/FBgn0263131/FBgn0051646/FBgn0038653/FBgn0085431/FBgn0085386/FBgn0011582/FBgn0038542/FBgn0034045/FBgn0051708/FBgn0085395/FBgn0032151/FBgn0261574/FBgn0000547/FBgn0039396/FBgn0038980/FBgn0033135/FBgn0016032/FBgn0004456/FBgn0038498/FBgn0264002/FBgn0037419/FBgn0051814/FBgn0037408/FBgn0036260/FBgn0033744/FBgn0264000/FBgn0035385/FBgn0266429/FBgn0030723/FBgn0266758/FBgn0260753/FBgn0019985/FBgn0264386/FBgn0004619/FBgn0011829/FBgn0039419/FBgn0052843
    ## GO:0005887             FBgn0040507/FBgn0036278/FBgn0000037/FBgn0032006/FBgn0263916/FBgn0037546/FBgn0051146/FBgn0003861/FBgn0261794/FBgn0015380/FBgn0037167/FBgn0052600/FBgn0259245/FBgn0010395/FBgn0037386/FBgn0036125/FBgn0263116/FBgn0260971/FBgn0035331/FBgn0040726/FBgn0005614/FBgn0028875/FBgn0050361/FBgn0038880/FBgn0036043/FBgn0038309/FBgn0040506/FBgn0024963/FBgn0040238/FBgn0260446/FBgn0022710/FBgn0263289/FBgn0265575/FBgn0028704/FBgn0039294/FBgn0031275/FBgn0265413/FBgn0034136/FBgn0004573/FBgn0033137/FBgn0004622/FBgn0004620/FBgn0085409/FBgn0053531/FBgn0053310/FBgn0261046/FBgn0000535/FBgn0033932/FBgn0265416/FBgn0015399/FBgn0038140/FBgn0050106/FBgn0032593/FBgn0024944/FBgn0032456/FBgn0085384/FBgn0035170/FBgn0266137/FBgn0263131/FBgn0051646/FBgn0038653/FBgn0085431/FBgn0085386/FBgn0011582/FBgn0038542/FBgn0034045/FBgn0051708/FBgn0085395/FBgn0032151/FBgn0261574/FBgn0000547/FBgn0039396/FBgn0038980/FBgn0033135/FBgn0016032/FBgn0004456/FBgn0038498/FBgn0264002/FBgn0037419/FBgn0051814/FBgn0037408/FBgn0036260/FBgn0033744/FBgn0264000/FBgn0035385/FBgn0266429/FBgn0030723/FBgn0266758/FBgn0260753/FBgn0019985/FBgn0264386/FBgn0004619/FBgn0011829/FBgn0039419/FBgn0052843
    ## GO:0060089                                                                                                                                                                                                                                                                                     FBgn0031055/FBgn0261929/FBgn0033043/FBgn0016696/FBgn0036278/FBgn0000037/FBgn0032006/FBgn0037546/FBgn0051146/FBgn0000464/FBgn0038840/FBgn0015380/FBgn0004369/FBgn0263116/FBgn0035331/FBgn0041243/FBgn0052693/FBgn0028875/FBgn0050361/FBgn0038880/FBgn0024963/FBgn0050340/FBgn0260446/FBgn0040321/FBgn0053639/FBgn0031275/FBgn0004842/FBgn0004573/FBgn0004622/FBgn0052547/FBgn0004620/FBgn0085409/FBgn0053531/FBgn0033932/FBgn0035382/FBgn0038140/FBgn0050106/FBgn0024944/FBgn0034012/FBgn0085384/FBgn0266137/FBgn0028956/FBgn0038653/FBgn0085431/FBgn0085386/FBgn0011582/FBgn0036150/FBgn0038542/FBgn0025680/FBgn0000119/FBgn0034013/FBgn0030437/FBgn0032151/FBgn0039396/FBgn0038980/FBgn0038902/FBgn0037324/FBgn0033404/FBgn0264002/FBgn0031770/FBgn0037408/FBgn0036260/FBgn0033744/FBgn0264000/FBgn0035385/FBgn0266429/FBgn0260753/FBgn0019985/FBgn0004619/FBgn0000489/FBgn0011829/FBgn0039419/FBgn0052843
    ## GO:0004888                                                                                                                                                                                                                                                                                                                                                                                     FBgn0031055/FBgn0261929/FBgn0033043/FBgn0036278/FBgn0000037/FBgn0032006/FBgn0037546/FBgn0000464/FBgn0038840/FBgn0015380/FBgn0004369/FBgn0263116/FBgn0035331/FBgn0041243/FBgn0052693/FBgn0028875/FBgn0050361/FBgn0038880/FBgn0024963/FBgn0050340/FBgn0260446/FBgn0053639/FBgn0031275/FBgn0004842/FBgn0004573/FBgn0004622/FBgn0052547/FBgn0004620/FBgn0085409/FBgn0053531/FBgn0033932/FBgn0035382/FBgn0038140/FBgn0050106/FBgn0024944/FBgn0085384/FBgn0266137/FBgn0028956/FBgn0038653/FBgn0011582/FBgn0036150/FBgn0038542/FBgn0025680/FBgn0000119/FBgn0034013/FBgn0030437/FBgn0032151/FBgn0039396/FBgn0038980/FBgn0037324/FBgn0033404/FBgn0264002/FBgn0031770/FBgn0037408/FBgn0036260/FBgn0033744/FBgn0264000/FBgn0035385/FBgn0266429/FBgn0260753/FBgn0019985/FBgn0004619/FBgn0011829/FBgn0039419/FBgn0052843
    ## GO:0007186                                                                                                                                                                                                                                                                                                                                                                                                                                                             FBgn0036278/FBgn0050054/FBgn0000037/FBgn0037546/FBgn0003861/FBgn0011581/FBgn0263116/FBgn0035331/FBgn0050361/FBgn0038880/FBgn0040506/FBgn0001104/FBgn0050340/FBgn0260446/FBgn0022710/FBgn0000253/FBgn0053639/FBgn0031275/FBgn0004842/FBgn0004573/FBgn0004622/FBgn0052547/FBgn0033932/FBgn0027794/FBgn0038140/FBgn0050106/FBgn0024944/FBgn0266137/FBgn0263131/FBgn0028956/FBgn0010223/FBgn0038653/FBgn0011582/FBgn0038542/FBgn0086704/FBgn0004784/FBgn0265959/FBgn0030437/FBgn0267252/FBgn0013767/FBgn0039396/FBgn0038980/FBgn0032048/FBgn0037976/FBgn0264002/FBgn0031770/FBgn0037408/FBgn0035092/FBgn0036260/FBgn0033744/FBgn0035385/FBgn0266429/FBgn0085413/FBgn0260753/FBgn0019985/FBgn0045038/FBgn0261549/FBgn0039419/FBgn0052843
    ## GO:0004930                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             FBgn0035331/FBgn0050361/FBgn0038880/FBgn0050340/FBgn0260446/FBgn0053639/FBgn0031275/FBgn0004842/FBgn0004573/FBgn0004622/FBgn0052547/FBgn0033932/FBgn0038140/FBgn0050106/FBgn0024944/FBgn0266137/FBgn0028956/FBgn0038653/FBgn0011582/FBgn0038542/FBgn0025680/FBgn0030437/FBgn0039396/FBgn0038980/FBgn0264002/FBgn0031770/FBgn0037408/FBgn0036260/FBgn0033744/FBgn0035385/FBgn0266429/FBgn0260753/FBgn0019985/FBgn0039419/FBgn0052843

Dotplot
-------

``` r
require(DOSE)
```

    ## Loading required package: DOSE

    ## Warning: package 'DOSE' was built under R version 3.5.2

    ## DOSE v3.8.2  For help: https://guangchuangyu.github.io/DOSE
    ## 
    ## If you use DOSE in published research, please cite:
    ## Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609

``` r
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-6-1.png)

Encrichment plot map:
---------------------

Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional modules.

``` r
emapplot(gse, showCategory = 10)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-7-1.png)

Category Netplot
----------------

The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network (helpful to see which genes are involved in enriched pathways and genes that may belong to multiple annotation categories).

``` r
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-8-1.png)

Ridgeplot
---------

Helpful to interpret up/down-regulated pathways.

``` r
ridgeplot(gse) + labs(x = "enrichment distribution")
```

    ## Picking joint bandwidth of 0.698

![](gsea-1_files/figure-markdown_github/unnamed-chunk-9-1.png)

GSEA Plot
---------

Traditional method for visualizing GSEA result.

Params:
**Gene Set** Integer. Corresponds to gene set in the gse object. The first gene set is 1, second gene set is 2, etc.

``` r
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-10-1.png)

PubMed trend of enriched terms
------------------------------

Plots the number/proportion of publications trend based on the query result from PubMed Central.

``` r
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-11-1.png)

KEGG Gene Set Enrichment Analysis
=================================

For KEGG pathway enrichment using the `gseKEGG()` function, we need to convert id types. We can use the `bitr` function for this (included in clusterProfiler). It is normal for this call to produce some messages / warnings.

In the `bitr` function, the param `fromType` should be the same as `keyType` from the `gseGO` function above (the annotation source). This param is used again in the next two steps: creating `dedup_ids` and `df2`.

`toType` in the `bitr` function has to be one of the available options from `keyTypes(org.Dm.eg.db)` and must map to one of 'kegg', 'ncbi-geneid', 'ncib-proteinid' or 'uniprot' because `gseKEGG()` only accepts one of these 4 options as it's `keytype` parameter. In the case of org.Dm.eg.db, none of those 4 types are available, but 'ENTREZID' are the same as ncbi-geneid for org.Dm.eg.db so we use this for `toType`.

As our intial input, we use `original_gene_list` which we created above.

Prepare Input
-------------

``` r
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## Warning in bitr(names(original_gene_list), fromType = "ENSEMBL", toType =
    ## "ENTREZID", : 21.84% of input gene IDs are fail to map...

``` r
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
```

Create gseKEGG object
---------------------

**organism** KEGG Organism Code: The full list is here: <https://www.genome.jp/kegg/catalog/org_list.html> (need the 3 letter code). I define this as `kegg_organism` first, because it is used again below when making the pathview plots.
**nPerm** permutation numbers, the higher the number of permutations you set, the more accurate your results is, but it will also cost longer time for running permutation.
**minGSSize** minimal size of each geneSet for analyzing.
**maxGSSize** maximal size of genes annotated for testing.
**pvalueCutoff** pvalue Cutoff.
**pAdjustMethod** one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
**keyType** one of 'kegg', 'ncbi-geneid', 'ncib-proteinid' or 'uniprot'.

``` r
kegg_organism = "dme"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
```

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in fgsea(pathways = geneSets, stats = geneList, nperm = nPerm, minSize = minGSSize, : There are ties in the preranked stats (7.69% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgsea(pathways = geneSets, stats = geneList, nperm = nPerm,
    ## minSize = minGSSize, : There are duplicate gene names, fgsea may produce
    ## unexpected results

    ## leading edge analysis...

    ## done...

``` r
head(kk2, 10)
```

    ##                ID                               Description setSize
    ## dme00053 dme00053         Ascorbate and aldarate metabolism      31
    ## dme04080 dme04080   Neuroactive ligand-receptor interaction      50
    ## dme04310 dme04310                     Wnt signaling pathway      87
    ## dme00511 dme00511                  Other glycan degradation      21
    ## dme04130 dme04130 SNARE interactions in vesicular transport      20
    ## dme00330 dme00330           Arginine and proline metabolism      49
    ## dme00380 dme00380                     Tryptophan metabolism      20
    ## dme00071 dme00071                    Fatty acid degradation      32
    ## dme00830 dme00830                        Retinol metabolism      31
    ## dme00040 dme00040  Pentose and glucuronate interconversions      45
    ##          enrichmentScore       NES       pvalue     p.adjust    qvalues
    ## dme00053      -0.6716021 -1.878235 0.0005132592 0.0005132592 0.06267165
    ## dme04080      -0.5751475 -1.761516 0.0014763780 0.0014763780 0.09013676
    ## dme04310       0.4423567  1.622515 0.0039738859 0.0039738859 0.13853800
    ## dme00511      -0.6767705 -1.743703 0.0045383138 0.0045383138 0.13853800
    ## dme04130       0.6310228  1.698490 0.0104675506 0.0104675506 0.22378283
    ## dme00330      -0.5293307 -1.616767 0.0109962252 0.0109962252 0.22378283
    ## dme00380      -0.6379120 -1.624255 0.0187620551 0.0187620551 0.32727795
    ## dme00071      -0.5635746 -1.583765 0.0225409836 0.0225409836 0.33975038
    ## dme00830      -0.5586698 -1.562403 0.0268605646 0.0268605646 0.33975038
    ## dme00040      -0.5039618 -1.511712 0.0296996848 0.0296996848 0.33975038
    ##          rank                   leading_edge
    ## dme00053   51  tags=16%, list=0%, signal=16%
    ## dme04080 1191  tags=36%, list=9%, signal=33%
    ## dme04310 2049 tags=29%, list=16%, signal=24%
    ## dme00511 1380 tags=57%, list=11%, signal=51%
    ## dme04130 2529 tags=50%, list=20%, signal=40%
    ## dme00330 1993 tags=24%, list=16%, signal=21%
    ## dme00380 1258 tags=25%, list=10%, signal=23%
    ## dme00071 2598 tags=41%, list=21%, signal=32%
    ## dme00830 1762 tags=23%, list=14%, signal=19%
    ## dme00040  117  tags=20%, list=1%, signal=20%
    ##                                                                                                                                                core_enrichment
    ## dme00053                                                                                                                               53507/35138/34256/35139
    ## dme04080                                                 32703/37892/40955/37191/42530/34878/33248/43669/43551/36601/41639/37004/43484/41726/36368/43838/36475
    ## dme04310 39596/34011/45307/45343/32188/37713/31151/48311/42692/31310/34367/31597/40090/44317/34819/33204/43769/43319/31659/39884/34010/39751/32132/38890/31014
    ## dme00511                                                                                     37211/34203/34437/37490/32726/41913/38528/34206/35989/41457/42178
    ## dme04130                                                                                          36080/40094/42215/36015/248102/38614/40373/38541/40724/39051
    ## dme00330                                                                                   326112/40026/46717/43620/34495/41167/40028/50001/43625/326188/34256
    ## dme00380                                                                                                                               33907/43689/35190/34256
    ## dme00071                                                                               34276/34315/34313/41311/37445/41480/40059/42364/31695/37217/35213/34256
    ## dme00830                                                                                                                   41311/38598/53584/53507/35138/35139
    ## dme00040                                                                                                       53513/39305/53584/53507/38463/35586/35138/35139

Dotplot
-------

``` r
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-15-1.png)

Encrichment plot map:
---------------------

Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional modules.

``` r
 emapplot(kk2)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-16-1.png)

Category Netplot:
-----------------

The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network (helpful to see which genes are involved in enriched pathways and genes that may belong to multiple annotation categories).

``` r
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-17-1.png)

Ridgeplot
---------

Helpful to interpret up/down-regulated pathways.

``` r
ridgeplot(kk2) + labs(x = "enrichment distribution")
```

    ## Picking joint bandwidth of 0.947

![](gsea-1_files/figure-markdown_github/unnamed-chunk-18-1.png)

GSEA Plot
=========

Traditional method for visualizing GSEA result.

Params:
**Gene Set** Integer. Corresponds to gene set in the gse object. The first gene set is 1, second gene set is 2, etc. Default: 1

``` r
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
```

![](gsea-1_files/figure-markdown_github/unnamed-chunk-19-1.png)

Pathview
========

This will create a PNG and *different* PDF of the enriched KEGG pathway.

Params:
**gene.data** This is `kegg_gene_list` created above
**pathway.id** The user needs to enter this. Enriched pathways + the pathway ID are provided in the gseKEGG output table (above).
**species** Same as `organism` above in `gseKEGG`, which we defined as `kegg_organism`

``` r
library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism, kegg.native = F)
```

``` r
knitr::include_graphics("dme04130.pathview.png")
```

<img src="dme04130.pathview.png" alt="KEGG Native Enriched Pathway Plot" width="100%" />
<p class="caption">
KEGG Native Enriched Pathway Plot
</p>
