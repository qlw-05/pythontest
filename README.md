scAPAmod
====================

Profiling APA complexity and dynamics on a global scale from scRNA-seq to provide single-cell-level insight into APA heterogeneity and pattern of APA usages.

About
====================
scAPAmod is developed for identification of patterns (or modalities) of APA usages . Given a poly(A) site expression matrix with each row denoting a poly(A) site and each column denoting a cell. Preprocessing was first performed to calculate the PUI value of each APA gene to generate a PUI matrix, with each row being a gene and each column being a cell. Then for each gene, the pattern of APA usage was obtained based on the GMM model -- unimodal, bimodal or multimodal. Unimodal means that PUI values of the gene across all cells are from the same component. Bimodal means that the PUI values consist of two components. Multimodal means that there are more than two components.

Installing scAPAmod
=============
Mandatory
---------

* R (>3.6). [R 4.0.0](https://www.r-project.org/) is recommended.

Required R Packages
---------
* attached base packages:<br>
  grid stats4 parallel stats graphics grDevices utils datasets methods base<br>
* other attached packages:<br>
clusterProfiler_3.14.3 org.Mm.eg.db_3.10.0 extrafont_0.17 ggstatsplot_0.6.6
  scAPAmod_0.1.0 ClusterR_1.2.2 gtools_3.8.2 movAPA_0.1.0 DEXSeq_1.32.0 
  DESeq2_1.26.0 SummarizedExperiment_1.16.1 DelayedArray_0.12.3 
  BiocParallel_1.20.1 matrixStats_0.57.0
 GenomicFeatures_1.38.2 AnnotationDbi_1.48.0
 Biobase_2.46.0 ggbio_1.34.0
 BSgenome_1.54.0 rtracklayer_1.46.0
 Biostrings_2.54.0 XVector_0.26.0
 ggplot2_3.3.2 data.table_1.13.2
 RColorBrewer_1.1-2 GenomicRanges_1.38.0
 GenomeInfoDb_1.22.1 IRanges_2.20.2
 S4Vectors_0.24.4 BiocGenerics_0.32.0
 reshape2_1.4.4 dplyr_1.0.2

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/scAPAmod")
library(scAPAmod)
```

Datasets
=============
* simulated data <br>
We used simulated data to evaluate the performance of scAPAmod. For unimodal, we constructed three groups of randomly generated data following one Gaussian distribution with mean values of 0, >0, and <0, respectively. For bimodal, we constructed five sets of mixed data following two Gaussian distributions with mean values of the two components being >0 & <0, >0 & >0, < 0 & <0, >0 & ~0, <0 & ~0, respectively. For multimodal, we constructed seven sets of mixed data following three Gaussian distributions with mean values of the three components being >0 & <0 & ~0, >0 & >0 & >0, <0 & <0 & <0, >0 & >0 & ~0, <0 & <0 & ~0, >0 & >0 & <0, <0 & <0 & >0, respectively. Here, the variance of the Gaussian distribution is set as 0.1 and we set the mean value ac-cording to the single-cell mouse spermatogenesis data. For uni-modal, the mean values for the three groups are 2 (>0), -2 (<0) and 0 (~0). For bimodal, the mean values for the two components in the five groups are: (2, -2), (2, 0.5), (-2, -0.5), (2, 0) and (-2, 0). For multimodal, the mean values for the three components in the seven groups are: (2, -2, 0), (2, 1, 0.5), (-2, -1, -0.5), (2, 1, 0), (-2, -1, 0), (2, 1, -2), and (-2, -1, 2).
For the 15 groups of data with the above three patterns, we iter-ate 100 times to obtain noise-free data (1500 poly(A) sites * 100 cells).

* single-cell mouse spermatogenesis data <br>
The single-cell mouse spermatogenesis data (GSE104556) was used, which contains three cell types according to their differentiation stages: the first stage is the differentiation of cells into spermatocytes (SC, 693 cells); the second stage is the differentiation of spermatocytes into round sperm cells (RS, 1140 cells); the third stage is differentiation into mature sperm cells (ES, 209 cells). We used scAPAtrap [17] to identify and quantify poly(A) sites in each single cell, which results in 43,395 poly(A) sites in 3’ UTR (including the 3’ UTR extension region within 1000 bp down-stream) and 70,108 sites in non-3’ UTR regions.



Using scAPAmod
=============
To get started, the user is recommended to use the example dataset which comes with the packages.

Preparations
---------
* PAC data of mouse sperm cells.
```
library(movAPA, warn.conflicts = FALSE, quietly=TRUE)

data(scPACds)
head(scPACds@counts[1:2,1:5])
AAACCTGAGAGGGCTT AAACCTGAGCTTATCG AAACCTGCATACGCCG AAACCTGGTTGAGTTC
PA3443 0 0 0 0
PA3446 0 0 0 0
AAACCTGTCAACGAAA
PA3443 0
PA3446 0

head(scPACds@anno, n=2)
chr strand coord peakID ftr gene_type ftr_start ftr_end
PA3443 chr12 - 100125475 peak3443 3UTR <NA> 100125452 100125605
PA3446 chr12 - 100549890 peak3446 3UTR <NA> 100549778 100551443
gene gene_start gene_end gene_stop_codon upstream_id
PA3443 ENSMUSG00000021179 100125452 100159653 100125606 <NA>
PA3446 ENSMUSG00000021180 100549778 100725028 100551444 <NA>
upstream_start upstream_end downstream_id downstream_start
PA3443 NA NA <NA> NA
PA3446 NA NA <NA> NA
downstream_end three_UTR_length three_extend
PA3443 NA 131 NA
PA3446 NA 1554 NA

head(scPACds@colData, n=2)
group celltype tsn1 tsn2
AAACCTGAGAGGGCTT AAACCTGAGAGGGCTT SC 22.54797966 4.077467845
AAACCTGAGCTTATCG AAACCTGAGCTTATCG RS 1.138437608 -32.9317999

levels(scPACds@colData$celltype)
[1] "ES" "RS" "SC"
```

* Preprocess
```
library(scAPAmod, warn.conflicts = FALSE, quietly=TRUE)

# 3' UTR
index <- which(scPACds@anno$ftr == "3UTR")
UTR_gene <- scPACds@anno$gene[index]
UTR_chr <- scPACds@anno$chr[index]
UTR_strand <- scPACds@anno$strand[index]
UTR_coord <- scPACds@anno$coord[index]
UTR_ftr_start <- scPACds@anno$ftr_start[index]
UTR_ftr_end <- scPACds@anno$ftr_end[index]
UTR_three_UTR_length <- scPACds@anno$three_UTR_length[index]
UTR_anno <- data.frame(chr = as.character(UTR_chr), 
                       strand = as.character(UTR_strand),
                       coord =as.integer(UTR_coord),
                       gene = as.character(UTR_gene),
                       ftr_start = as.integer(UTR_ftr_start),
                       ftr_end = as.integer(UTR_ftr_end),
                       three_UTR_length = as.integer(UTR_three_UTR_length))
UTR_counts <- scPACds@counts[,index]
ct1 <- which(scPACds@colData$celltype[index] == "SC")
result1 <- extrPairPA(UTR_counts[,ct1],
                      as.character(UTR_anno$gene),UTR_anno)
8 PACs

# non 3' UTR
ct <- which(scPACds@colData$celltype == "SC")
results <- exnon3UTRPA(scPACds@counts[,ct],
                       scPACds@anno$gene, scPACds@anno, scPACds@anno$ftr,
                       gn = 1, cn = 1)
114 PACs
```

Analyses of APA dynamics
---------
* Identifying modalities in 3’ UTR.
```
mod <- getMod(result1$PUI)
mod$modalities
[1] "Multimodal" "Unimodal" "Unimodal" "Bimodal"
```

* if you want to see the modalities directly in 3’ UTR, you can use UTRmod.
```
nonmod <- getMod(results$PUI)
nonmod$modalities
[1] "Unimodal" "Unimodal" "Bimodal" "Multimodal" "Bimodal"
[6] "Bimodal" "Multimodal" "Bimodal" "Multimodal" "Unimodal"
[11] "Multimodal" "Bimodal" "Bimodal" "Unimodal" "Multimodal"
[16] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[21] "Multimodal" "Bimodal" "Bimodal" "Multimodal" "Multimodal"
[26] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[31] "Multimodal" "Multimodal" "Bimodal" "Multimodal" "Unimodal"
[36] "Unimodal" "Multimodal" "Multimodal" "Bimodal" "Multimodal"
[41] "Multimodal" "Bimodal" "Multimodal" "Multimodal" "Multimodal"
[46] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[51] "Multimodal" "Multimodal" "Bimodal" "Bimodal" "Multimodal"
[56] "Multimodal" "Bimodal" "Multimodal" "Multimodal" "Multimodal"
[61] "Bimodal" "Bimodal" "Bimodal" "Multimodal" "Multimodal"
[66] "Bimodal" "Multimodal" "Bimodal" "Multimodal" "Multimodal"
[71] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[76] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Bimodal"
[81] "Multimodal" "Bimodal"
# use chi-square test to test Bimodal
ind <- which(nonmod$modalities == "Bimodal")
bigene <- results$gene[ind]
label <- lapply(c(1:length(nonmod$results)), function(y){
la <- nonmod$results[[y]][[2]][["cluster_labels"]]
if (length(which(is.na(results$PUI[y,])==TRUE))>0) {
dat.tmp <- results$PUI[y,][-which(is.na(results$PUI[y,]))]
}else{
dat.tmp <- results$PUI[y,]
}
names(la) <- names(dat.tmp)
return(la)})
bilabel <- label[ind]
# if(length(which(is.na(bigene)))>0){
# id <- which(is.na(bigene))
# bigene <- bigene[-id]
# bilabel <- bilabel[-id]}
pval1 <- chisqtest(results$filter.data, results$gene, bigene, bilabel)
pval1 <- p.adjust(pval1, method = "BH")
# use KS test to test Bimodal
pval2 <- KStest(results$PUI, results$gene, results$ftr, bigene, bilabel)
pval2 <- p.adjust(pval2, method = "BH")
```

* if you want to see the modalities directly in non 3’ UTR, you can use nonUTRmod.

```
nonmod <- nonUTRmod(scPACds,ct,gn = 1, cn = 1)
114 PACs
nonmod$modalities
[1] "Unimodal" "Unimodal" "Bimodal" "Multimodal" "Bimodal"
[6] "Bimodal" "Multimodal" "Bimodal" "Multimodal" "Unimodal"
[11] "Multimodal" "Bimodal" "Bimodal" "Unimodal" "Multimodal"
[16] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[21] "Multimodal" "Bimodal" "Bimodal" "Multimodal" "Multimodal"
[26] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[31] "Multimodal" "Multimodal" "Bimodal" "Multimodal" "Unimodal"
[36] "Unimodal" "Multimodal" "Multimodal" "Bimodal" "Multimodal"
[41] "Multimodal" "Bimodal" "Multimodal" "Multimodal" "Multimodal"
[46] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[51] "Multimodal" "Multimodal" "Bimodal" "Bimodal" "Multimodal"
[56] "Multimodal" "Bimodal" "Multimodal" "Multimodal" "Multimodal"
[61] "Bimodal" "Bimodal" "Bimodal" "Multimodal" "Multimodal"
[66] "Bimodal" "Multimodal" "Bimodal" "Multimodal" "Multimodal"
[71] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Multimodal"
[76] "Multimodal" "Multimodal" "Multimodal" "Multimodal" "Bimodal"
[81] "Multimodal" "Bimodal"
```

* Research APA preferences<br>
There are two types of APA preferences, one is the major PA, and the other is the minor PA. Here,
the ratio value of the APA data is calculated, and the largest is extracted as the major PA, and
the smallest is the minor PA. Analyze the modalities of APA usage with diﬀerent preferences and
the distribution of APA usage modalities in diﬀerent regions.
```
# the major PA
mmod <- MAMIMod(scPACds,"SC","MajorPA")
677 PACs
mmod$modalities
[1] "Multimodal" "Multimodal" "Bimodal" "Multimodal"
# the minor PA
mimod <- MAMIMod(scPACds,"SC","MinorPA")
677 PACs
mimod$modalities
[1] "Unimodal" "Multimodal" "Bimodal" "Unimodal"
```
If you want to see the ratio value specifcally, you can use exMajorPA, or if you want to directly
model the recognition modalities, you can use getMMod which is a bit diﬀerent from PUI data
modeling.
```
# the major PA and minor PA
mresult <- exMajorPA(scPACds,"SC")
677 PACs
# modalities of majorPA
mmod <- getMMod(mresult$PAmax,"PAmax")
mmod$modalities
[1] "Multimodal" "Multimodal" "Bimodal" "Multimodal"
# modalities of minorPA
mimod <- getMMod(mresult$PAmin,"PAmin")
mimod$modalities
[1] "Unimodal" "Multimodal" "Bimodal" "Unimodal"
```

Statistics of modalities
----------------------------
* Statistics on modalities distribution of diﬀerent cell types

```
# cell type of RS
ct2 <- which(scPACds@colData$celltype[index] == "RS")
result2 <- extrPairPA(UTR_counts[,ct2],
as.character(UTR_anno$gene),UTR_anno)
4 PACs
mod2 <- getMod(result2$PUI)
mod2$modalities
[1] "Bimodal" "Bimodal"
# cell type of ES
ct3 <- which(scPACds@colData$celltype[index] == "ES")
result3 <- extrPairPA(UTR_counts[,ct3],
as.character(UTR_anno$gene),UTR_anno)
4 PACs
mod3 <- getMod(result3$PUI)
mod3$modalities
[1] "Unimodal" "Unimodal"
# set the cell type
celltype <- c(rep("SC", 3), rep("RS", 3), rep("ES", 3))
data <- data.frame(celltype)
data$modality <- c("Bimodal","Multimodal","Unimodal")
data$number <- c(table(mod$modalities),table(mod2$modalities),table(mod3$modalities))
ggplot(data, aes(x=celltype, y=number)) +
ggplot2::geom_bar(stat = "identity", position = "dodge", aes(fill=modality))
```
![image](https://github.com/qlw-05/pythontest/blob/master/pictures/Rplot1.png)
* Distribution of usage modalties at diﬀerent stages of diﬀerentiation
```
library(ggplot2)
library(ggstatsplot)
#extrafont::loadfonts()
data("tUTRModalChange")
ggstatsplot::ggbarstats(data = staChange, x = condition, y = celltype,
title = "Exchange of modalities from different cell type",
ylab = "% of total common gene of PA",
# ggstatsplot.layer = FALSE,
sampling.plan = "jointMulti",
ggtheme = hrbrthemes::theme_ipsum_pub(),
legend.title = "condition", messages = F, palette = "Set2")
```
![image](https://github.com/qlw-05/pythontest/blob/master/pictures/Rplot2.png)
* Changes in the modalities of diﬀerent cell diﬀerentiation stages
```
library(extrafont)
data("tUTRModalChangeDetail")
ggstatsplot::ggbarstats(data = detChange, x = condition, y = celltype,
title = "Exchange of modalities from different cell type",
ylab = "% of total common gene of PA",
# ggstatsplot.layer = FALSE,
sampling.plan = "jointMulti",
ggtheme = hrbrthemes::theme_ipsum_pub(),
legend.title = "condition", messages = F, palette = "Set2")
```
![image](https://github.com/qlw-05/pythontest/blob/master/pictures/Rplot3.png)
* Visualize the distribution of PA expression according to the components
```
# heatmap
library(grid)
library(org.Mm.eg.db)
PUI <- result1$PUI
tUTR.pair.cd.tmp <- result1$filter.data
tUTR.gene <- result1$gene
tUTR.gene.pui <- rownames(PUI)
which(mod$modalities == "Bimodal")
integer(0)
id2 <- tUTR.gene.pui[3]
genename <- select(org.Mm.eg.db, keys = id2,
                   columns = c("SYMBOL","ENTREZID","GENENAME"),
                   keytype = "ENSEMBL")
tUTR.gene <- select(org.Mm.eg.db, keys = tUTR.gene,
                    columns = c("SYMBOL","ENTREZID","GENENAME"),
                    keytype = "ENSEMBL")
plotGenePACount3(org.Mm.eg.db, genename$ENTREZID, tUTR.pair.cd.tmp,
                 tUTR.gene$ENTREZID,label[[3]])

There was a problem when running diffusion map. Trying PCA instead...
The standard deviations of PC1: 1.230856
```
![image](https://github.com/qlw-05/pythontest/blob/master/pictures/Rplot4.png)
* GO analysis
```
library(clusterProfiler)
scego <- enrichGO(OrgDb="org.Mm.eg.db", gene = rownames(result1$PUI),
                  keyType = "ENSEMBL", ont = "ALL", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable= TRUE)
barplot(scego,showCategory = 10)
```
![image](https://github.com/qlw-05/pythontest/blob/master/pictures/Rplot5.png)

Session Information
---------------
The session information records the versions of all the packages used in the generation of the present
document.
```
R version 4.0.0 (2020-04-24)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] movAPA_0.1.0                DEXSeq_1.36.0               BiocParallel_1.24.1         DESeq2_1.30.1              
 [5] SummarizedExperiment_1.20.0 MatrixGenerics_1.2.1        matrixStats_0.61.0          GenomicFeatures_1.42.3     
 [9] AnnotationDbi_1.52.0        Biobase_2.50.0              ggbio_1.38.0                BSgenome_1.58.0            
[13] rtracklayer_1.49.5          Biostrings_2.58.0           XVector_0.30.0              ggplot2_3.3.5              
[17] data.table_1.14.2           RColorBrewer_1.1-2          GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
[21] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1         reshape2_1.4.4             
[25] dplyr_1.0.7                

loaded via a namespace (and not attached):
  [1] utf8_1.2.2               tidyselect_1.1.1         RSQLite_2.2.8            htmlwidgets_1.5.4       
  [5] grid_4.0.0               gmp_0.6-2                devtools_2.4.2           munsell_0.5.0           
  [9] statmod_1.4.36           withr_2.4.2              colorspace_2.0-2         OrganismDbi_1.32.0      
 [13] knitr_1.36               rstudioapi_0.13          Rttf2pt1_1.3.9           GenomeInfoDbData_1.2.4  
 [17] hwriter_1.3.2            bit64_4.0.5              datawizard_0.2.1         rprojroot_2.0.2         
 [21] vctrs_0.3.8              generics_0.1.1           xfun_0.26                biovizBase_1.38.0       
 [25] BiocFileCache_1.14.0     BWStest_0.2.2            R6_2.5.1                 locfit_1.5-9.4          
 [29] AnnotationFilter_1.14.0  bitops_1.0-7             cachem_1.0.6             reshape_0.8.8           
 [33] DelayedArray_0.16.3      assertthat_0.2.1         scales_1.1.1             nnet_7.3-16             
 [37] gtable_0.3.0             multcompView_0.1-8       processx_3.5.2           ensembldb_2.14.1        
 [41] rlang_0.4.12             zeallot_0.1.0            genefilter_1.72.1        systemfonts_1.0.3       
 [45] PMCMRplus_1.9.2          splines_4.0.0            extrafontdb_1.0          lazyeval_0.2.2          
 [49] dichromat_2.0-0          checkmate_2.0.0          BiocManager_1.30.16      backports_1.3.0         
 [53] Hmisc_4.6-0              RBGL_1.66.0              extrafont_0.17           usethis_2.1.3           
 [57] tools_4.0.0              ellipsis_0.3.2           WRS2_1.1-3               sessioninfo_1.1.1       
 [61] Rcpp_1.0.7               plyr_1.8.6               base64enc_0.1-3          progress_1.2.2          
 [65] zlibbioc_1.36.0          purrr_0.3.4              RCurl_1.98-1.5           ps_1.6.0                
 [69] prettyunits_1.1.1        rpart_4.1-15             openssl_1.4.5            correlation_0.7.1       
 [73] cluster_2.1.2            hrbrthemes_0.8.0         fs_1.5.0                 magrittr_2.0.1          
 [77] mvtnorm_1.1-3            ProtGenerics_1.22.0      pkgload_1.2.3            hms_1.1.1               
 [81] patchwork_1.1.1          evaluate_0.14            xtable_1.8-4             XML_3.99-0.8            
 [85] jpeg_0.1-9               gridExtra_2.3            testthat_3.1.0           compiler_4.0.0          
 [89] biomaRt_2.46.3           tibble_3.1.5             ggstatsplot_0.9.0        crayon_1.4.2            
 [93] htmltools_0.5.2          mc2d_0.1-21              Formula_1.2-4            geneplotter_1.68.0      
 [97] DBI_1.1.1                SuppDists_1.1-9.5        kSamples_1.2-9           dbplyr_2.1.1            
[101] MASS_7.3-54              rappdirs_0.3.3           Matrix_1.3-4             cli_3.1.0               
[105] insight_0.14.5           pkgconfig_2.0.3          GenomicAlignments_1.26.0 statsExpressions_1.2.0  
[109] foreign_0.8-81           xml2_1.3.2               paletteer_1.4.0          annotate_1.68.0         
[113] stringr_1.4.0            VariantAnnotation_1.36.0 callr_3.7.0              digest_0.6.28           
[117] parameters_0.15.0        graph_1.68.0             rmarkdown_2.11           htmlTable_2.3.0         
[121] gdtools_0.2.3            curl_4.3.2               Rsamtools_2.6.0          lifecycle_1.0.1         
[125] desc_1.4.0               askpass_1.1              fansi_0.5.0              pillar_1.6.4            
[129] lattice_0.20-45          GGally_2.1.2             fastmap_1.1.0            httr_1.4.2              
[133] pkgbuild_1.2.0           survival_3.2-13          glue_1.4.2               remotes_2.4.1           
[137] bayestestR_0.11.5        png_0.1-7                bit_4.0.4                stringi_1.7.5           
[141] performance_0.8.0        rematch2_2.1.2           blob_1.2.2               latticeExtra_0.6-29     
[145] memoise_2.0.0            Rmpfr_0.8-7  
```


Citation
---------
If you are using scAPAmod, please cite: [scAPAmod: profiling alternative polyadenyla-tion modalities in single cells from single-cell RNA-seq data, , , .](https://github.com/BMILAB/scAPAmod)
