
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `muleaData` - GMT Datasets for the [`mulea`](https://github.com/ELTEbioinformatics/mulea) Package

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/ELTEbioinformatics/muleaData)](https://github.com/ELTEbioinformatics/muleaData/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/ELTEbioinformatics/muleaData)](https://github.com/ELTEbioinformatics/muleaData/pulls)

<!-- badges: end -->

`muleaData` is an ExperimentHubData Bioconductor package providing
pre-processed gene set data for use with the
[`mulea`](https://github.com/ELTEbioinformatics/mulea) R package, a
comprehensive tool for overrepresentation and functional enrichment
analysis. `mulea` leverages ontologies (gene and protein sets) stored in
the standardized Gene Matrix Transposed (GMT) format.

We provide these *GMT* files for 27 different model organisms, ranging
from *Escherichia coli* to human. These files are compiled from publicly
available sources and include various gene and protein identifiers like
*UniProt* protein IDs, *Entrez*, *Gene Symbol*, and *Ensembl* gene IDs.
The GMT files and the scripts we applied to create them are available at
the
[GMT_files_for_mulea](https://github.com/ELTEbioinformatics/GMT_files_for_mulea)
repository. For the `muleaData` we read these *GMT* files with the
`mulea::read_gmt()` function and saved them to **.rds** files with the
standard R `saveRDS()` function.

List of species `muleaData` covers:

- *Arabidopsis thaliana*
- *Bacillus subtilis*
- *Bacteroides thetaiotaomicron VPI-5482*
- *Bifidobacterium longum*
- *Bos taurus*
- *Caenorhabditis elegans*
- *Chlamydomonas reinhardtii*
- *Danio rerio*
- *Daphnia pulex*
- *Dictyostelium discoideum*
- *Drosophila melanogaster*
- *Drosophila simulans*
- *Escherichia coli*
- *Gallus gallus*
- *Homo sapiens*
- *Macaca mulatta*
- *Mus musculus*
- *Mycobacterium tuberculosis*
- *Neurospora crassa*
- *Pan troglodytes*
- *Rattus norvegicus*
- *Saccharomyces cerevisiae*
- *Salmonella enterica subsp. enterica serovar Typhimurium str. LT2*
- *Schizosaccharomyces pombe*
- *Tetrahymena thermophila*
- *Xenopus tropicalis*
- *Zea mays*

Type, name, link and citation of the databases `muleaData` covers:

|                                     |                                                                                        |                                                                                                                                                                    |                                                                                                                                                                                                |
|-------------------------------------|:--------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| **Ontology category**               |                                   **Ontology name**                                    |                                                                  **Short description of content**                                                                  |                                                                                         **Reference**                                                                                          |
| **Gene expression**                 |                          [FlyAtlas](http://www.flyatlas.org/)                          |                                                   Tissue-specific expression data for *Drosophila melanogaster*.                                                   |                      Chintapalli,V.R. *et al.* (2007) Using FlyAtlas to identify better *Drosophila melanogaster* models of human disease. *Nat Genet*, **39**, 715–720.                       |
|                                     |                        [ModEncode](http://data.modencode.org/)                         | Functional characterization (cell line, temporal expression, tissue expression, treatment) of elements for *Caenorhabditis elegans* and *Drosophila melanogaster*. |                The Modencode Consortium *et al.* (2010) Identification of functional elements and regulatory circuits by *Drosophila* modENCODE. *Science*, **330**, 1787–1797.                |
| **Genomic location**                |                                   Chromosomal Bands                                    |                                                                Location of genes on the chromosome.                                                                |                                                       Martin,F.J. *et al.* (2023) Ensembl 2023. *Nucleic Acids Res,* **51**, D933–D941.                                                        |
|                                     |                                   Consecutive genes                                    |                                                              *n* consecutive genes on the chromosome.                                                              |                                                                                                                                                                                                |
| **miRNA regulation**                | [miRTarBase](https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php) |                                                       Experimentally validated miRNA - target interactions.                                                        |            Huang,H.-Y. et al. (2022) miRTarBase update 2022: an informative resource for experimentally validated miRNA–target interactions. Nucleic Acids Res, **50**, D222–D230.             |
| **Gene Ontology**                   |                            [GO](https://geneontology.org/)                             |                                            Gene Ontology (GO) categorizes genes into unified categories and attributes.                                            |                                      The Gene Ontology Consortium *et al.* (2023) The Gene Ontology knowledgebase in 2023. *Genetics*, **224**, iyad031.                                       |
| **Pathway**                         |                   [Pathway Commons](https://www.pathwaycommons.org/)                   |                                                       Collection of biological pathway and interaction data.                                                       |                    Rodchenkov,I. et al. (2020) Pathway Commons 2019 Update: integration, analysis and exploration of pathway data. *Nucleic Acids Res*, **48**, D489–D497.                     |
|                                     |                           [Reactome](https://reactome.org/)                            |                                                       Collection of biological pathway and interaction data.                                                       |                                             Jassal,B. *et al.* (2020) The reactome pathway knowledgebase. *Nucleic Acids Res*, **48**, D498–D503.                                              |
|                                     |                           [Signalink](http://signalink.org/)                           |                                              Interaction database focussing on pathways and interactions of pathways.                                              |                     Csabai,L. *et al.* (2022) SignaLink3: a multi-layered resource to uncover tissue-specific signaling networks. *Nucleic Acids Res*, **50**, D701–D709.                      |
|                                     |                     [Wikipathways](https://www.wikipathways.org/)                      |                                                       Collection of biological pathway and interaction data.                                                       |                                            Martens,M. *et al.* (2021) WikiPathways: connecting communities. *Nucleic Acids Res*, **49**, D613–D621.                                            |
| **Protein domain**                  |                             [PFAM](http://pfam.xfam.org/)                              |                                                                 Protein domain structure database.                                                                 |                                         Mistry,J. *et al.* (2021) Pfam: The protein families database in 2021. *Nucleic Acids Res*, **49**, D412–D419.                                         |
| **Transcription factor regulation** |                            [ATRM](http://atrm.gao-lab.org/)                            |                                            Transcription factor - target gene interactions for *Arabidopsis thaliana*.                                             | Jin,J. et al. (2015) An *Arabidopsis* transcriptional regulatory map reveals distinct functional and evolutionary features of novel transcription factors. *Mol Biol Evol*, **32**, 1767–1773. |
|                                     |                    [dorothEA](https://saezlab.github.io/dorothea/)                     |                                                Transcription factor - target gene interactions for human and mouse.                                                |             Garcia-Alonso,L. *et al.* (2019) Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Res*, **29**, 1363–1375.              |
|                                     |                      [RegulonDB](https://regulondb.ccg.unam.mx/)                       |                                          Transcription factor - target gene interactions for *Escherichia coli* bacteria.                                          |        Tierrafría,V.H. *et al.* (2022) RegulonDB 11.0: Comprehensive high-throughput datasets on transcriptional regulation in *Escherichia coli* K-12. *Microb Genom*, **8**, 000833.         |
|                                     |                             [TFLink](https://tflink.net/)                              |                              Small- and large-scale transcription factor - target gene interactions for human and 6 model organisms.                               |              Liska,O. *et al.* (2022) TFLink: an integrated gateway to access transcription factor–target gene interactions for multiple species. *Database*, **2022**, baac083.               |
|                                     |                       [TRRUST](https://www.grnpedia.org/trrust/)                       |                                                     Transcription factor - target gene interactions for human.                                                     |              Han,H. *et al.* (2018) TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. *Nucleic Acids Res*, **46**, D380–D386.              |
|                                     |                         [Yeastract](http://www.yeastract.com/)                         |                                          Transcription factor - target gene interactions for *Saccharomyces cerevisiae*.                                           |    Teixeira,M.C. *et al.* (2018) YEASTRACT: an upgraded database for the analysis of transcription regulatory networks in Saccharomyces cerevisiae. *Nucleic Acids Res*, **46**, D348–D353.    |

# Installation Instructions

Install the developmental version of `R` from
[CRAN](https://cran.r-project.org/sources.html). Then install the
developmental version of [Bioconductor](http://bioconductor.org/) and
the `ExperimentHub` library using the following code:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ExperimentHub")
BiocManager::install("muleaData")
```

# Example

This is a basic example which shows you how to use the `muleaData`:

``` r

# Calling the ExperimentHub library.
library(ExperimentHub)
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
#>     tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: AnnotationHub
#> Loading required package: BiocFileCache
#> Loading required package: dbplyr

# Downloading the metadata from ExperimentHub.
eh <- ExperimentHub()

# Creating the muleaData variable.
muleaData <- query(eh, "muleaData")

# Checking the muleaData variable.
muleaData
#> ExperimentHub with 879 records
#> # snapshotDate(): 2024-03-06
#> # $dataprovider: muleaData
#> # $species: Drosophila melanogaster, Homo sapiens, Mus musculus, Caenorhabdi...
#> # $rdataclass: data.frame
#> # additional mcols(): taxonomyid, genome, description,
#> #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#> #   rdatapath, sourceurl, sourcetype 
#> # retrieve records with, e.g., 'object[["EH8571"]]' 
#> 
#>            title                                                              
#>   EH8571 | Genomic_location_Ensembl_Arabidopsis_thaliana_10genes_EnsemblID.rds
#>   EH8572 | Genomic_location_Ensembl_Arabidopsis_thaliana_10genes_EntrezID.rds 
#>   EH8573 | Genomic_location_Ensembl_Arabidopsis_thaliana_10genes_GeneSymbol...
#>   EH8574 | Genomic_location_Ensembl_Arabidopsis_thaliana_10genes_UniprotID.rds
#>   EH8575 | Genomic_location_Ensembl_Arabidopsis_thaliana_20genes_EnsemblID.rds
#>   ...      ...                                                                
#>   EH9445 | Genomic_location_Ensembl_Zea_mays_5genes_UniprotID.rds             
#>   EH9446 | Protein_domain_PFAM_Zea_mays_EnsemblID.rds                         
#>   EH9447 | Protein_domain_PFAM_Zea_mays_EntrezID.rds                          
#>   EH9448 | Protein_domain_PFAM_Zea_mays_GeneSymbol.rds                        
#>   EH9449 | Protein_domain_PFAM_Zea_mays_UniprotID.rds

# Looking for the ExperimentalHub ID of i.e. target genes of transcription
# factors from TFLink in Caenorhabditis elegans.
mcols(muleaData) %>% 
  as.data.frame() %>% 
  dplyr::filter(species == "Caenorhabditis elegans" & 
                  sourceurl == "https://tflink.net/")
#>                                                                        title
#> EH8727  Transcription_factor_TFLink_Caenorhabditis_elegans_All_EnsemblID.rds
#> EH8728   Transcription_factor_TFLink_Caenorhabditis_elegans_All_EntrezID.rds
#> EH8729 Transcription_factor_TFLink_Caenorhabditis_elegans_All_GeneSymbol.rds
#> EH8730  Transcription_factor_TFLink_Caenorhabditis_elegans_All_UniprotID.rds
#> EH8731   Transcription_factor_TFLink_Caenorhabditis_elegans_LS_EnsemblID.rds
#> EH8732    Transcription_factor_TFLink_Caenorhabditis_elegans_LS_EntrezID.rds
#> EH8733  Transcription_factor_TFLink_Caenorhabditis_elegans_LS_GeneSymbol.rds
#> EH8734   Transcription_factor_TFLink_Caenorhabditis_elegans_LS_UniprotID.rds
#> EH8735   Transcription_factor_TFLink_Caenorhabditis_elegans_SS_EnsemblID.rds
#> EH8736    Transcription_factor_TFLink_Caenorhabditis_elegans_SS_EntrezID.rds
#> EH8737  Transcription_factor_TFLink_Caenorhabditis_elegans_SS_GeneSymbol.rds
#> EH8738   Transcription_factor_TFLink_Caenorhabditis_elegans_SS_UniprotID.rds
#>        dataprovider                species taxonomyid genome
#> EH8727    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8728    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8729    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8730    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8731    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8732    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8733    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8734    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8735    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8736    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8737    muleaData Caenorhabditis elegans       6239   <NA>
#> EH8738    muleaData Caenorhabditis elegans       6239   <NA>
#>                                                                                                                                                                                                                                                                                                                                                                         description
#> EH8727 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8728 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8729 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8730 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8731 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8732 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8733 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8734 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8735 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8736 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8737 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#> EH8738 ExperimentHubData package for the 'mulea' comprehensive overrepresentation and functional enrichment analyser R package. Here we provide ontologies (gene sets) in a data.frame for 27 different organisms, ranging from Escherichia coli to human, all acquired from publicly available data sources. Each ontology is provided with multiple gene and protein identifiers.
#>        coordinate_1_based                     maintainer rdatadateadded
#> EH8727                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8728                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8729                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8730                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8731                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8732                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8733                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8734                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8735                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8736                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8737                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#> EH8738                  1 Eszter Ari arieszter@gmail.com     2024-02-07
#>        preparerclass         tags rdataclass
#> EH8727     muleaData Arabidop.... data.frame
#> EH8728     muleaData Arabidop.... data.frame
#> EH8729     muleaData Arabidop.... data.frame
#> EH8730     muleaData Arabidop.... data.frame
#> EH8731     muleaData Arabidop.... data.frame
#> EH8732     muleaData Arabidop.... data.frame
#> EH8733     muleaData Arabidop.... data.frame
#> EH8734     muleaData Arabidop.... data.frame
#> EH8735     muleaData Arabidop.... data.frame
#> EH8736     muleaData Arabidop.... data.frame
#> EH8737     muleaData Arabidop.... data.frame
#> EH8738     muleaData Arabidop.... data.frame
#>                                                                              rdatapath
#> EH8727  muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_All_EnsemblID.rds
#> EH8728   muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_All_EntrezID.rds
#> EH8729 muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_All_GeneSymbol.rds
#> EH8730  muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_All_UniprotID.rds
#> EH8731   muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_LS_EnsemblID.rds
#> EH8732    muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_LS_EntrezID.rds
#> EH8733  muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_LS_GeneSymbol.rds
#> EH8734   muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_LS_UniprotID.rds
#> EH8735   muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_SS_EnsemblID.rds
#> EH8736    muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_SS_EntrezID.rds
#> EH8737  muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_SS_GeneSymbol.rds
#> EH8738   muleaData/Transcription_factor_TFLink_Caenorhabditis_elegans_SS_UniprotID.rds
#>                  sourceurl sourcetype
#> EH8727 https://tflink.net/   Multiple
#> EH8728 https://tflink.net/   Multiple
#> EH8729 https://tflink.net/   Multiple
#> EH8730 https://tflink.net/   Multiple
#> EH8731 https://tflink.net/   Multiple
#> EH8732 https://tflink.net/   Multiple
#> EH8733 https://tflink.net/   Multiple
#> EH8734 https://tflink.net/   Multiple
#> EH8735 https://tflink.net/   Multiple
#> EH8736 https://tflink.net/   Multiple
#> EH8737 https://tflink.net/   Multiple
#> EH8738 https://tflink.net/   Multiple

# Creating a variable for the GMT data.frame of EH8735.
# EH8735 contains small-scale measurement results, where the target genes are
# coded with Ensembl ID-s
Transcription_factor_TFLink_Caenorhabditis_elegans_SS_EnsemblID <- muleaData[["EH8735"]]
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://cloud.r-project.org
#> Bioconductor version 3.19 (BiocManager 1.30.22), R Under development (unstable)
#>   (2024-02-21 r85967)
#> Installing package(s) 'muleaData'
#> Warning: package 'muleaData' is not available for Bioconductor version '3.19'
#> 
#> A version of this package for your version of R might be available elsewhere,
#> see the ideas at
#> https://cran.r-project.org/doc/manuals/r-devel/R-admin.html#Installing-packages
#> Installation paths not writeable, unable to update packages
#>   path: /usr/local/lib/R/library
#>   packages:
#>     boot
#> loading from cache
```

# Session Info

``` r
sessionInfo()
#> R Under development (unstable) (2024-02-21 r85967)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Etc/UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ExperimentHub_2.11.1 AnnotationHub_3.11.1 BiocFileCache_2.11.1
#> [4] dbplyr_2.4.0         BiocGenerics_0.49.1 
#> 
#> loaded via a namespace (and not attached):
#>  [1] rappdirs_0.3.3          utf8_1.2.4              generics_0.1.3         
#>  [4] BiocVersion_3.19.1      RSQLite_2.3.5           digest_0.6.34          
#>  [7] magrittr_2.0.3          evaluate_0.23           fastmap_1.1.1          
#> [10] blob_1.2.4              AnnotationDbi_1.65.2    GenomeInfoDb_1.39.7    
#> [13] DBI_1.2.2               BiocManager_1.30.22     httr_1.4.7             
#> [16] purrr_1.0.2             fansi_1.0.6             Biostrings_2.71.2      
#> [19] cli_3.6.2               rlang_1.1.3             crayon_1.5.2           
#> [22] XVector_0.43.1          Biobase_2.63.0          bit64_4.0.5            
#> [25] withr_3.0.0             cachem_1.0.8            yaml_2.3.8             
#> [28] tools_4.4.0             memoise_2.0.1           dplyr_1.1.4            
#> [31] GenomeInfoDbData_1.2.11 filelock_1.0.3          curl_5.2.1             
#> [34] mime_0.12               vctrs_0.6.5             R6_2.5.1               
#> [37] png_0.1-8               stats4_4.4.0            lifecycle_1.0.4        
#> [40] zlibbioc_1.49.0         KEGGREST_1.43.0         S4Vectors_0.41.4       
#> [43] IRanges_2.37.1          bit_4.0.5               pkgconfig_2.0.3        
#> [46] pillar_1.9.0            glue_1.7.0              xfun_0.42              
#> [49] tibble_3.2.1            tidyselect_1.2.0        rstudioapi_0.15.0      
#> [52] knitr_1.45              htmltools_0.5.7         rmarkdown_2.26         
#> [55] compiler_4.4.0
```

# Citation

To cite package `muleaData` in publications use:

C. Turek, M. Olbei, T. Stirling, G. Fekete, E. Tasnadi, L. Gul, B.
Bohar, B. Papp, W. Jurkowski, E. Ari: mulea - an R package for
enrichment analysis using multiple ontologies and empirical FDR
correction. *bioRxiv* (2024),
[doi:10.1101/2024.02.28.582444](https://doi.org/10.1101/2024.02.28.582444).
