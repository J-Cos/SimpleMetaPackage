# SimpleMetaPackage [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7990341.svg)](https://doi.org/10.5281/zenodo.7990341) [![codecov](https://codecov.io/gh/J-Cos/SimpleMetaPackage/graph/badge.svg?token=5U5WYFJ96R)](https://codecov.io/gh/J-Cos/SimpleMetaPackage)


An R package to support the conversion of outputs from SimpleMetaPipeline (https://github.com/J-Cos/SimpleMetaPipeline) to phyloseq objects and support multi-algorithm agreement tests.


## How to use:
#### Install and load the package
    devtools::install_github("J-Cos/SimpleMetaPackage")
    library(SimpleMetaPackage)

#### Convert your pipeline output to phyloseq object
    ps<-SeqDataTable2Phyloseq(SeqDataTablePath= [Path to your SeqDataTable.RDS here] )

#### Convert your pipeline output to phyloseq object requiring agreement between clusters and taxonomic assignments
If you want to specify a given threshold of agreement between clustering and assignment algorithms this can be done with the 'ClusterAssignment' argument. For example the below requires 85% of reads in a cluster to recieve the same assignment in order for that cluster to receive the assignment.
    
    ps<-SeqDataTable2Phyloseq(SeqDataTablePath= [Path to your SeqDataTable.RDS here], clustering="OTU", ClusterAssignment=0.85)

#### Require agreement between Idtaxa and Blast assignments in a sequence data table
    data(SDT) #load a small example sequence data table
    newSDT<-AssignmentAgreement(SDT)

#### Quick example multiple rarefying for richness estimation. 
In real analysis you would not use such a shallow depth or small number of replicates.

    data(ps) #load a small example phyloseq object
    GetMultiRarefiedRichnessEstimates(ps, RarefyDepth=100, NumReplicates=5)
