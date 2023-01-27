# BioinformaticsPackage
A trial bioinformatics package.


## How to use:
### install and load the package
    devtools::install_github("J-Cos/BioinformaticsPackage")
    library(BioinformaticsPackage)

### convert your data to phyloseq object

    ps<-SeqDataTable2Phyloseq(SeqDataTablePath= [YOUR PATH HERE] )

### if you don't have any data you can load some example data using
    data(ps)

### Here is a very small and quick to run example of multiple rarefying for richness estimation. 
### In real analysis you would use such a shallow depth or small number of replicates
    GetMultiRarefiedRichnessEstimates(ps, RarefyDepth=100, NumReplicates=5)
