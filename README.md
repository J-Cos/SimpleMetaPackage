# BioinformaticsPackage ![image](https://github.com/J-Cos/BioinformaticsPackage/assets/72385600/1f3a9aef-3408-4b6a-b0f8-a88e75818e26)

An R package to support the conversion of outputs from SimpleMetaPipeline (https://github.com/J-Cos/SimpleMetaPipeline) to phyloseq objects and support multi-algorithm agreement tests.


## How to use:
#### Install and load the package
    devtools::install_github("J-Cos/BioinformaticsPackage")
    library(BioinformaticsPackage)

#### Convert your data to phyloseq object
    ps<-SeqDataTable2Phyloseq(SeqDataTablePath= [YOUR PATH HERE] )

#### If you don't have any data you can load some example data
    data(ps)

#### Here is a very small and quick to run example of multiple rarefying for richness estimation. In real analysis you would not use such a shallow depth or small number of replicates
    GetMultiRarefiedRichnessEstimates(ps, RarefyDepth=100, NumReplicates=5)
