#' Create phyloseq object.
#'
#' This function takes the output of the bioinformatic pipeline into phyloseq. The output will include the run identifier for each sample in the sample_data() slot. You will want to add this column to your metadata to ensure that you can test/control for batch effects. The new metadata can then be added back to the phyloseq object as follows: YourPhyloseq@sam_data <- sample_data(YourMetadata). (Remember your sample names need to match.)
#'
#' @param SeqDataTablePath Path to the input SeqDataTable.
#' @param clustering Type of sequence clustering to be used in the phyloseq object, can be ESV, curatedESV, OTU, or curatedOTU.
#' @param assignment Taxonomic assignments from which to generate the taxa_table(). Can be "BLAST", "Idtaxa" or FALSE.
#' @param ClusterAssignment How clusters will recive taxonomic assignments. Can be "RepresentativeSequence" or between 0 and 1 to represent the proportion of reads that must share an assignment for an OTU to receive that assignment.
#' @param StandardFastqNaming Whether fastqs were named as standard convention ("[yoursamplename]_XX_L001_R1/2_001.fastq") or some other way. Default TRUE.
#' @return A phyloseq object generated from the input SeqDataTable
#' @examples
#' SeqDataTable2Phyloseq(SeqDataTablePath="23s_test_SeqDataTable.RDS", clustering="curatedOTU", assignment="BLAST", ClusterAssignment="RepresentativeSequence")
#' SeqDataTable2Phyloseq(SeqDataTablePath="~/Dropbox/BioinformaticPipeline_Env/Results/16s_multirun_test_SeqDataTable.RDS", clustering="ESV", assignment="Idtaxa")
#' @export

SeqDataTable2Phyloseq<-function(SeqDataTablePath, clustering="ESV", assignment="Idtaxa", BLASTThreshold=NULL, ClusterAssignment="RepresentativeSequence", StandardFastqNaming=TRUE){

    require(DECIPHER)
    require(phyloseq)
    require(tidyverse)

    #internal functions
        reformatSampleNames<-function(list_item) {
            string<-unlist(strsplit(list_item, "_"))
            newstring<-string[c(-1, (-length(string)+1):-length(string))]
            newstring<-paste(newstring, collapse="_")
            return(newstring)
        }
        
        HandleDuplicateSamplesInMultipleRuns<-function(otumat, StandardFastqNaming){            
            #1) pick deepest version of any duplicate samples
                StrippedSampleNames_original<-otumat %>% colnames %>%
                        strsplit(., "__") %>%
                        unlist  %>%
                        `[`(c(TRUE, FALSE))
                if (StandardFastqNaming){
                    StrippedSampleNames_original<-lapply(StrippedSampleNames_original, reformatSampleNames) %>% unlist
                }

                DuplicatedSamples<-which(table(StrippedSampleNames_original)>1) %>% names
                for (DuplicatedSample in DuplicatedSamples) {
                    DuplicateSampleRunNames<-colnames(otumat)[grep(DuplicatedSample, colnames(otumat),fixed=TRUE)]
                    SizeOfDuplicates<-otumat[,DuplicateSampleRunNames] %>% colSums
                    DuplicateToKeep<-sample(names(SizeOfDuplicates[SizeOfDuplicates==max(SizeOfDuplicates)]), 1) #keep biggest duplicate sample in case there are two equally large)
                    DuplicatesToRemove<-DuplicateSampleRunNames[DuplicateSampleRunNames!=DuplicateToKeep] #
                    otumat<-otumat[,-which(colnames(otumat)%in%DuplicatesToRemove)] #remove unwanted duplicate colums from otumat
                }
            
            # 2) strip '__RunX' from end of sample name and record 'RunX' in metadata 
                SampsAndRuns<-otumat %>% colnames %>%
                            strsplit(., "__") %>%
                            unlist 
                StrippedSampleNames<-SampsAndRuns %>%
                            `[`(c(TRUE, FALSE))
                RunIdentifier<-SampsAndRuns %>%
                            `[`(!c(TRUE, FALSE))
                
                colnames(otumat)<-StrippedSampleNames

                
                Metadata<-data.frame(row.names=colnames(otumat),RunIdentifier) %>% sample_data

            return(list(otumat, Metadata))
        }

    #check clustering valid, break if not
        validclusterings<-c("ESV", "curatedESV", "OTU", "curatedOTU")
        if (!clustering  %in%  validclusterings) { print(paste0(clustering, " clustering choice invalid please choose one of: ", validclusterings)) ; return(NULL)}
        validassignments<-c("BLAST", "Idtaxa", "None")
        if (!assignment  %in%  validassignments) { print(paste0( assignment, " assignment choice invalid please choose one of: ", validassignments)) ; return(NULL) }

    #preparatory steps

        SeqDataTable<-readRDS(SeqDataTablePath)

        #rename this column to match all others - should be fixed upstream in the pipeline later
        names(SeqDataTable)[names(SeqDataTable)=="OTUrepresentativeSequence"]<-"OTURepresentativeSequence"


        #get seq datatable and reorder columns
            clusteringcolumns<- c(  "ESV",
                                    "Sequence",
                                    "curatedESV",                      
                                    "CuratedESVRepresentativeSequence",
                                    "OTU",
                                    "OTURepresentativeSequence",       
                                    "curatedOTU",                       
                                    "CuratedOTURepresentativeSequence")
            
            #ensure standard order to clustering cols as above
            SeqDataTable<-SeqDataTable %>% relocate(clusteringcolumns)

        #recode clustering variables to work with tidyverse as arguments
            symclustering<-sym(clustering)
            RepresentativeSequenceColumnName<-paste0(toupper(substr(clustering, 1, 1)), substr(clustering, 2, nchar(clustering)), sep="", "RepresentativeSequence")

        ## get indices
            MainIndices<-which(names(SeqDataTable) %in% clusteringcolumns)
            SampleIndices<-grep("Sample_", colnames(SeqDataTable))
            if (assignment =="BLAST") {
                AssignmentIndices<-( grep("Blast_query_coverage", colnames(SeqDataTable)) +1 ) : dim(SeqDataTable)[2]
                
                #filter out blast top hits below threshold
                if (!is.null(BLASTThreshold)) {
                    ESVsBelowBlastThreshold<-SeqDataTable$Blast_percentIdentical < BLASTThreshold | SeqDataTable$Blast_query_coverage < BLASTThreshold
                    ESVsBelowBlastThresholdOrNa<-ESVsBelowBlastThreshold
                    ESVsBelowBlastThresholdOrNa[is.na(ESVsBelowBlastThresholdOrNa)] <-TRUE
                    SeqDataTable[ESVsBelowBlastThresholdOrNa,AssignmentIndices]<-NA
                }

            } else if (assignment== "Idtaxa") {
                AssignmentIndices <- ( max(SampleIndices) + 1 ) : ( grep("_confidence", ignore.case=TRUE, colnames(SeqDataTable))[1] -1 )
            }


    #create input matrices
        #otumat - with sample names reformating to match standard in metadata
            #create
                if (clustering=="ESV"){
                    otumat<-SeqDataTable %>% 
                        select("ESV", SampleIndices) %>%
                        column_to_rownames("ESV")%>% 
                        as.matrix()

                } else {        
                    otumat<-SeqDataTable %>% 
                        group_by( !!symclustering ) %>% 
                        summarise_at(colnames(SeqDataTable) [SampleIndices],sum) %>% 
                        column_to_rownames(clustering)%>% 
                        as.matrix()
                }

    #check pipeline version used
        RunAppendedToSample<-otumat %>% colnames %>% grepl( "__Run", ., fixed = TRUE) %>% all

            if (RunAppendedToSample) {
                HandledDuplicates<-HandleDuplicateSamplesInMultipleRuns(otumat, StandardFastqNaming)
                otumat<-HandledDuplicates[[1]]
                Metadata<-HandledDuplicates[[2]]

                #remove "Sample_" and added to front of sample names by pipeline and "001" added by sequencin
                if (StandardFastqNaming){
                    colnames(otumat)<-lapply(colnames(otumat), reformatSampleNames)
                    rownames(Metadata)<-lapply(rownames(Metadata), reformatSampleNames)
                } else { #strip the 'Sample_' part only 
                    colnames(otumat)<-sub("^.*?_", "", colnames(otumat))
                    rownames(Metadata)<-sub("^.*?_", "", rownames(Metadata))
                }

            } else {
                print("Warning: Your SeqDataTable was generated with a early version of the bioinformatic pipeline with suboptimal handling of multiple sequencing runs. If your data includes multiple sequencing runs consider rerunnning your bioinformatics.")
                colnames(otumat)<-lapply(colnames(otumat), reformatSampleNames)
                Metadata<-NULL
            }


        #taxmat
        if (clustering=="ESV") {
              taxmat<-SeqDataTable[,AssignmentIndices]
                taxmat$cluster<-SeqDataTable$ESV
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
        } else {

            if (ClusterAssignment=="RepresentativeSequence") { #take assignments from representative seqs
                taxmat<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,][,AssignmentIndices]
                taxmat$cluster<-SeqDataTable[[clustering]][SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE]
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
            } else { # take assignments based on most common assignments to ESVs within cluster (with threshold set as function argument)
                taxmat<-c()
                for (cluster in rownames(otumat)) { #loop through clusters
                    ClusterDataTable<-SeqDataTable[SeqDataTable[[symclustering]]==cluster,]
                    ClusterAssignments<-cbind(  ClusterDataTable[,AssignmentIndices],
                            Proportions=rowSums(ClusterDataTable[,SampleIndices])/sum(ClusterDataTable[,SampleIndices]) )
                    
                    clusterRow<-cluster
                    Ranks<-colnames(ClusterAssignments)[-length(colnames(ClusterAssignments))]
                    for (col in Ranks ) { #loop through ranks for the cluster
                        if (length(unique(ClusterAssignments[[col]]))==1 ) {
                            clusterRow <- c( clusterRow,unique(ClusterAssignments[[col]]) )
                        } else {
                            AbundancePerAssignment<-ClusterAssignments %>%
                                                        group_by(!!sym(col)) %>%
                                                        summarise(Proportions=sum(Proportions))
                            if (max(AbundancePerAssignment$Proportions) >ClusterAssignment) {
                                #colAssignment<-as.character(AbundancePerAssignment[AbundancePerAssignment$Proportions==max(AbundancePerAssignment$Proportions),1])
                                colAssignment<-as.character(AbundancePerAssignment[order(-AbundancePerAssignment$Proportions),][1,1]) #avoids problem of two equally abundant options
                                clusterRow <- c( clusterRow,colAssignment )
                            } else {
                                clusterRow<- c( clusterRow,NA )
                            }
                        }
                    }
                names(clusterRow)<-c("clusterName",Ranks)
                taxmat<-as.data.frame(rbind(taxmat, clusterRow))
                }
            rownames(taxmat)<-NULL
            taxmat<-column_to_rownames(taxmat,"clusterName" )
            taxmat<-as.matrix(taxmat)
            }
        }

        #refseqs
        if (clustering=="ESV") {
            refseqs_df<-SeqDataTable[,c("Sequence","ESV")]
            refseqs_df<- arrange(refseqs_df, ESV ) %>%
            column_to_rownames("ESV")
            refseqs <- DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        } else {
            refseqs_df<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,c("Sequence",clustering)]
            refseqs_df<- arrange(refseqs_df, clustering ) %>%
            column_to_rownames(clustering)
            refseqs <- DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        }


    #create phyloseq object depending on which matrices present
        if (assignment=="None" & is.null(Metadata)) {
            ps<-phyloseq(
                        #sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        #tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        } else if (assignment!="None" & is.null(Metadata)) {
            ps<-phyloseq(
                        #sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        } else if (assignment=="None" & !is.null(Metadata)) {
            ps<-phyloseq(
                        sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        #tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        } else if (assignment!="None" & !is.null(Metadata)) {
            ps<-phyloseq(
                        sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        }

    return(ps)
}
