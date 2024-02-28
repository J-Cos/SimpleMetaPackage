#' Create phyloseq object.
#'
#' This function takes the output of the bioinformatic pipeline into phyloseq. The output will include the run identifier for each sample in the sample_data() slot. You will want to add this column to your metadata to ensure that you can test/control for batch effects. The new metadata can then be added back to the phyloseq object as follows: YourPhyloseq@sam_data <- sample_data(YourMetadata). (Remember your sample names need to match.)
#'
#' @param SeqDataTablePath Path to the input SeqDataTable.RDS.
#' @param clustering Type of sequence clustering to be used in the phyloseq object, can be "ESV", "curatedESV", "OTU", or "curatedOTU". Default = "ESV"
#' @param assignment Taxonomic assignments from which to generate the taxa_table(). Can be "BLAST", "Idtaxa" or FALSE (no taxa_table created). Default = "Idtaxa".
#' @param BLASTThreshold Whether to filter out BLAST assignments below some threshold percent idential and query coverage, value can range between 0 and 100. If NULL then no assignments are filtered. Default = NULL.
#' @param ClusterAssignment How clusters will recive taxonomic assignments. Can be "RepresentativeSequence" or between 0 and 1 to represent the proportion of reads that must share an assignment for an OTU to receive that assignment. Default = "RepresentativeSequence".
#' @param StandardFastqNaming Whether fastqs were named as standard convention ("yoursamplename_XX_L001_R1/2_001.fastq") or some other way. Default = TRUE, you are unlikely to need to chnage this.
#' @return A phyloseq object generated from the input SeqDataTable
#' @examples
#' SeqDataTable2Phyloseq(SeqDataTablePath="YOURPATH/YOURDATANAME_SeqDataTable.RDS", clustering="curatedOTU", assignment="BLAST", ClusterAssignment="RepresentativeSequence")
#' SeqDataTable2Phyloseq(SeqDataTablePath="YOURPATH/YOURDATANAME_SeqDataTable.RDS", clustering="ESV", assignment="Idtaxa")
#' @export

SeqDataTable2Phyloseq<-function(SeqDataTablePath, clustering="ESV", assignment="Idtaxa", BLASTThreshold=NULL, ClusterAssignment="RepresentativeSequence", StandardFastqNaming=TRUE){


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
                        unlist
                
                if (length(StrippedSampleNames_original)==length(colnames(otumat))*2) { # if runs were appended length is doubled
                StrippedSampleNames_original<- StrippedSampleNames_original %>%    `[`(c(TRUE, FALSE))
                }

                        
                if (StandardFastqNaming){
                    StrippedSampleNames_original<-lapply(StrippedSampleNames_original, reformatSampleNames) %>% unlist %>% paste0("_", . ,"_")  #surrounding underscores added to ensure that nested names aren't deduplicated incorrectly
                } else if (!StandardFastqNaming) {
                   StrippedSampleNames_original <-StrippedSampleNames_original %>% paste0( . ,"_")  #  #only trailing underscore added as no name reformatting for non-standard named fastqs - thus they start with "Sample_"
                }

                DuplicatedSamples<-which(table(StrippedSampleNames_original)>1) %>% names(.)
                for (DuplicatedSample in DuplicatedSamples) {
                    DuplicateSampleRunNames<-colnames(otumat)[grep(DuplicatedSample, colnames(otumat),fixed=TRUE)]
                    SizeOfDuplicates<-otumat[,DuplicateSampleRunNames] %>% colSums
                    DuplicateToKeep<-sample(names(SizeOfDuplicates[SizeOfDuplicates==max(SizeOfDuplicates)]), 1) #keep biggest duplicate sample in case there are two equally large)
                    DuplicatesToRemove<-DuplicateSampleRunNames[DuplicateSampleRunNames!=DuplicateToKeep] #
                    otumat<-otumat[,-which(colnames(otumat)%in%DuplicatesToRemove)] #remove unwanted duplicate colums from otumat
                }
            
            # 2) strip '__RunX' from end of sample name and record 'RunX' in metadata 
                SampsAndRuns<-otumat %>% 
                            colnames(.) %>%
                            strsplit(., "__") %>%
                            unlist(.)
                if (length(SampsAndRuns)==length(colnames(otumat))*2) { # if runs were appended length is doubled
                    StrippedSampleNames<-SampsAndRuns %>%
                                `[`(c(TRUE, FALSE))
                    RunIdentifier<-SampsAndRuns %>%
                                `[`(!c(TRUE, FALSE))
                    
                    colnames(otumat)<-StrippedSampleNames

                    
                    Metadata<-data.frame(row.names=colnames(otumat),RunIdentifier) %>% phyloseq::sample_data(.)
                } else {Metadata<-NULL}
            return(list(otumat, Metadata))
        }

    #check clustering valid, break if not
        validclusterings<-c("ESV", "curatedESV", "OTU", "curatedOTU")
        if (!clustering  %in%  validclusterings) { stop(paste0(clustering, " clustering choice invalid please choose one of: ", paste(validclusterings)))}
        validassignments<-c("BLAST", "Idtaxa", "None")
        if (!assignment  %in%  validassignments) { stop(paste0( assignment, " assignment choice invalid please choose one of: ", paste(validassignments)))}
    #check ClusterAssignment within bounds if numeric
        if (!ClusterAssignment=="RepresentativeSequence" & ( ClusterAssignment<0 | ClusterAssignment>1) ) {stop ("ClusterAssignment selection invalid")}

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
            SeqDataTable<-SeqDataTable %>% dplyr::relocate(., dplyr::all_of(clusteringcolumns))

        #recode clustering variables to work with tidyverse as arguments
            symclustering<-dplyr::sym(clustering)
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
                        dplyr::select(., "ESV", dplyr::all_of(SampleIndices)) %>%
                        tibble::column_to_rownames(., "ESV")%>% 
                        as.matrix(.)

                } else {        
                    otumat<-SeqDataTable %>% 
                        dplyr::group_by(.,  !!symclustering ) %>% 
                        dplyr::summarise_at(., colnames(SeqDataTable) [SampleIndices],sum) %>% 
                        tibble::column_to_rownames(., clustering)%>% 
                        as.matrix()
                }

    #check pipeline version used
        RunAppendedToSample<-otumat %>% colnames(.) %>% grepl("__Run", ., fixed = TRUE) %>% all(.)

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
                warning("Warning: Your SeqDataTable was generated with an early version of the bioinformatic pipeline with suboptimal handling of multiple sequencing runs. If your data includes multiple sequencing runs consider rerunnning your bioinformatics.")
                HandledDuplicates<-HandleDuplicateSamplesInMultipleRuns(otumat, StandardFastqNaming)
                otumat<-HandledDuplicates[[1]]

                #remove "Sample_" and added to front of sample names by pipeline and "001" added by sequencin
                if (StandardFastqNaming){
                #colnames(otumat)<-lapply(colnames(otumat), reformatSampleNames) # hashed out as it can never be called as the combination of 
                                                                                # standardfastq naming and deprecated sample names throws an error in the 
                                                                                # HandleDuplicateSamplesInMultipleRuns function.
                } else { #strip the 'Sample_' part only 
                    colnames(otumat)<-sub("^.*?_", "", colnames(otumat))
                }
                Metadata<-NULL
            }


        #taxmat
        if (assignment != "None") {
            if (clustering=="ESV") {
                taxmat<-SeqDataTable[,AssignmentIndices]
                    taxmat$cluster<-SeqDataTable$ESV
                    taxmat<-dplyr::arrange(taxmat, cluster) %>%
                        tibble::column_to_rownames(., "cluster") %>%
                        as.matrix()
            } else {

                if (ClusterAssignment=="RepresentativeSequence") { #take assignments from representative seqs
                    taxmat<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,][,AssignmentIndices]
                    taxmat$cluster<-SeqDataTable[[clustering]][SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE]
                    taxmat<-dplyr::arrange(taxmat, cluster) %>%
                        tibble::column_to_rownames(., "cluster") %>%
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
                                                            dplyr::group_by(., !!dplyr::sym(col)) %>%
                                                            dplyr::summarise(., Proportions=sum(Proportions))
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
                taxmat<-tibble::column_to_rownames(taxmat,"clusterName" )
                taxmat<-as.matrix(taxmat)
                }
            }
        }

        #refseqs
        if (clustering=="ESV") {
            refseqs_df<-SeqDataTable[,c("Sequence","ESV")]
            refseqs_df<- dplyr::arrange(refseqs_df, ESV ) %>%
                tibble::column_to_rownames(., "ESV")
            refseqs <- Biostrings::DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        } else {
            refseqs_df<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,c("Sequence",clustering)]
            refseqs_df<- dplyr::arrange(refseqs_df, clustering ) %>%
                tibble::column_to_rownames(., clustering)
            refseqs <- Biostrings::DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        }


    #create phyloseq object depending on which matrices present
        if (assignment=="None" & is.null(Metadata)) {
            ps<-phyloseq::phyloseq(
                        #sample_data(Metadata),
                        phyloseq::otu_table(otumat, taxa_are_rows=TRUE),
                        #tax_table(taxmat),
                        #phy_tree(),
                        phyloseq::refseq(refseqs) 
                    )
        } else if (assignment!="None" & is.null(Metadata)) {
            ps<-phyloseq::phyloseq(
                        #sample_data(Metadata),
                        phyloseq::otu_table(otumat, taxa_are_rows=TRUE),
                        phyloseq::tax_table(taxmat),
                        #phy_tree(),
                        phyloseq::refseq(refseqs) 
                    )
        } else if (assignment=="None" & !is.null(Metadata)) {
            ps<-phyloseq::phyloseq(
                        phyloseq::sample_data(Metadata),
                        phyloseq::otu_table(otumat, taxa_are_rows=TRUE),
                        #tax_table(taxmat),
                        #phy_tree(),
                        phyloseq::refseq(refseqs) 
                    )
        } else if (assignment!="None" & !is.null(Metadata)) {
            ps<-phyloseq::phyloseq(
                        phyloseq::sample_data(Metadata),
                        phyloseq::otu_table(otumat, taxa_are_rows=TRUE),
                        phyloseq::tax_table(taxmat),
                        #phy_tree(),
                        phyloseq::refseq(refseqs) 
                    )
        }

    return(ps)
}
