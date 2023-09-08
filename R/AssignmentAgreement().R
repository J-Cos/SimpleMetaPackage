#' Get agreed assignments
#'
#' This function tests for agreement between Idtaxa and Blast assignments and returns only those assignments where algorithms agree.
#'
#' @param SDT A sequence data table with both Idtaxa and Blast assignments output by SimpleMetaPipeline.
#' @return An updated sequence data table contianing only agreeing assignments. Note that assignment quality metrics are stripped out by ths function, so any filtering based on these metrics should be done first.
#' @examples
#' AssignmentAgreement(SDT=SDT)
#' @export

AssignmentAgreement<-function(SDT){ 

    require(tidyverse)

    #get correct columns
        NumCols<-length(names(SDT))
        LastNonTaxaColumn<-which(names(SDT)=="CuratedOTURepresentativeSequence")
        NumRanks<-(NumCols-LastNonTaxaColumn-3)/3

    #extract dataframes of alternative assignments
        IdtaxaAssignments<-SDT[,(LastNonTaxaColumn+1):(LastNonTaxaColumn+NumRanks)]
        BlastAssignments<-SDT[,(NumCols-NumRanks+1):NumCols]

    #create matrix shoiwng agreements between algorithms
    MatchMatrix<-IdtaxaAssignments ==BlastAssignments

    # create new single dataframe wit only agreed assignments
        rankMatrix<-as.matrix(IdtaxaAssignments)
        Agreements <- matrix(rankMatrix[ifelse(MatchMatrix, TRUE, NA)], ncol = ncol(rankMatrix))
        Agreements<-as.data.frame(Agreements)
        colnames(Agreements) <- colnames(BlastAssignments)
    
    #replace old assignment columns with new agreed assignment columns
        newSDT<-cbind( SDT[, 1:LastNonTaxaColumn],
            Agreements)
    
    return(newSDT)
}