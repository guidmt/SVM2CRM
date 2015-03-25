################################################################################################
################################################################################################
###
### METHOD: findFeatureOverlap
### CLASS: -------
###
### Find overlap between between the output of cisREfind and a database that validated 
### cis-regulatory elements. 
###
################################################################################################
################################################################################################
###
### INPUT
###
###     - query: database of validatated cis-regulatory elements
###
###     - subject: output cisREfind 
###
###     - select: parameter of findOverlaps functions (default all)
###          
### OUTPUT
###
###     - Overlapping regions between query and subject 
###
################################################################################################
################################################################################################

findFeatureOverlap <- function(query, subject, select="all") {
    
    subjectGRanges <- GRanges(subject[, "chromosome"], ranges = IRanges(subject[, "start"], subject[, "end"]))
    
    queryGRanges <- GRanges(query[, "chromosome"], ranges = IRanges(query[, "start"], query[, "end"]))
    
    if (select=="all"){

    overlap <- findOverlaps(query = queryGRanges, subject = subjectGRanges, select = "all")
    
    matrix.overlap <- as.matrix(overlap)
        
    row.data <- matrix.overlap[, 1]
    
    row.DB <- matrix.overlap[, 2]
    
    table.correlation.analysis.subsetovp<-query[row.data,]    

    colnames(table.correlation.analysis.subsetovp)[c(2,4,5)]<-c("chromosome.feat","start.feat","end.feat")

    dbvalidation.select.ovp <- subject[row.DB, ]
    
    table.final.ovp <- data.frame(table.correlation.analysis.subsetovp, dbvalidation.select.ovp)
    }

    if (select=="first"){

    overlap <- findOverlaps(query = queryGRanges, subject = subjectGRanges, select = "first")
    row.data <- which(!is.na(overlap))
    table.correlation.analysis.subsetovp <- query[row.data,] 
    table.final.ovp <- table.correlation.analysis.subsetovp
    }
    return(table.final.ovp)
    
} 

