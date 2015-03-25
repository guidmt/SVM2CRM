\name{findFeatureOverlap}
\alias{findFeatureOverlap}
\alias{findFeatureOverlap}
\title{
    Find the overlap of genomic regions between the output of cisREfind and a database of validated
    cis-regulatory elements.
}
\description{
    Find the overlap of genomic regions between the output of cisREfind and a database of validated
    cis-regulatory elements.
}
\usage{
    findFeatureOverlap(query, subject, select)
}
\arguments{
    \item{query}{
         database of validatated cis-regulatory elements (e.g. p300 binding sites), GRanges class.
    }
    \item{subject}{
         output cisREfindbed function, GRanges class.
    }
    \item{select}{
     parameters of findOverlaps
}
}
\details{
    See documentation of GenomicRanges for details.
}
\value{
    A data frame with the genomic regions that overlap 
}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfind, findOverlaps, GenomicRanges}

\examples{
    
    library("SVM2CRMdata")
    library("GenomicRanges")

    setwd(system.file("data",package="SVM2CRMdata"))
    load("CD4_matrixInputSVMbin100window1000.rda")    
    completeTABLE<-CD4_matrixInputSVMbin100window1000
    new.strings<-gsub(x=colnames(completeTABLE[,c(6:ncol(completeTABLE))]),pattern="CD4.",replacement="")
    new.strings<-gsub(new.strings,pattern=".norm.w100.bed",replacement="")
    colnames(completeTABLE)[c(6:ncol(completeTABLE))]<-new.strings
    
    setwd(system.file("extdata",package="SVM2CRMdata"))
    data_level2 <- read.table(file = "GSM393946.distal.p300fromTSS.txt",sep = "\t", stringsAsFactors = FALSE)
    data_level2<-data_level2[data_level2[,1]=="chr1",]

    DB <- data_level2[, c(1:3)]
    colnames(DB)<-c("chromosome","start","end")

    label <- "p300"

    table.final.overlap<-findFeatureOverlap(query=completeTABLE,subject=DB)
}
\keyword{overlap, positive examples, negative examples}
