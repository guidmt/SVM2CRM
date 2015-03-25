\name{getSignal}
\alias{getSignal}
\alias{getSignal}
\title{
   Model the signals of each histone marks around genomic features (e.g. enhancers, not_enhancers). 
}
\description{
   This function simply model the signal of each histone marks around the features used in the input files and considering the bin.size and window size defined  during the pre-processing step.
}
\usage{
  getSignal(bedfilelist,chr,reference,win.size,bin.size,label1="enhancers")

}
\arguments{
    \item{bedfilelist}{
         test_set produced for svm model
}
    \item{chr}{
         a vector containin the list of chromsome that you want use during the analysis (e.g."chr1")
}
   \item{reference}{
         file with the reference position of the features. The genomic coordinates of positive and negative examples (e.g. enhancers, not_enhancers)
}
   \item{win.size}{
         windows size used to smooth the signal
}
   \item{bin.size}{
         original bin size used
}
  \item{label1}{
    class of reference (e.g. enhancers or not_enhancers)
}

}
\details{
    Some detailled description
}
\value{
   A data.frame with the signals where in the column there are the signals of the histone marks and in the rows the cis-regulatory elements.

}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfindbed}

\examples{
    library("SVM2CRMdata")
    library("GenomicRanges")

    setwd(system.file("data",package="SVM2CRMdata"))
    load("CD4_matrixInputSVMbin100window1000.rda")
    completeTABLE<-CD4_matrixInputSVMbin100window1000

    new.strings<-gsub(x=colnames(completeTABLE[,c(6:ncol(completeTABLE))]),pattern="CD4.",replacement="")
    new.strings<-gsub(new.strings,pattern=".norm.w100.bed",replacement="")
    colnames(completeTABLE)[c(6:ncol(completeTABLE))]<-new.strings

    #list_file<-grep(dir(),pattern=".sort.txt",value=TRUE)
    #train_positive<-getSignal(list_file,chr="chr1",reference="p300.distal.fromTSS.txt",win.size=500,bin.size=100,label1="enhancers")
    #train_negative<-getSignal(list_file,chr="chr1",reference="random.region.hg18.nop300.txt",win.size=500,bin.size=100,label1="not_enhancers")
    setwd(system.file("data",package="SVM2CRMdata"))
    load("train_positive.rda")
    load("train_negative.rda")
    training_set<-rbind(train_positive,train_negative)
    colnames(training_set)[c(5:ncol(training_set))]<-gsub(x=gsub(x=colnames(training_set[,c(5:ncol(training_set))]),pattern="sort.txt.",replacement=""),pattern="CD4.",replacement="")

}
\keyword{training-set, positive examples, negative examples,input}
