\name{createSVMinput}
\alias{createSVMinput}
\alias{createSVMinput}
\title{
     Take the output of findFeatureOverlap and then create a positive and negative set of cis-regulatory elements
}
\description{
     Take the output of findFeatureOverlap and then create a positive and negative set of cis-regulatory elements
}
\usage{
     createSVMinput(inputpos,inputfull,label1,label2)
}
\arguments{
    \item{inputpos}{
           output of findFeatureOverlap (see documentation findFeatureOverlap)
    }
    \item{inputfull}{
     output of cisREfindbed (see documentation cisREfinbed)
   }
    \item{label1}{
          a string to define the first class (e.g. enhancers)
    }
    \item{label2}{
         a string to define the second class (e.g. not_enhancers)

}
}
\details{
    Some detailled description
}
\value{
    A data.frame with the signals of the histone modifications for positive (e.g. enhancers) and negative (e.g. not_enhancers) examples.
}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{findFeatureOverlap, cisREfindbed}

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
 
    setwd(system.file("extdata",package="SVM2CRMdata"))
    data_level2 <- read.table(file = "GSM393946.distal.p300fromTSS.txt",sep = "\t", stringsAsFactors = FALSE)
    data_level2<-data_level2[data_level2[,1]=="chr1",]

    DB <- data_level2[, c(1:3)]
    colnames(DB)<-c("chromosome","start","end")

    label <- "p300"

    table.final.overlap<-findFeatureOverlap(query=completeTABLE,subject=DB)
    data_enhancer_svm<-createSVMinput(inputpos=table.final.overlap,inputfull=completeTABLE,label1="enhancers",label2="not_enhancers")
}

\keyword{inputSVM,positive examples, negativexamples, overlap}
