\name{smoothInputFS}
\alias{smoothInputFS}
\alias{smoothInputFS}
\title{
Smooth the signals of the histone marks to prepare the input for feature selection analysis
}
\description{
Give the matrix obtained using getSignal this functions smooth the signals of each histone marks using a particular window (if bin=100).To size of smooth is  bin*k (e.g. a parameter k equal to 2 means thatthe signal is smooth every 200bp).
}
\usage{
   
smoothInputFS(input_ann,k,listcolnames)

}
\arguments{
    \item{input_ann}{
         the data.frame with the training set
}
      \item{k}{
         the size of smooth in bp 
}

      \item{listcolnames}{
         the names of column in which perform the smoothing. A vector with the list of histone marks. 
}

}
\details{
    The smoothing is perfomed using the median
}
\value{
   A data.frame with the smoothed signals of histone marks
}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfindbed, featSelectionWithKmeans, tuningParametersCombROC}

\examples{

   library("GenomicRanges")
   library("SVM2CRMdata")

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


   setwd(system.file("extdata", package = "SVM2CRMdata"))
   data_level2 <- read.table(file = "GSM393946.distal.p300fromTSS.txt",sep = "\t", stringsAsFactors = FALSE)
   data_level2<-data_level2[data_level2[,1]=="chr1",]

   DB <- data_level2[, c(1:3)]
   colnames(DB)<-c("chromosome","start","end")

   label <- "p300"

   table.final.overlap<-findFeatureOverlap(query=completeTABLE,subject=DB,select="all")

   data_enhancer_svm<-createSVMinput(inputpos=table.final.overlap,inputfull=completeTABLE,label1="enhancers",label2="not_enhancers")
   colnames(data_enhancer_svm)[c(5:ncol(data_enhancer_svm))]<-gsub(gsub(x=colnames(data_enhancer_svm[,c(5:ncol(data_enhancer_svm))]),pattern="CD4.",replacement=""),pattern=".norm.w100.bed",replacement="")

   listcolnames<-c("H2AK5ac","H2AK9ac","H3K23ac","H3K27ac","H3K27me3","H3K4me1","H3K4me3")

   dftotann<-smoothInputFS(train_positive[,c(6:ncol(train_positive))],listcolnames,k=20)

}
\keyword{features selection, smooth signal,histone marks, reduction of complexity}

