\name{createBed}
\alias{createBed}
\alias{createBed}
\title{
     Create bed file of predictions using svm.
}
\description{
 With this function the user can create a bed files with the regions of predicted cis-regulatory elements. Internal function.

}
\usage{
  createBed(test_set,label1,pred,outputfile)
}
\arguments{
    \item{test_set}{
     test_set produced for svm model
}
    \item{label1}{
     class1   
}
   \item{pred}{
     p$prediction object
}
   \item{outputfile}{
     name of bed file
}

}
\details{
    Some detailled description
}
\value{
A bed files with the genomic coordinates of the cis-regulatory elements predicted using SVM2CRM.
}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfindbed}

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
    #the colnames of the training set should be the same of data_enhancer_svm
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


    results<-featSelectionWithKmeans(dftotann,5)

    resultsFS<-results[[7]]

    resultsFSfilter<-resultsFS[which(resultsFS[,2]>median(resultsFS[,2])),]

    resultsFSfilterICRR<-resultsFSfilter[which(resultsFSfilter[,3]<0.50),]

    listHM<-resultsFSfilterICRR[,1]
    listHM<-gsub(gsub(listHM,pattern="_.",replacement=""),pattern="CD4.",replacement="")

    selectFeature<-grep(x=colnames(training_set[,c(6:ncol(training_set))]),pattern=paste(listHM,collapse="|"),value=TRUE)

    colSelect<-c("chromosome","start","end","label",selectFeature)
    training_set<-training_set[,colSelect]

    vecS <- c(2:length(listHM))
    typeSVM <- c(0, 6, 7)[1]
    costV <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)[6]
    wlabel <- c("not_enhancer", "enhancer")
    infofile<-data.frame(a=c(paste(listHM,"signal",sep=".")))
    infofile[,1]<-gsub(gsub(x=infofile[,1],pattern="CD4.",replacement=""),pattern=".sort.bed",replacement="")

    tuningTAB <- tuningParametersCombROC(training_set = training_set, typeSVM = typeSVM, costV = costV,different.weight="TRUE", vecS = vecS[1],pcClass=100,ncClass=400,infofile)

    tuningTABfilter<-tuningTAB[tuningTAB$fscore<0.95,]
    #row_max_fscore<-which.max(tuningTABfilter[tuningTABfilter$nHM >2,"fscore"])
    row_max_fscore<-which.max(tuningTABfilter[,"fscore"])
    listHM_prediction<-gsub(tuningTABfilter[row_max_fscore,4],pattern="//",replacement="|")

    columnPR<-grep(colnames(training_set),pattern=paste(listHM_prediction,collapse="|"),value=TRUE)

    predictionGW(training_set=training_set,data_enhancer_svm=data_enhancer_svm, listHM=columnPR,pcClass.string="enhancers",nClass.string="not_enhancers",pcClass=100,ncClas=400,cost=100,type=0,"prediction_enhancers_CD4_results_cost=100_type=0")

}
\keyword{create bed, save output, results}

