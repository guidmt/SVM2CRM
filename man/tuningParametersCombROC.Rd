\name{tuningParametersCombROC}
\alias{tuningParametersCombROC}
\alias{tuningParametersCombROC}
\title{
 Test different models using different kernel of SVM, values of cost functions, the number of histone marks. 
}
\description{
The function tuningParametersCombROC allow high flexibility: the user can set the type of kernel, the cost parameter of SVM, the number of histone marks. tuningParametersCombROC use performanceSVM and the output is a data.frame where for each model there are the parameters compute with performanceSVM. 
To help the user to discriminate how to discriminate the best model SVM2CRM implement several function to plot the results of performanceSVM.
}
\usage{
 tuningParametersCombROC(training_set, typeSVM, costV,different.weight=TRUE, vecS,infofile,pcClass.string="enhancers",nClass.string="not_enhancers",pcClass,ncClass) 
}
\arguments{
    \item{training_set}{
         training set (data.frame)
}
    \item{typeSVM}{
         a vector with the kinds of Kernel

}
   \item{costV}{
        a vector with a list of cost parameters 
}
   \item{different.weight}{
        the data are unbalanced (default TRUE)

}
   \item{vecS}{
       a vector with the number of histone marks (e.g. from 2 to x)
}

   \item{infofile}{
       a data.frame where in the column "a" there are the histone marks, while in the column "b" a vectors of letters.
}

\item{pcClass.string}{
      label of the first class (e.g. "enhancers")

}

\item{nClass.string}{
      label of the second class (e.g. "not_enhancers")
}

\item{pcClass}{
       number of positive class in the training set
}

   \item{ncClass}{
       number of negative class in the training set
}


}

\details{
    Some detailled description
}
\value{
   A data.frame with the values from perfomanceSVM for each trained model.

}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfindbed, perfomanceSVM, plotFscore, plotROC}

\examples{
   library("GenomicRanges")
   library("SVM2CRMdata")

   setwd(system.file("data",package="SVM2CRMdata"))
   load("CD4_matrixInputSVMbin100window1000.rda")
   completeTABLE<-CD4_matrixInputSVMbin100window1000

   new.strings<-gsub(x=colnames(completeTABLE[,c(6:ncol(completeTABLE))]),pattern="CD4.",replacement="")
   new.strings<-gsub(new.strings,pattern=".norm.w100.bed",replacement="")
   colnames(completeTABLE)[c(6:ncol(completeTABLE))]<-new.strings

   #Create a vector that contain the list of the bed files that you want use during the analysis
   #list_file<-grep(dir(),pattern=".sort.txt",value=TRUE)
   #print(list_file)
   #Here we used a data.frame that contain the genomic coordinates of p300 binding sites
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
}
\keyword{tuningparameters,brute-force,cost,type,SVM,models,combination,histone-marks}

