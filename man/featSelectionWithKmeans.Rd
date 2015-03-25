\name{featSelectionWithKmeans}
\alias{featSelectionWithKmeans}
\alias{featSelectionWithKmeans}
\title{
  This function select the most meaningful variables in a matrix of ChIP-seq data using k-means and ICRR.
}
\description{
In the research of promoters and enhancers the users start with large dataset of histone modification. The number of histone marks is a variable that influence the prediction of enhancers. 
Several papers investigated what is the best combination of histone marks to predict enhancers. There are not a consensus about the optimal number of histone mark in the prediction of cis-regulatory elements. Moreover, in the genome there are a number of histone marks that is greaters of 50. The biological roles of all of these is not already clear. However from the computational aspect the admit of a huge number of variables (histone marks) can insert redundant infomation (signal of histone marks) during the step of prediction. Therefore here we introduced a step of features selection. The step of variablesselection consist of two step. In the first one the function perfom a smoothing of the signals of each histone mark. For example, if the data were binned every 100bp, a window size of smoothing equal to 2 means that the signal of the histone mark is smoothed every 200 bp. In the second step, the function using a k number of clusters define by the user gather the signal of each histone mark and estimate the mean of signals inside each group. Finally is estimated the median of all means above computed and all histone marks with a signal less than of this are filtered out.
We introduced  also the concept of index coverage of the regulatory regions (ICRR). In the first step the function measure the total coverage of the cis-regulatory elements in the list of positive class and considering w, the windows size used to model the signal of the histone marks around a particular features (positive class, negative class). Then is inspect in those enhancers the signals (S) of a particular histone marks is less or greater than the global mean (M) of the signal (S). Next the function compute the coverage only for these two fractions of enhancers and estimate the ICRR. This values assume values from 0 to 1. A  value close to 0 means that the difference between the coverage of the two classes of enhancers is little, in contrast, if this value is close to 1 this means that there is a diversity between the two fraction of the cis-regulatory elements. The ICRR values can be used to select the histone marks to use during the analysis. Next using tuningParametersCombROC the user can tune the bestoptimal set of svm parameters and  histone marks. 
The function featSelectionWithKmeans return a list containing the results of this analysis. The first element contain a matrix with the results of the selection analysis, the parameters used durint the analysis and the filtered histone marks. Finally this function create a report with three plots. In the first one the x-axis contain the labels of the histone marks while in y-axis reports the median of mean of each group. The second plot represent the same thing but using a scatterplot. The third page contain the same graph in the first page and a plot that contain the index of coverage between enhancers and not enhancers regions.

}
\usage{
  featSelectionWithKmeans(dftotann, nk, autoK="no",w=1000, gmax=7,outputplot="feature.selection.pdf")
}
\arguments{
    \item{dftotann}{
         a data.frame created with smoothInputFS

}
      \item{nk}{
        number of cluster to use for k-means algorithm
}
      \item{autoK}{
        logical, if autoK=yes, the function estimate automatically the number of cluster to use during k-means (autoK="no")
}
      \item{w}{
        the windows size used to build the model (default=1000bp)
}
      \item{gmax}{
       the parameter gmax used during the analysis (gmax=7)
}
      \item{outputplot}{
       the name of the report pdf file (default="feature.selection.pdf")
}

}
\details{
input: the results of smoothinputFS. This function clusterized every column (histone marks) using a K number of cluster user defined. Next the mean of the signals in each group are estimated.
The function featSelectionWithKmeans return a list containing the results of this analysis. The first element contain a matrix with the results of the selection analysis, the parameters used durint the analysis and the filtered histone marks. Finally this function create a report with three plots. In the first one the x-axis contain the labels of the histone marks while in y-axis reports the median of mean of each group. The second plot represent the same thing but using a scatterplot. The third page contain the same graph in the first page and a plot that contain the index of coverage between enhancers and not enhancers regions.
Definition of k-clusters:The user can set manually the parameter k otherwise featSelectionWithKmeans implement automatically the investigation of k using a Bayesian Information Criterion for EM initialized. 
}

\value{
        - a list where the first element is the matrix of output of feature selection analysis
        - mom: the vertical mean of the matrix in the element 1 
        - nk: number of clusters used during the clustering
        - gmax: the parameter gmax used during the analysis
        - which histone marks have a the signal that is greater than the median global signal
        - which histone marks have a signal where the histone are greater than the mean of the global signal. 
        - a data.frame where the first column contain the results of feature selection analysis with k-means 
          and in the second column the results of ICRR (index coverage regulatory regions).

}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfindbed, smoothinputFS, tuningParametersCombROC}

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
    load("train_positive.rda")
    load("train_negative.rda")

    training_set<-rbind(train_positive,train_negative)

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

}
\keyword{featureselection, k-means, relevant variables, ICRR}


