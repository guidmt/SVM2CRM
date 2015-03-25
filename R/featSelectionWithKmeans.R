################################################################################################
################################################################################################
###
### METHOD: featSelectionWithKmeans
### CLASS: ------
### 
### This function select the most significant variables in a matrix of ChIP-seq data using
### k means
### 
################################################################################################
################################################################################################
###
### INPUT
###     
###     - dftotann: a data.frame created with smoothInputFS
###     - nk: number of cluster to use for k-means algorithm
###     - autooK=logical, if autoK=yes, the function estimate automatically the number of cluster to use during k-means
###     - w: the windows size used to build the model (default=1000bp)
###     - gmax: the number max of cluster to use if w is set to "yes"
###     - vecS:a vector with the number of histone marks (e.g. from 2 to x)
###     - outputplot: the name of the report pdf file
###
### OUTPUT
###
###     - a list where the first element is the matrix of output of feature selection analysis
###     - mom: the vertical mean of the matrix in the element 1 
###     - nk: number of cluster used during the clustering
###     - gmax: the parameter gmax used during the analysis
###     - which histone marks have a the signal that is greater than the median global signal
###     - which histone marks have a signal where the histone are greater than the mean of the global signal. 
###     - a data.frame where in the first column there are the results of feature selection analysis with k-means 
###       and in the second column the results of ICRR (index coverage regulatory regions).
################################################################################################
################################################################################################

featSelectionWithKmeans<-function(dftotann, nk, autoK="no",w=1000, gmax=7,outputplot="feature.selection.pdf"){
 
    featureSelDF<-list()
  
    gc()
 
    if(isTRUE(autoK=="yes")){

    dfK<-data.frame(nrow=1:gmax)

    } else {
    nk<-nk  
    dfK<-data.frame(nrow=1:nk)
  
    }
  
    namecol<-colnames(dftotann)[1:ncol(dftotann)]
      
    listHM<-NULL

    for(i in namecol){

    print(i)
  
    if(autoK=="yes"){
  
    selectMatrix<-as.matrix(dftotann[,i])
    resultsMclust<-Mclust(selectMatrix,G=2:gmax)
    k.best <- dim(resultsMclust$z)[2]
    cat("model-based optimal number of clusters:", k.best, "\n")
    clusteringHM<-kmeans(selectMatrix,k.best)
  
    } else {
   
    selectMatrix<-as.matrix(dftotann[,i])
    clusteringHM<-kmeans(selectMatrix,nk)
  
    }
  
    clusterDF<-data.frame(HM=dftotann[,i],cluster=clusteringHM$cluster)

    clusterDF<-clusterDF[order(clusterDF$cluster),]
    if(autoK=="yes"){

    vectorS<-seq(1:gmax)
  
    } else {vectorS<-seq(1:nk)}
  

    meanGroups<-NULL
    
    for(c in vectorS){
    clusterDFmean<-mean(clusterDF[clusterDF$cluster==c,1])
    meanGroups<-c(meanGroups,clusterDFmean)
  
    }
  
    listHM<-c(listHM,i)
  
    meanGroupsDF<-data.frame(meanGroups)
    
    dfK<-cbind(dfK,meanGroups)
    
    colnames(dfK)[2:ncol(dfK)]<-listHM
    
    }
    
    dfK<-dfK[,c(2:ncol(dfK))]
  
    mom<-apply(dfK,2,mean, na.rm=TRUE)

    d<-dist(mom,method="euclidean")
    fit <- hclust(d, method="average") 
    mom<-mom[fit$order]

    w<-w

    meanCol<-apply(dftotann,2,mean)
    ntotCol<-ncol(dftotann)

    dfBalance<-NULL
    stringList<-NULL

    for(i in 1:ntotCol){

    dfTotsubselect<-dftotann[,i]
    totCov<-nrow(dftotann)*1000
    stringName<-colnames(dftotann)[i]
    print(colnames(dftotann)[i])

    lengthTRUE<-length(which(dfTotsubselect > meanCol[i]))
    lengthTRUEw<-lengthTRUE*w

    lengthFALSE<-length(which(dfTotsubselect < meanCol[i]))
    lengthFALSEw<-lengthFALSE*w


    frTRUEwcov<-lengthTRUEw/totCov
    frFALSEwcov<-lengthFALSEw/totCov

    ratioBM=abs(frFALSEwcov-frTRUEwcov)

    dfBalance<-c(dfBalance,ratioBM)
    stringList<-c(stringList, stringName)
    names(dfBalance)<-stringList

    }

    dfBalance<-data.frame(ID=names(dfBalance),ICRR=dfBalance,stringsAsFactors=FALSE)
    results<-data.frame(ID=names(mom),values=mom)

    results.merge<-merge(results,dfBalance,by="ID",stringsAsFactors=FALSE)


    pdf(outputplot)

    par(mfrow=c(2,1))
    plot(results.merge$values,xaxt="n",pch=20,col="red",type="o",las=2,ylab="Normalize signal",xlab="",main="Results of Feature selection")
    axis(side=1, at=1:length(results.merge$ID), labels=results.merge$ID,las=2,cex.axis=0.60)
    abline(h=median(results.merge$values))

    plot(results.merge$ICRR,xaxt="n",pch=20,col="red",las=2,type="o",ylab="ICRR",xlab="",main="ICRR of the histone marks")
    axis(side=1, at=1:length(results.merge$ID), labels=results.merge$ID,las=2,cex.axis=0.60)
    abline(h=0.5)
    dev.off()

    featureSelDF[[1]]<-dfK
    featureSelDF[[2]]<-mom
    featureSelDF[[3]]<-nk
    featureSelDF[[4]]<-gmax
    featureSelDF[[5]]<-mom[which(mom>=median(mom))]
    featureSelDF[[6]]<-mom[which(mom>=mean(mom))]
    featureSelDF[[7]]<-results.merge

   return(featureSelDF)

}
