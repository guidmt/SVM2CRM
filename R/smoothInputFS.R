################################################################################################
################################################################################################
###
### METHOD: smoothInputFS
### CLASS: ------
### Give the matrix obtained using getSignal this functions smooth the signals of each histone marks using a particular window (if bin=100).To size of smooth is  bin*k (e.g. a parameter k equal to 2 means thatthe signal is smooth every 200bp).
################################################################################################
################################################################################################
###
### INPUT
###
###     - input_ann:          The data.frame with the training set
###     - k:  The size to set to smooth
###     - listcolnames: The names of column in which perform the smoothing. A vector with the list of histone marks.
### OUTPUT
###
###     - Give the matrix obtained using getSignal this functions smooth the signals of each histone marks using a particular window (if bin=100).To size of smooth is  bin*k (e.g. a parameter k equal to 2 means thatthe signal is smooth every 200bp).
###       
################################################################################################
################################################################################################
smoothInputFS<-function(input_ann,k,listcolnames){

    listHM<-unique(listcolnames)

    dftotann<-data.frame(rownames=1:nrow(input_ann))

    for (i in listHM) {

    print(i)

    selectHM<-grep(x=colnames(input_ann),pattern=i,value=TRUE)
  
    x<-1:length(selectHM)
  
    re.chunk<-split(x,as.numeric(gl(length(x),k,length(x)))) 
    
    lengthREchunk<-length(re.chunk)
  
    dfPartial<-data.frame(row.names=1:nrow(input_ann))
  
    for(n in 1:lengthREchunk){
  
    input_annSHM<-input_ann[,selectHM]
  
    selectIndexcol<-re.chunk[[n]]
  
    dfHMmedian<-apply(input_annSHM[,selectIndexcol],1,mean)
  
    dfPartial<-cbind(dfPartial,dfHMmedian)      
  
    }
  
    print("rename")
    newstring<-paste(rep(i,n),seq(1,ncol(dfPartial)),sep="_")
    colnames(dfPartial)<-newstring
  
    dftotann<-cbind(dftotann,dfPartial,stringsAsFactors=FALSE)
 
    }

    dftotann<-dftotann[,c(2:ncol(dftotann))]
    dftotann<-na.omit(dftotann)

    return(dftotann)

}
