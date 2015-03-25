################################################################################################
################################################################################################
###
### METHOD: cisREfindbed
### CLASS: ------
###
### This function require a sorted .bed file in format chr, start, end, signalNorm. The data  
### should be normalize external using others method. For e.g the ChIP-seq data should be
### partitioned in 100bp and then normalize by the library size. This function simply
### load each bed file of one histone mark and then, partion the genome in n not overlapping windows of size (w) . In
### particular, build a matrix where the number of columns for each histone mark depends on the size of the window and
### the bin size used during the preprocessing of the bed file. This function allow to smooth the signal
### of the bed file every n bin. The default function of smoothing is median. This is the suggested
### function to model the signal. Other function have not been tested.
################################################################################################
################################################################################################
###
### INPUT
###	
###	- list.file: the list of bed files (e.g. character)
###	- chr: what chromosomes consider in the analysis
###	- bin.size: the size of bin used to preprocessed the bed files
###	- windows: the size of the window (default=5000bp)
###     - windows.smooth: the size of the window to smooth the bin (default=200)
###   	- smoothing: logical, default is FALSE, if TRUE smooth the signal
###	- function.smoothing: set the function to smooth the signal (default is median)
### OUTPUT
###
###     - A data.frame where the number of column depends on the windows size and bin.size and obviously the number of HM. 
###       
################################################################################################
################################################################################################

cisREfindbed<-function(list_file,chr="chr1", bin.size, windows, window.smooth=200,smoothing="FALSE",function.smoothing="median"){

    completeTABLE<-data.frame(data.frame(ncol = 2 * length(list_file)))

    windows<-windows+windows
    bin.col.table<-windows/bin.size


    for (i in list_file) {

    print(chr)
    print(i)
    print(windows)
    print(bin.size)

    print(getwd())

    reads <- import.bed(con=i,extraCols = "valuesN")


    reads<-as(reads,"RangedData")
    reads.select<-reads[chr]

    lengthVector<-length(reads.select@ranges@unlistData@start)

    IRange.reads=GRanges(seqnames=rep(chr,lengthVector),
    ranges=IRanges(sort(reads.select@ranges@unlistData@start), width=reads.select@ranges@unlistData@width),value=reads.select@values@unlistData@listData)


    if(isTRUE(smoothing=="FALSE")){
    print("processing of signals: not smoothing")

    s<-as.numeric(reads.select@values@unlistData@listData$NA.)
    matrixsignal=matrix(s, ncol=bin.col.table, byrow=TRUE)
    print(dim(matrixsignal))
    matrixsignal<-data.frame(matrixsignal)
    stringcolnames<-rep(i,bin.col.table)
    colnames(matrixsignal)<-paste(stringcolnames,"signal",sep=".")
    chromosome.column<-rep(chr,nrow(matrixsignal))

    list.start<-start(IRange.reads)

    lengthV<-length(start(IRange.reads))

    index.start<-seq(1,lengthV,by=bin.col.table)
    start.int<-list.start[index.start]
    end.int<-start.int+windows



    resultsfull<-data.frame(chromosome=chromosome.column,index=index.start,start=start.int,end=end.int,matrixsignal)

    } else {


    s<-as.numeric(reads.select@values@unlistData@listData$NA.)

    step<-window.smooth/bin.size
    bin.col.table<- windows/window.smooth



    print("processing of signals: smoothing")
    smoothsignal<-rollapply(s,width=step,FUN=function.smoothing,na.rm=FALSE)
    matrixsignal=matrix(smoothsignal, ncol=bin.col.table, byrow=TRUE)
    matrixsignal<-data.frame(matrixsignal)
    stringcolnames<-rep(i,bin.col.table)
    colnames(matrixsignal)<-paste(stringcolnames,"signal",sep=".")
    chromosome.column<-rep(chr,nrow(matrixsignal))

    list.start<-start(IRange.reads)
    
    lengthV<-length(start(IRange.reads))
    
    index.start<-seq(1,lengthV,by=bin.col.table)
    start.int<-list.start[index.start]
    end.int<-start.int+windows

    resultsfull<-data.frame(chromosome=chromosome.column,index=index.start,start=start.int,end=end.int,matrixsignal)
    }

    completeTABLE<-data.frame(completeTABLE,resultsfull)

    matchlist<-c("chromosome","start","end","index")
    
    colremove<-unique(grep(paste(matchlist,collapse="|"), colnames(completeTABLE)[5:ncol(completeTABLE)], value=TRUE,invert=TRUE))

    completeTABLE<-data.frame(completeTABLE[,c(1:5)],completeTABLE[,colremove],stringsAsFactors=FALSE)

}

return(completeTABLE)

}
