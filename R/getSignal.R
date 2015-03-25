################################################################################################
################################################################################################
###
### METHOD: getSignal
### CLASS: ------
###
### This function require as input a list of bed files of each histone marks, and two reference files: one containg positive sites of enhancers and another one of background regions. 
### This function allow to generate a training set where in the rows there are the genomic coordinate and in the columns the signals of the histone marks.
###  The number of columns depends on the window.size. getSignal take the signals of the histone marks considering the windows size (w) around the start’s of the genomic regions provided in the analysis.
###
################################################################################################
################################################################################################
###
### INPUT
###	
###	- bedfilelist: the list of bed files (e.g. character)
###	- chr: what chromosomes consider in the analysis
###	- reference: file with the reference position
###	- win.size: the size of the window (default=5000bp)
###	- bin.size:          original bin size used 
###     - label1:          class of the reference (e.g. enhancers)
### OUTPUT
###
###     - A data.frame where the number of column depends on the windows size and bin.size and obviously the number of HM. 
###       
################################################################################################
################################################################################################
getSignal<-function(bedfilelist,chr,reference,win.size,bin.size,label1="enhancers"){

    Rle<-NULL

    ref<-read.table(file=reference,sep="\t",header=FALSE,stringsAsFactors=FALSE)

    tabrid<-ref[ref[,1] %in% chr,]


    matrixPeak<-data.frame(1:nrow(tabrid))


    db <- GRanges(tabrid[, 1], ranges = IRanges(tabrid[,2], tabrid[, 3]))

    TSS=start(db)

    TSS.GRanges=GRanges(seqnames=Rle(tabrid[,1]), ranges=IRanges(start=TSS-win.size, end=TSS+win.size))


    for(i in bedfilelist) {

    print(i)

    reads<-import.bed(con=i)
    reads<-as(reads,"RangedData")

    reads.select<-reads[chr]

    lengthVector<-length(reads.select@ranges@unlistData@start)

    IRange.reads=GRanges(seqnames=rep(chr,lengthVector),
    ranges=IRanges(reads.select@ranges@unlistData@start, width=reads.select@ranges@unlistData@width))

    win.size=win.size
    winsize=bin.size

    TSS=start(db)

    TSS.GRanges=GRanges(seqnames=Rle(tabrid[,2]), ranges=IRanges(start=TSS-win.size, end=TSS+win.size))

    allss=sapply(TSS, function(x) seq(x-win.size, x+win.size-1, by=winsize))

    ncol.matrix=(win.size*2/winsize)

    allchrs=as.character(seqnames(db))
    allseq=Rle(allchrs,rep(win.size/winsize*2,length(allchrs)))
    allrange=IRanges(start=allss, end=allss+winsize-1)
    TSS.GRanges=GRanges(seqnames=allseq, ranges=allrange)
    tmp=countOverlaps(TSS.GRanges, IRange.reads)
    results=matrix(tmp, ncol=ncol.matrix, byrow=TRUE)

    nreads<-length(reads@ranges@unlistData@start)
    nreads<-length(reads@ranges@unlistData@start)
   
    results<-(results/nreads)
 
    stringcolnames<-rep(i,ncol.matrix)
    results<-matrix(results,ncol=ncol.matrix,byrow=TRUE)
    colnames(results)<-paste(stringcolnames,"signal",sep=".")
    print(nrow(results))
    matrixPeak<-cbind(matrixPeak,results)

}


    list.start<-TSS
    chromosome.column<-rep(chr,nrow(results))
    end.int<-list.start+win.size*2
    collabel<-rep(label1,nrow(results))

    matrixPeak2<-data.frame(chromosome=chromosome.column,label=collabel,start=list.start,end=end.int,matrixPeak)
 
    return(matrixPeak2)

}
