\name{cisREfindbed}
\alias{cisREfindbed}
\alias{cisREfindbed}
\title{
     Create and preprocessing of the input for SVM2CRM analysis
}
\description{
This function require as input a bed file in format chr, start, end, signalNorm. The genomic coordinates must be sort.
The data  should be normalized using others methods. For e.g the ChIP-seq data should be partitioned in 100bp and then normalize by the library size. 
This function simply load each bed file of one histone mark and then, partion the genome in n not overlapping windows of size (w). In
particular, build a matrix where the number of columns for each histone mark depends on the size of the window and
the bin size used during the preprocessing of the bed file. This function allow to smooth the signal
of the bed file every n bin. The default function of smoothing is median. This is the suggested
function to model the signal. Other function have not been tested.
}
\usage{
 cisREfindbed(list_file,chr="chr1", bin.size, windows, window.smooth=200,smoothing="FALSE",function.smoothing="median")
}
\arguments{
    \item{list_file}{
         The list of bed files (e.g. character)
}
    \item{chr}{
         A vector containing the list of chromosomes to consider during the analysis
}
   \item{bin.size}{
         The size of bin used to preprocessed the bed files
}
   \item{windows}{
         The size of the window (default=5000bp)
}
   \item{window.smooth}{
         The size of the window to smooth the bin (default=200)

}
  \item{smoothing}{
    logical, default is FALSE, if TRUE smooth the signal
}
 \item{function.smoothing}{
    Set the function to smooth the signal (default is median)
}

}
\details{
    Some detailed description
}
\value{
    A data.frame where the number of columns depends on the windows size, bin.size and the number of histone marks considered for the prediction.

}
\author{
    Guidantonio Malagoli Tagliazucchi guidantonio.malagolitagliazucchi@unimore.it
}
\seealso{cisREfindbed, mclapply}

\examples{
   library(rtracklayer)
   chr<-"chr1"
   #the values use to binarized the genome (in bp)
   bin.size<-100
   #the size of the windows (in bp)
   windows<-500
   #the step of smoothing (in bp): e.g. smooth the signal every 200bp
   window.smooth<-200
   #do you want apply a smoothing
   smoothing<-"FALSE"
   #what kind of function do you want use to smooth the signal
   function.smoothing<-"median"
   #list of file: format bed
   list_file<-"CD4-H3K14ac.norm.w100.bed"

#completeTABLE<-cisREfindbed(list_file[1],chr=chr,bin.size=bin.size,windows=windows,window.smooth=window.smooth,smoothing="FALSE",function.smoothing="median")

}
\keyword{preprocessing}

