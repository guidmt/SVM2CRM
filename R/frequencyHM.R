################################################################################################
################################################################################################
###
### METHOD: frequencyHM
### CLASS: --------
###
### Compute the frequency of the histone marks after tuning of parameters (tuningParametersComb)
################################################################################################
################################################################################################
###
### INPUT
###
###     - tuningTABfilter: a data.frame from tuningParametersComb
### OUTPUT
###
###     - A vector which the frequencies of the histone marks
###     
###
################################################################################################
################################################################################################
frequencyHM<-function(tuningTABfilter){

    frequencyHM<- NULL
    for (i in unique(tuningTABfilter[,1])){
    print(i)
    select<-tuningTABfilter[which(tuningTABfilter[,1]==i ),]
    string3<-as.character(gsub(select[,4],pattern="//",replacement=","))
    STRING<-strsplit(string3,split=",")
    string<-unlist(STRING)
    frequencyHM<-c(frequencyHM,string)
    }

    return(frequencyHM)

}
