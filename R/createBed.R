################################################################################################
################################################################################################
###
### METHOD: createBed
### CLASS: ------
###
### Create bed file format of output SVM.
### 
###
################################################################################################
################################################################################################
###
### INPUT
###	
###	- test_set: test_set produce using LibLinear or every kind of histone marks 
###     - label1: label of positive classes
###     - pred: p$prediction object obtained from predict function
###	- outputfile: name of bed file 
###
### OUTPUT
###
###     -  bed format 
### 
###     
###
################################################################################################
################################################################################################

createBed<-function(test_set,label1,pred,outputfile){

    fullDF<-test_set
    bedoutput<-data.frame(fullDF[,c(1,3:4)],pred$predictions,rep(0,nrow(fullDF)),rep(".",nrow(fullDF)),fullDF[,c(3:4)],prediction2=pred$predictions)
    bedoutput$prediction2<-as.character(bedoutput$prediction2)
    bedoutput<-bedoutput[bedoutput$prediction2==label1,]
    bedoutput$prediction2<-gsub(x=bedoutput$prediction2,pattern=label1,replacement="255,0,0",fixed=TRUE)
    write.table(bedoutput,file=outputfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

}
