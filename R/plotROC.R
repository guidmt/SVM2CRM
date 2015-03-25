################################################################################################
################################################################################################
###
### METHOD: plotROC
### CLASS: --------
###
### Plot the ROC curve of the best model
################################################################################################
################################################################################################
###
### INPUT
###
###     - data: xTraining_set
###     - y: yTraining_set
###     - k: number of iterations for cross-correlations plots
###     - output: namefile output (.pdf)        
###     - different.weight=the data are unbalanced         
### OUTPUT
###
###     - A pdf with the k-cross-correlation plots for: 1) k-cross-validation, 2) ROC horizontal 3) ROC vertical
###     
###
################################################################################################
################################################################################################

plotROC<-function(datatrain,k,y,different.weight,type,cost,output){
  
    index<-cvsegments(nrow(datatrain),k) 

    ROC.probabilities<-data.frame()
    ROC.labels<-data.frame()
    PR<-data.frame()
    CORdf<-NULL
    Pvalcolumn<-NULL

    for(i in 1:length(index)) 

    { 

    x<-index[i][[1]]
    testset<-datatrain[x,] 
    trainset<-datatrain[-x,] 

    testlab<-y[x]
    trainlab<-y[-x]

    indexLab<-table(trainlab)

    indexpos<-indexLab[1]
    indexneg<-indexLab[2]


    if (indexpos<indexneg){
    indexselect<-sample(indexpos:indexneg,indexpos)
    corK<-cor(c(as.matrix(trainset[1:indexpos,])),c(as.matrix(trainset[indexselect,])))
    } else {

    indexselect<-sample(1:indexpos,indexneg)
    corK<-cor(c(as.matrix(trainset[indexselect,])),c(as.matrix(trainset[indexpos:(nrow(trainset)-1),])))

    }

    if(different.weight=="TRUE"){

    w<-(table(trainlab))
    w2<-w[2]/w[1]
    w<-round(c(1,w2))
    names(w)<-names((table(trainlab)))

    } else {

    wlabel<-names(table(trainlab))
    w<-c(1,1)
    names(w)<-wlabel}

    s<-scale(trainset,center=TRUE,scale=TRUE)

    m=LiblineaR(data=s,target=trainlab,type=type,cost=cost,verbose=FALSE, wi = w,bias=FALSE)
     
    s2<-scale(testset,attr(s,"scaled:center"),attr(s,"scaled:scale"))

    pred = predict(m, s2, proba=TRUE)

    predictionColumn<-pred$probabilities[,"enhancers"]

    lab<-testlab
    lab<-gsub(lab,pattern="not_enhancers",replacement=0)
    lab<-gsub(lab,pattern="enhancers",replacement=1)
    lab<-as.numeric(lab)

    ROC.probabilities<-rbind(ROC.probabilities,predictionColumn)
    ROC.labels<-rbind(ROC.labels,lab)

    pval<-roc.area(lab,predictionColumn)

    res <- table(prediction.svm = pred$predictions, true.values = testlab)

    if (nrow(res) > 1) {

    results <- performanceSVM(res)

    PR<-rbind(PR,res=results)
    CORdf<-c(CORdf,corK)
    Pvalcolumn<-c(Pvalcolumn,pval$p.value)
    }


    }
 
    ROC.probabilities <- t(data.frame(ROC.probabilities))
    ROC.labels<- t(data.frame(ROC.labels))
 
    ROC.object<-list()
    ROC.object$predictions<-ROC.probabilities
    ROC.object$labels<-ROC.labels
 
    pred2 <- prediction(ROC.object$predictions, ROC.object$labels)
    perf <- performance(pred2,"tpr","fpr")

    namefile<-paste(output,sep=".","pdf")

    pdf(namefile)
    aucv<- performance(pred2,"auc")


    plot(perf,col="red",lty=3)

    plot(perf,col="red",avg="vertical")


    plot(perf,col="red",avg="horizontal")
    dev.off()

    return(list(aucv,PR,CORdf,Pvalcolumn))

}
