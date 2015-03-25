################################################################################################
################################################################################################
###
### METHOD: predictionGW
### CLASS: ------
### Perform prediction of cis-regulatory elements genome-wide
################################################################################################
################################################################################################
###
### INPUT
###     
###     - training_set: The data.frame with the training set
###     - data_enhancer_svm:  The signals of all histone marks along the genome
###     - listHM: a vector with the histone marks that you want to use perform the prediction
###     - pcClass.string: label of the first class ("enhancer")
###     - ncClass.string: label of the second class ("not_enhancers")
###     - pcClass: number of positive class in the test set
###     - ncClass: number of negative class in the test set
###     - cost: parameter of svm (default=100)
###     - type: type of kernel (default=0)
###     - output.file: name of the bed file of output
### OUTPUT
### - The performance of prediction
### - A bed file with the list of enhancers and not_enhancers regions
################################################################################################
predictionGW<-function(training_set,data_enhancer_svm,listHM,pcClass.string="enhancer",nClass.string="not_enhancers",pcClass,ncClass,cost=100,type=0,output.file){

    test_set<-data_enhancer_svm

    pc <- sample(which(training_set$label == pcClass.string), pcClass)
    nc <- sample(which(training_set$label == nClass.string), ncClass)
    rowtset <- c(pc, nc)

    xTraining_set <- training_set[, listHM]

    yTraininig_set <- factor(training_set[, "label"])

    wl<-(table(yTraininig_set))
    w2<-wl[2]/wl[1]
    wl<-round(c(1,w2))
    names(wl)<-names((table(yTraininig_set)))

    aucv<-plotROC(datatrain=xTraining_set,y=yTraininig_set,different.weight=TRUE,k=3,type=type,cost=cost,output=output.file)

    test_set[is.na(test_set)] <- 0

    xTest_set <- test_set[, listHM]

    yTest_set <- factor(test_set[, "label"])

    s <- scale(xTraining_set, center = TRUE, scale = TRUE)

    m = LiblineaR(data = s, target = yTraininig_set, type =type, cost =cost,wi=wl)

    s2 = scale(xTest_set, attr(s, "scaled:center"), attr(s, "scaled:scale"))

    p <- predict(m, s2, decisionValues =FALSE, proba =TRUE)

    res <- table(prediction.svm = p$predictions, true.values = yTest_set)


    fileoutput<-paste(output.file,".bed",sep="")
    createBed(test_set,label1="enhancer",p,fileoutput)

    return(aucv)

}


