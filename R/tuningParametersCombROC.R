################################################################################################
################################################################################################
###
### METHOD: tuningParametersCombROC
### CLASS: ------
### 
### Test different models using changing the kernel of SVM, the values of cost functions, the
### number of histone marks. 
### 
################################################################################################
################################################################################################
###
### INPUT
###     
###     - training_set: a data.frame created with buildTrainingSet
###     - typeSVM: a vector with the kinds of Kernel
###     - costV: a vector with a list of cost parameters
###     - different.weight: the data are unbalanced (default TRUE)
###     - vecS:a vector with the number of histone marks (e.g. from 2 to x)
###     - infofile: a data.frame where in the column "a" there are the kinds of histone marks, while in the column "b" a vectors of letters.
###     - pcClass.string: label of the first class ("enhancer")
###     - ncClass.string: label of the second class ("not_enhancers")
###     - pcClass: number of positive class in the test set
###     - ncClass: number of negative class in the test set
###
### OUTPUT
###
###     - a data.frame with the values from perfomanceSVM for each trained model.
###       
###
################################################################################################
################################################################################################

tuningParametersCombROC <- function(training_set, typeSVM, costV,different.weight=TRUE, vecS,infofile,pcClass.string="enhancers",nClass.string="not_enhancers",pcClass,ncClass) {
        
    resultsDFtotHM <- data.frame()
    
    for (v in vecS) {
        
    if(v==max(vecS)){

    sub.info<-infofile[,1]

    combinationHM <- combn(as.character(sub.info), v)

    tcombinationHM <- t(combinationHM)

    print("last comparison")

    print(dim(tcombinationHM))

     } else {
                   
     sub.info<-infofile[,1]
        
     combinationHM <- combn(as.character(sub.info), v)
     print(dim(combinationHM))
     tcombinationHM <- t(combinationHM)

     total<-nrow(tcombinationHM) 
     pb   <- txtProgressBar(1, total, style=3)
     }

     print(dim(tcombinationHM))
       
     print("k-fold-validation in all combination")                                    


    for (n in 1:nrow(tcombinationHM)) {

        if(v<max(vecS)){

        Sys.sleep(0.0000001)
        setTxtProgressBar(pb, n)

    }

    selectHM <- tcombinationHM[n, ]
    print(selectHM)     
    resultsDFtot <- data.frame()
            
        for (t in typeSVM) {
                               
        cost <- costV
                
        resultsDF <- data.frame()
                
        for (c in cost) {

        matchlist<-gsub(x=gsub(paste(x=selectHM,collapse="@"),pattern="-",replacement="."),pattern=".signal",replacement="")

        Totcolnames<-colnames(training_set[, c(5:ncol(training_set))]) 
               
        z<-strsplit(matchlist,split="@")[[1]]

        selectFULL<-NULL

        for(i in 1:length(z)){

        selectHM2<-grep(x=Totcolnames,pattern=z[i],value=TRUE)
        selectFULL<-c(selectFULL,selectHM2)

        }

        training_set$label<-as.character(training_set$label)

        pc <- sample(which(training_set$label == pcClass.string), pcClass)
        nc <- sample(which(training_set$label == nClass.string), ncClass)
        rowtset <- c(pc, nc)

        xTraining_set <- training_set[rowtset,selectFULL]

        yTraining_set <- factor(training_set[rowtset, "label"])

        aucv<-plotROC(datatrain=xTraining_set,y=yTraining_set,different.weight,k=3,type=t,cost=c,output="roctest_training_set_comparison_FS_test")
        aucvAUC<-aucv[[1]]
        medianAUC<-median(unlist(aucvAUC@y.values))
        PR<-aucv[[2]]
        statisticsPSVM<-apply(PR,2,median)
        CORdf<-aucv[[3]]
        statisticsCOR<-median(CORdf)               
        Pvalcolumn<-aucv[[4]]
        pvalestimate<-median(Pvalcolumn)

        stringHM <- gsub(pattern = "signal.output_", paste(selectHM,collapse = "//"), replacement = "")
                    
        resultsDFrow <- data.frame(nHM=v,type = t, cost = c, HM = paste(stringHM, sep = "//", collapse = ""), AUC_k_mean=medianAUC,t(statisticsPSVM),cor.training=statisticsCOR,p.val=pvalestimate)

        resultsDF <- rbind(resultsDF, resultsDFrow)
                                        
            
        }
                
        resultsDFtot <- rbind(resultsDFtot, resultsDF)
                
        }
        resultsDFtotHM <- rbind(resultsDFtotHM, resultsDFtot)
            
        }

    }     
    return(resultsDFtotHM)
}
