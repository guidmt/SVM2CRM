################################################################################################
################################################################################################
###
### METHOD: createSVMinput
### CLASS: ----
###
### Take the output of findFeatureOverlap and the create a positive set of cis-regulatory elements
### 
### 
### 
################################################################################################
################################################################################################
###
### INPUT
###
###     - input = output of findFetureOverlap
###     - label1 = string to define the positive class (default: enhancers)
###     - label2= string to define the negative class (default: not_enhancers)
###     
### OUTPUT
###
###     - Training_Set for support vector machine 
###     
################################################################################################
################################################################################################

createSVMinput<-function(inputpos,inputfull,label1="enhancers",label2="not_enhancers")

  {

    namecol <- colnames(inputpos)
    namecolSelect <- grep(namecol, pattern = "signal", value = TRUE)

    table.final.overlap.subset <- inputpos[, c("chromosome.feat","start.feat","end.feat", namecolSelect)]

    positive.subset <- data.frame(chromosome=table.final.overlap.subset[,c("chromosome.feat")],label = rep(label1, nrow(table.final.overlap.subset)), 
    table.final.overlap.subset[,c(2:ncol(table.final.overlap.subset))], stringsAsFactors = FALSE)


    indexdelete <- which(inputpos$start.feat %in% inputfull$start)

    colnames(positive.subset)[1:5]<-gsub(colnames(positive.subset)[1:5],pattern=".feat",replacement="")

    negative.subset <- inputfull[-c(indexdelete), c("chromosome","start","end", namecolSelect)]
    negative.subset <- data.frame(chromosome=negative.subset[,"chromosome"],label = rep(label2, nrow(negative.subset)), 
    negative.subset[,c(2:ncol(negative.subset))], stringsAsFactors = FALSE)

    data_enhancer_svm <- rbind(positive.subset, negative.subset)

    return(data_enhancer_svm)

}
