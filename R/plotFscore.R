################################################################################################
################################################################################################
###
### METHOD: plotFscore
### CLASS: --------
###
### Plot the fscore in relation with sensitivity and specificity
################################################################################################
################################################################################################
###
### INPUT
###
###     - tuning_grid: the output of tuningParametersComb.R                   
### OUTPUT
###
###     - A plot where in the x axis there are the sensitivity, in y axis there are the specificity, z-values: the fscore.
###     
###
################################################################################################
################################################################################################
plotFscore<-function(tuningTAB){

    map<-makecmap(tuningTAB$fscore,colFn=jet)
    plot(x=tuningTAB$tpr.sensitivity,y=tuningTAB$spc.specificity,pch = 16, 
    col = cmap(tuningTAB$fscore, map = map),xlab="sensitivity",ylab="specificity")
    hkey(map, 'fscore')
}
