################################################################################################
################################################################################################
###
### METHOD: performancerSVM
### CLASS: ----
###
### Take a contingency table and then estimate: tpr.sensitivity, fpr, acc, spc.specificity, precision,
###  ppv, npv, fdr, fscore. Useful to estimate the performance of a model
###
################################################################################################
################################################################################################
###
### INPUT
###
###     - res: confusion matrix derived using predict() function
###
### OUTPUT
###
###     - A data.frame with tpr.sensitivity, fpr, acc, spc.specificity, precision, ppv, npv, fdr, fscore.
###
################################################################################################
################################################################################################
performanceSVM <- function(res) {
    
    true.positive <- res[1, 1]
    false.negative <- res[2, 1]
    false.positive <- res[1, 2]
    true.negative <- res[2, 2]
    
    # sensitivity
    tpr.sensitivity <- true.positive/(true.positive + false.negative)
    
    # falsepositive.rate
    fpr <- false.positive/(false.positive + true.negative)
    
    # accuracy
    acc <- (true.positive + true.negative)/(true.positive + false.negative + true.negative + 
        false.positive)
    
    # specificity
    spc.specificity <- true.negative/(false.positive + true.negative)
    spc.specificity <- 1 - fpr
    
    # precision
    precision <- true.positive/(true.positive + false.positive)
    
    # positive predictive values
    ppv <- true.positive/(true.positive + false.positive)
    
    # negative predictive values
    npv <- true.negative/(true.negative + false.negative)
    
    # false positive rate
    fdr <- false.positive/(false.positive + true.positive)
    
    # fscore:harmonic mean between precision and recall
    fscore <- 2 * spc.specificity * tpr.sensitivity/(spc.specificity + tpr.sensitivity)
    confusionresults <- data.frame(tpr.sensitivity = tpr.sensitivity, fpr = fpr, 
    acc = acc, fscore = fscore, spc.specificity = spc.specificity, ppv = ppv, 
    npv = npv, fdr = fdr)
    
    return(confusionresults)
    
} 

