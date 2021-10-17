
## print method for the crsc_cv class
print.rsc_cv <-function(x, ...){
    cat("\n")
    cat("================================================\n")
    cat("   Expected Normalized Squared Frobenius Loss   \n")
    cat("================================================\n")
    print(x$loss)
    cat("================================================\n")
    cat("\n")
}



## ## ## plot method for the crsc_cv class
plot.rsc_cv <- function(x, ...){

    ## add check object
    
    tstar <- which(x$loss$Flag == "minimum")
    hstar <- x$loss$Threshold[tstar]
    
    inf_loss <- x$loss$Average - x$loss$SE
    sup_loss <- x$loss$Average + x$loss$SE
    a <- inf_loss[tstar]
    b <- sup_loss[tstar]
    hstar1se <- max(x$loss$Threshold[which(x$loss$Flag == "*")])

    Ylim     <- range(c(inf_loss, sup_loss))
    plot(x$loss$Threshold, x$loss$Average, t='b', ylim = Ylim ,
         pch=20 , col= "#0052A5", lwd = 2,
         main = "RSC Optimal Threshold Selection",
         xlab = "Threshold",
         ylab = "Average loss", ...)
    arrows(x$loss$Threshold, inf_loss, x$loss$Threshold, sup_loss ,
           length=.05, angle=90, code=3, col="#0052A5")
    abline(v = hstar1se,     col = "#31A853", lty=2, lwd=2)
    abline(v = hstar,        col = "#E0162B", lty=2, lwd=2)
}


