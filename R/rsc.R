rsc <- function(cv, threshold = "minimum"){

    ## inputs
    ## cv         = u   ## a class cv_rsc or any other correlation matrix
    ## threshold  = "minimum" ## "minimum", "minimum1se" or numeric in (0,1)

    if(class(cv) == "rsc_cv"){

        ## check threshold
        if(is.numeric(threshold)){
            if(length(threshold)>1){
                stop("if a specific value for 'threshold' is chosen, this must be a single numeric value in (0,1)")
            }else if(threshold <=0 | threshold >=1){
                stop("if a specific value for 'threshold' is chosen, this must be a single numeric value in (0,1)")
            }
        }else{
            if({threshold != "minimum"} & {threshold != "minimum1se"}){
                stop("'threshold' must be one of the following: 'minimum', 'minimum1se', a numeric value in (0,1).")
            }

            if(threshold == "minimum"){
                threshold <- cv$minimum
            }else if(threshold == "minimum1se"){
                threshold <- cv$minimum1se
            }
        }


        ## threshold the rmadvec 
        cv$rmadvec[ abs(cv$rmadvec) < threshold ] <- 0

        nc <- length(cv$rmadvec)
        p  <- {1 + sqrt( 1 + 8 * nc ) } / 2
        R  <- Matrix(1, nrow = p, ncol = p, sparse = TRUE)

        R[lower.tri(R , diag = FALSE)]  <- cv$rmadvec
        R <- forceSymmetric(R , uplo="L")

        ## attach dimnames if needed 
        if(!is.null(cv$varnames)){
            dimnames(R)[[1]] <- dimnames(R)[[2]] <- cv$varnames
        }
    }else{
        stop("'cv' must be a an object of class 'rsc_cv' obtained from 'rsc::rsc_cv'")
        }


    return(R)
}
