\name{rsc_cv}

\alias{rsc_cv}

\title{Optimal threshold selection for the RSC estimator}

\description{
   Perform cross-validation to select an adaptive optimal threshold for
   the RSC estimator proposed in Serra et al. (2018).
}


\usage{
  rsc_cv(x, cv.type = "kfold", R = 10, K = 10, threshold = seq(0.05, 0.95, by = 0.025),
         even.correction = FALSE, na.rm = FALSE, ncores = NULL, monitor = TRUE) 
}


\arguments{
  \item{x}{
    A matrix or a data.frame. Rows of \code{x} correspond to sample units
    and columns correspond to variables. Categorical variables are not
    allowed.
  }
  \item{cv.type}{
    A character string indicating the cross-validation algorithm. Possible
    values are \code{"kfold"} for repeated K-fold cross-validation, and
    \code{"random"} for random cross-validation (see \emph{Details}).
  }
  \item{R}{
    An integer corresponding to the number of repeated foldings  when
    \code{cv.type = "kfold"}. When \code{cv.type = "random"} \code{R}
    defines the number of random splits (see \emph{Details}).
  }
  \item{K}{
    An integer corresponding to the number of \emph{folds} in K-fold
    cross-validation. Therefore this argument is not relevant when
    \code{cv.type = "random"}.
  }
  \item{threshold}{
    A sequence of reals taken onto the interval (0,1) defining the
    threshold values at which the loss is estimated. 
  }
  \item{even.correction}{
    A logical value. It sets the parameter  \code{even.correction} in
    each of the underlying RMAD computations (see \emph{Details} in
    \code{\link{rmad}}). 
  }
  \item{na.rm}{
    A logical value, it defines the treatment of missing values in 
    each of the underlying RMAD computations (see \emph{Details}).
  }
  \item{ncores}{
    An integer value defining the number of cores used for parallel
    computing. When \code{ncores=NULL} (default), the number \code{r} of
    available cores is detected, and \code{(r-1)} of them are used
    (see \emph{Details}).
  }
  \item{monitor}{
    A logical value. If \code{TRUE} progress  messages are
    printed on screen.
  }
}


\details{
  The \code{rsc_cv} function performs cross-validation to estimate the
  expected Frobenius loss proposed in Bickel and Levina (2008). The
  original contribution of Bickel and Levina (2008), and its extension
  in Serra et al. (2018), is based on a random
  cross-validation  algorithm where the training/test size depends on
  the sample size \emph{n}. The latter is implemented selecting
  \code{cv.type = "ramdom"}, and fixing an appropriate number \code{R} of random
  train/test splits. \code{R} should be as large as possible, but
  in practice this impacts the computing time strongly for
  high-dimensional data sets.  

  Although Serra et al. (2018) showed that the random cross-validation
  of Bickel and Levina (2008) works well for the RSC estimator,
  subsequent experiments suggested that repeated K-fold cross-validation
  on average produces better results. Repeated K-fold cross-validation
  is implemented with the default \code{cv.type = "kfold"}. In this case
  \code{K} defines the number of \emph{folds}, while \code{R} defines
  the number of times that the K-fold cross-validation is repeated with
  \code{R} independent shuffles of the original data. Selecting
  \code{R=1} and \code{K=10} one performs the standard 10-fold
  cross-validation. Ten replicates (\code{R=10}) of the K-fold
  cross-validation are generally sufficient to obtain reasonable
  estimates of the underlying loss, but for extremely high-dimensional
  data \code{R} may be varied to speed up calculations. 

  On multi-core hardware the cross-validation is executed in parallel
  setting \code{ncores}. The parallelism is implemented on the
  total number of data splits, that is \code{R} for the random
  cross-validation, and \code{R*K} for the repeated K-fold
  cross-validation. The software is optimized so that generally the
  total computing time scales almost linearly with the number of
  available computer cores (\code{ncores}). 

  For both the random and the K-fold cross-validation it is computed the
  normalized version of the expected squared Frobenius loss proposed in
  Bickel and Levina (2008). The normalization is such
  that the squared Frobenius norm of the identity matrix equals to 1
  whatever is its dimension.

  Two optimal threshold selection types are reported with flags (see
  \emph{Value} section below): \code{"minimum"} and
  \code{"minimum1se"}. The flag \code{"minimum"} denotes the threshold
  value that minimizes the average loss. The flag \code{"minimum1se"}
  implements the so called
  \emph{1-SE rule}: this is the maximum threshold value such that the
  corresponding average loss is within \emph{1-standard-error} with
  respect to the threshold that minimizes the average loss
  (that is the one corresponding to the \code{"minimum"} flag). 

  Since unbiased standard errors for the K-fold cross-validation are
  impossible to compute (see Bengio and Grandvalet, 2004), when
  \code{cv.type="kfold"} the reported standard errors have to be
  considered as a downward biased approximation.
}



\value{
  An S3 object of class \code{'cv_rsc'} with the following components:
  \item{rmadvec}{
    A vector containing the lower triangle of the underlying RMAD
    matrix.
  }
  \item{varnames}{
    A character vector if variable names are available for the input
    data set \code{x}. Otherwise this is \code{NULL}.
  }
  \item{loss}{
    A data.frame reporting cross-validation estimates. Columns of
    \code{loss} are as follows:   \code{loss$Threshold} is the threshold value;
    \code{loss$Average} is averaged loss;   \code{loss$SE} is the standard error
    for the average loss; \code{loss$Flag="minimum"} denotes the  threshold
    achieving the minimum average loss; \code{loss$Flag = "*"} denotes threshold
    values such that the average loss is within \emph{1-standard-error}
    with respect to the \code{"minimum"} solution.
  }
  \item{minimum}{
    A numeric value. This is the minimum of the average loss. This
    corresponds to the flag \code{"minimum"} in the  loss component
    above (see \emph{Details}).
  }
  \item{minimum1se}{
     A numeric value. This is the largest threshold such that the
     corresponding \code{flag = "*"}. In practice this selects the
     optimal threshold based on the \emph{1-SE rule} discussed in the
     \emph{Details} Section above.
  }
}




\section{References}{
   Bengio, Y., and Grandvalet, Y. (2004). No unbiased estimator of the
   variance of k-fold cross-validation. \emph{Journal of Machine Learning
   Research}, 5(Sep), 1089-1105.

   Bickel, P. J., and Levina, E. (2008). Covariance regularization by
   thresholding. The \emph{Annals of Statistics}, 36(6), 2577-2604.
   doi:10.1214/08-AOS600
   
   Serra, A., Coretto, P., Fratello, M., and Tagliaferri, R. (2018).
   Robust and sparsecorrelation matrix estimation for the analysis of
   high-dimensional genomics data. \emph{Bioinformatics}, 34(4),
   625-634.  doi:10.1093/bioinformatics/btx642
}



\seealso{
   \code{rsc}, \code{plot.rsc_cv}
}







\examples{
## simulate a random sample from a multivariate Cauchy distribution
## note: example in high-dimension are obtained increasing p
set.seed(1)
n   <- 100  # sample size
p   <- 10   # dimension
dat <- matrix(rt(n*p, df = 1), nrow = n, ncol = p)
colnames(dat) <- paste0("Var", 1:p)

   
## perform 10-fold cross-validation repeated R=10 times
## note: for multi-core machines experiment with 'ncores'
set.seed(2)
a <- rsc_cv(x = dat, R = 10, K = 10, ncores = 1)
a

   
## threshold selection: note that here, knowing the sampling designs,
## we would like to threshold any correlation larger than zero in
## absolute value
## 
a$minimum        ## "minimum"    flagged solution 
a$minimum1se     ## "minimum1se" flagged solution

## plot the cross-validation estimates
plot(a)

## to obtain the RSC matrix we pass 'a' to the rsc() function
b <- rsc(cv = a, threshold = "minimum")
b

d <- rsc(cv = a, threshold = "minimum1se")
d

## since the object 'a' stores the RMAD underlying estimator, we can
## apply thresholding at any level without re-estimating the RMAD 
## matrix
e <- rsc(cv = a, threshold = 0.5)
e
}

