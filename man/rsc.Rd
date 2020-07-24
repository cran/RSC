\name{rsc}

\alias{rsc}

\title{Robust and Sparse Correlation Matrix Estimator}

\description{
  Compute the Robust and Sparse Correlation Matrix (RSC) estimator 
  proposed in Serra et al. (2018).
}


\usage{
  rsc(cv, threshold = "minimum") 
}


\arguments{
  \item{cv}{
    An S3 object of class \code{"rsc_cv"} (see \code{\link{rsc_cv}}).
  }
  \item{threshold}{
    Threshold parameter to compute the RSC estimate. This
    is a numeric value taken onto the interval (0,1), or it is
    equal to  \code{"minimum"} or \code{"minimum1se"} for selecting the
    optimal threshold according to the selection performed in 
    \code{\link{rsc_cv}}.
   }
}


\details{
   The setting \code{threshold = "minimum"} or \code{threshold =
   "minimum1se"} applies thresholding according to the criteria
   discussed in the \emph{Details} section in \code{\link{rsc_cv}}.
   When \code{cv} is obtained using \code{\link{rsc_cv}} with
   \code{cv.type = "random"}, the default settings for \code{\link{rsc}}
   implements exactly the RSC estimator proposed in Serra et al.,
   (2018).
  
   Although \code{threshold = "minimum"} is the default choice, in
   high-dimensional situations \code{threshold = "minimum1se"} usually
   provides a more parsimonious representation of the correlation
   structure. Since the underlying RMAD matrix is passed through the
   \code{cv} input, any other hand-tuned threshold to the RMAD matrix
   can be applied without significant additional computational
   costs. The latter can be done setting \code{threshold} to any value
   onto the (0,1) interval.

   The software is optimized to handle high-dimensional data sets,
   therefore, the output RSC matrix is packed into a storage efficient
   sparse format using the \code{"dsCMatrix"} S4 class from the
   \code{\link{Matrix}} package. The latter is specifically designed for
   sparse real symmetric matrices.
}



\value{
  Returns a sparse correlaiton matrix of class \code{"dsCMatrix"}
  (S4 class object) as defined in the \code{\link{Matrix}} package. 
}




\section{References}{
   Serra, A., Coretto, P., Fratello, M., and Tagliaferri, R. (2018).
   Robust and sparsecorrelation matrix estimation for the analysis of
   high-dimensional genomics data. \emph{Bioinformatics}, 34(4),
   625-634.  doi:10.1093/bioinformatics/btx642
}



\seealso{
  \code{\link{rsc_cv}}
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

## obtain the RSC matrix with "minimum" flagged solution 
b <- rsc(cv = a, threshold = "minimum")
b
   
## obtain the RSC matrix with "minimum1se" flagged solution 
d <- rsc(cv = a, threshold = "minimum1se")
d

## since the object 'a' stores the RMAD underlying estimator, we can
## apply thresholding at any level without re-estimating the RMAD 
## matrix
e <- rsc(cv = a, threshold = 0.5)
e
}

