\name{plot.rsc_cv}
\alias{plot.rsc_cv}

\title{
  Plot method for rsc_cv objects
}

\description{
  Plot the cross-validation estimates of the Frobenius loss.
}


\usage{
   \method{plot}{rsc_cv}(x, \dots)
}




\arguments{
  \item{x}{
    Output from \code{\link{rsc_cv}}, that is an S3 object of class \code{"rsc_cv"}.
  }
  \item{\dots}{
    additional arguments passed to  \code{\link[graphics]{plot.default}}.
  }
}


\value{
   Plot the Frobenius loss estimated via cross-validation (y-axis) vs
   threshold values (x-axis). The dotted blue line represents the average
   expected normalized Frobenius loss, while the vertical segments
   around the average  are  \emph{1-standard-error} error bars
   (see \emph{Details} in \code{\link{rsc_cv}}.

   The vertical dashed red line identifies the minimum of the average
   loss, that is the optimal threshold flagged as \code{"minimum"}.  The
   vertical dashed green line identifies the optimal selection flagged
   as \code{"minimum1se"} in the output of \code{\link{rsc_cv}} (see
   \emph{Details} in \code{\link{rsc_cv}}).
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
\donttest{
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

## plot the cross-validation estimates
plot(a)
   
## pass additional parameters to graphics::plot
plot(a , cex = 2)
}
}
