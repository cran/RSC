\name{rmad}

\alias{rmad}

\title{RMAD correlation matrix}

\description{
   Compute the RMAD robust correlation matrix proposed in Serra et
   al. (2018) based on the robust correlation coefficient proposed in
   Pasman and Shevlyakov (1987).
}


\usage{
  rmad(x , y = NULL, na.rm = FALSE , even.correction = FALSE, num.threads = "half-max")   
}


\arguments{
  \item{x}{
    A numeric vector, a matrix or a data.frame. If \code{x} is a matrix
    or a data.frame, rows of  \code{x} correspond to sample units
    and columns correspond to variables. If \code{x} is a numerical
    vector, and \code{y} is not \code{NULL}, the RMAD correlation
    coefficient between \code{x} and \code{y} is computed. Categorical
    variables are not allowed.
  }
  \item{y}{
    A numerical vector if not \code{NULL}. If both \code{x} and \code{y}
    are numerical vectors, the RMAD correlation coefficient between
    \code{x} and \code{y} is computed.
  }
  \item{na.rm}{
    A logical value, if \code{TRUE} sample observation
    containing \code{NA} values are excluded (see \emph{Details}).
  }
  \item{even.correction}{
    A logical value, if \code{TRUE} a correction
    for the calculation of the medians is applied to reduce the bias
    when the number of samples even (see \emph{Details}).
  }
  \item{num.threads}{
     An integer value or the string \code{"half-max"} (default), specifying the number of threads for parallel execution (see \emph{Details}).
  }
}


\details{
   The \code{rmad} function computes the correlation matrix based on the
   pairwise robust correlation coefficient of Pasman and Shevlyakov
   (1987). This correlation coefficient is based on repeated median
   calculations for all pairs of variables. This is a computational
   intensive task when the number of variables (that is \code{ncol(x)})
   is large.
   
   The software is optimized for large dimensional data sets, the median
   is approximated as the central observation obtained based on the
   \emph{find} algorithm of Hoare (1961) (also known as \emph{quickselect})
   implemented in C language. For small samples this may be a crude
   approximation, however, it makes the computational cost feasible for
   high-dimensional data sets. With the option \code{even.correction
   = TRUE} a correction is applied to reduce the bias for data sets with
   an even number of samples. Although \code{even.correction = TRUE}
   has a small computational cost for each pair of variables, it is
   suggested to use the default \code{even.correction = FALSE} for large
   dimensional data sets.

   The function can handle a data matrix with missing values (\code{NA}
   records). If \code{na.rm = TRUE} then missing values are handled by
   casewise deletion (and if there are no complete cases, an error is
   returned). In practice, if \code{na.rm = TRUE} all rows of
   \code{x} that contain at least an \code{NA} are removed.

   Since the software is optimized to work with high-dimensional data sets,
   the output RMAD matrix is packed into a storage efficient format
   using the \code{"dspMatrix"} S4 class from the \code{\link{Matrix}}
   package. The latter is  specifically designed for dense real symmetric
   matrices. A sparse correlation matrix can be obtained applying
   thresholding using the \code{\link{rsc_cv}} and \code{\link{rsc}}.

   \code{rmad} function supports parallel execution.
   This is provided via \emph{openmp} (http://www.openmp.org), which must be already available on the system at installation time;
   otherwise, falls back to single-core execution.
   For later installation of openmp, the RSC package needs to be re-installed (re-compiled) to provide multi-threads execution. 
   If \code{num.threads > 0}, function is executed using \code{min(num.threads, max.threads)} threads, where \code{max.threads} is the maximum number of available threads. That is, if positive the specified number of threads (up to the maximum available) are used.
   If \code{num.threads < 0}, function is executed using \code{max(max.threads - num.threads, 1)} threads, i.e. when negative \code{num.threads} indicates the number of threads not to use (at least one thread is used).
   If \code{num.threads == 0}, a single thread is used (equivalent to \code{num.threads = 1}).
   If \code{num.threads == "half-max"}, function is executed using half of the available threads (\code{max(max.threads/2, 1)}). This is the default.
}



\value{
  \item{If \code{x} is a matrix or a data.frame}{
    Returns a correlation matrix of class \code{"dspMatrix"} (S4 class object)
    as defined in the \code{\link{Matrix}} package.
  }
  \item{If \code{x} and \code{y} are numerical vectors}{
    Returns a numerical value, that is the RMAD correlation coefficient
    between \code{x} and \code{y}.    
  }
}



\section{References}{
  Hoare, C. A. (1961). Algorithm 65: find.
  \emph{Communications of the ACM}, 4(7), 321-322.

  Musser, D. R. (1997). Introspective sorting and selection algorithms.
  \emph{Software: Practice and Experience}, 27(8), 983-993.
  
  Pasman,V. and Shevlyakov,G. (1987). Robust methods of estimation of
  correlation coefficient. \emph{Automation Remote Control}, 48, 332-340.

  Serra, A., Coretto, P., Fratello, M., and Tagliaferri, R. (2018).
  Robust and sparsecorrelation matrix estimation for the analysis of
  high-dimensional genomics data. \emph{Bioinformatics}, 34(4), 625-634.
  doi: 10.1093/bioinformatics/btx642 
}



\seealso{
  \code{rsc_cv}, \code{rsc}
}







\examples{
## simulate a random sample from a multivariate Cauchy distribution
set.seed(1)
n   <- 100    # sample size
p   <- 7      # dimension
dat <- matrix(rt(n*p, df = 1), nrow = n, ncol = p)
colnames(dat) <- paste0("Var", 1:p)

   
## compute the rmad correlation coefficient between dat[,1] and dat[,2]
a <- rmad(x = dat[,1], y = dat[,2])


## compute the RMAD correlaiton matrix   
b <- rmad(x = dat)
b
}

