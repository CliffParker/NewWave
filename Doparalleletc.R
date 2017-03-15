#' The foreach and %do%/%dopar% operators provide a looping construct that can be viewed as a hybrid of the 
#' standard for loop and lapply function. It looks similar to the for loop, and it evaluates an expression,
#'  rather than a function (as in lapply), but it's purpose is to return a value (a list, by default), rather 
#'  than to cause side-effects. This faciliates parallelization, but looks more natural to people that prefer 
#'  for loops to lapply.
# 
#' The %:% operator is the nesting operator, used for creating nested foreach loops. Type vignette("nested") 
#' at the R prompt for more details.
# 
#' Parallel computation depends upon a parallel backend that must be registered before performing the 
#' computation. The parallel backends available will be system-specific, but include doParallel, which uses 
#' R's built-in parallel package, doMC, which uses the multicore package, and doSNOW. Each parallel backend has 
#' a specific registration function, such as registerDoParallel or registerDoSNOW.
# 
#' The times function is a simple convenience function that calls foreach. It is useful for evaluating 
#' an R expression multiple times when there are no varying arguments. This can be convenient for resampling,
#'  for example.

#' #' # equivalent to rnorm(3)
#' times(3) %do% rnorm(1)
#' 
#' # equivalent to lapply(1:3, sqrt)
#' foreach(i=1:3) %do%
#'   sqrt(i)
#' 
#' # equivalent to colMeans(m)
#' m <- matrix(rnorm(9), 3, 3)
#' foreach(i=1:ncol(m), .combine=c) %do%
#'   mean(m[,i])
#' 
#' # normalize the rows of a matrix in parallel, with parenthesis used to
#' # force proper operator precedence
#' # Need to register a parallel backend before this example will run
#' # in parallel
#' foreach(i=1:nrow(m), .combine=rbind) %dopar%
#'   (m[i,] / mean(m[i,]))
#' 
#' # simple (and inefficient) parallel matrix multiply
#' library(iterators)
#' a <- matrix(1:16, 4, 4)
#' b <- t(a)
#' foreach(b=iter(b, by='col'), .combine=cbind) %dopar%
#'   (a %*% b)
#' 
#' # split a data frame by row, and put them back together again without
#' # changing anything
#' d <- data.frame(x=1:10, y=rnorm(10))
#' s <- foreach(d=iter(d, by='row'), .combine=rbind) %dopar% d
#' identical(s, d)
#' 
#' # a quick sort function
#' qsort <- function(x) {
#'   n <- length(x)
#'   if (n == 0) {
#'     x
#'   } else {
#'     p <- sample(n, 1)
#'     smaller <- foreach(y=x[-p], .combine=c) %:% when(y <= x[p]) %do% y
#'     larger  <- foreach(y=x[-p], .combine=c) %:% when(y >  x[p]) %do% y
#'     c(qsort(smaller), x[p], qsort(larger))
#'   }
#' }
#' qsort(runif(12))
#' 


library(doParallel)
registerDoParallel(cores=2)
foreach(i=1:5,
        .combine=rbind, 
        .inorder = T,
        .multicombine = T,
        .errorhandling = "remove",
        .packages = c("ggplot2"),
        .verbose = T) %dopar% {
          x=sqrt(i)
          y=sin(i)
          c(x,y)
  }


foreach(..., .combine, .init, .final=NULL, .inorder=TRUE,
        .multicombine=FALSE,
        .maxcombine=if (.multicombine) 100 else 2,
        .errorhandling=c('stop', 'remove', 'pass'),
        .packages=NULL, .export=NULL, .noexport=NULL,
        .verbose=FALSE) %dopar% ex


