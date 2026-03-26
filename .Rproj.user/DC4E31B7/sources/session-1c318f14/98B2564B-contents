alfaknnimp.tune <- function(x, k = 2:10, a = seq(-1, 1, by = 0.1), type = "kl", R = 50, graph = FALSE ) {

  if ( min(x, na.rm = TRUE) == 0 )  a <- a[a > 0]  ## checks for any zeros in the data
  pou <- which( is.na(x), arr.ind = TRUE)
  ind <- unique( pou[, 1] )
  xf <- x[-ind, ]
  xm <- x[ind, ]
  nf <- dim(xf)[1]   ;   nm <- dim(xm)[1]
  names <- paste("poia", 1:nm)
  poia <- sapply(names, function(x) NULL)
  for ( i in 1:nm )  poia[[ i ]] <- which( is.na( x[ind[i], ] ) )

  per <- matrix(0, nrow = length(a), ncol = length(k) )
  yf <- xf

  runtime <- proc.time()
  for ( i in 1:R ) {
    yf <- xf
    index <- sample(1:nf, nm)
    for ( j in 1:nm )  yf[index[j], poia[[ j ]]] <- NA
    mod <- CompositionalNAimp::alfa.knnimp(yf, k = k, a = a, econ = TRUE)

    if ( type == "kl" )  {
      for ( vim in 1:length(a) ) {
        for ( j in 1:length(k) ) {
          per[vim, j] <- per[vim, j] + sum( xf[index, ] * log( xf[index, ] / mod[[ 2 ]][[ vim ]][[ j ]] ), na.rm = TRUE ) / nm
        }
      }
    } else if ( type == "js" ) {
      for ( vim in 1:length(a) ) {
        for ( j in 1:length(k) ) {
          per[vim, j] <- per[vim, j] + sum( Compositional::esova( xf[index, ], mod[[ 2 ]][[ vim ]][[ j ]] ) ) / nm
        }
      }
    }
  } ##  end  for (i in 1:R) {
  runtime <- proc.time() - runtime

  perf <- per / R
  colnames(perf) <- paste("k=", k, sep = "")
  rownames(perf) <- paste("alfa=", a, sep = "")

  if ( graph )  filled.contour( a, k, perf, ylab = "k nearest-neighbours", cex.lab = 1.2, cex.axis = 1.2,
                                xlab = expression(paste(alpha, " values") ) )
  opt <- min(perf)
  poia <- as.vector( which(perf == opt, arr.ind = TRUE)[1, ] )
  list( perf = perf, performance = opt, best_a = a[ poia[1] ], best_k = poia[2] + 1, runtime = runtime )

}
