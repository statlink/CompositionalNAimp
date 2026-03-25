alfa.knnimp <- function(x, k = 2:10, a = seq(0.1, 1, by = 0.1), econ = TRUE ) {

  p <- dim(x)[2]  ## number of components
  sx <- Rfast::rowsums(x)
  pou <- which( is.na(sx) )  ## where are the missing values
  x1 <- x[-pou, , drop = FALSE]  ## rows with good data, i.e. no missing values
  ina <- 1:dim(x1)[1]
  x2 <- x[pou, , drop = FALSE]  ## rows with missing values (bad data)
  tot <- 1 - Rfast::rowsums(x2, na.rm = TRUE)  ## 1 - sum of the rows of the bad data
  if ( min(x, na.rm = TRUE) == 0 )  a <- a[a >0 ]
  lena <- length(a)
  lenk <- length(k)

  if ( econ ) {

    names <- paste("kappa=", k, sep = "")
    y <- sapply(names, function(x) NULL)
    for ( i in 1:lenk )  y[[ i ]] <- x2
    names <- paste("alfa=", a, sep = "")
    xa <- sapply(names, function(x) NULL)

    for ( i in 1:length(pou) ) {
      id <- which( !is.na(x2[i, ]) )
      w1 <- x1[, id, drop = FALSE] / Rfast::rowsums( x1[, id, drop = FALSE] )
      w2 <- x2[i, id, drop = FALSE] / sum(x2[i, id], drop = FALSE)
      NN <- Rnanoflann::nn(w1, w2, k = nrow(w1), method = "jensen_shannon", sorted = TRUE)$indices
      for ( vim in 1:lena ) {
        for ( j in 1:lenk ) {
          nn <- NN[ 1:k[j] ]
          mu <- Compositional::frechet(x1[nn, , drop = FALSE], a[vim])
          mu <- replace(mu, id, NA)
          est <- mu / sum(mu, na.rm = TRUE) * tot[i]
          y[[ j ]][i, -id] <- est[-id]
        }  ##  end  for ( j in 1:length(k) ) {
        xa[[ vim ]] <- y
      }  ##  end  for ( vim in 1:length(a) ) {
    }  ##  end  for ( i in 1:length(pou) ) {
    res <- list(index = pou, xa = xa)

  } else {

    names <- paste("kappa=", k, sep = "")
    y <- sapply(names, function(x) NULL)
    for ( i in 1:lenk )  y[[ i ]] <- x
    names <- paste("alfa=", a, sep = "")
    xa <- sapply(names, function(x) NULL)

    for ( i in 1:length(pou) ) {
      id <- which( !is.na(x2[i, ]) )
      w1 <- x1[, id, drop = FALSE] / Rfast::rowsums( x1[, id, drop = FALSE] )
      w2 <- x2[i, id, drop = FALSE] / sum(x2[i, id], drop = FALSE)
      NN <- Rnanoflann::nn(w1, w2, k = nrow(w1), method = "jensen_shannon", sorted = TRUE)$indices
      for ( vim in 1:lena ) {
        for ( j in 1:lenk ) {
          nn <- NN[ 1:k[j] ]
          mu <- Compositional::frechet(x1[nn, , drop = FALSE], a[vim])
          mu <- replace(mu, id, NA)
          est <- mu / sum(mu, na.rm = TRUE) * tot[i]
          y[[ j ]][pou[i], -id] <- est[-id]
        }  ##  end  for ( j in 1:length(k) ) {
        xa[[ vim ]] <- y
      }  ##  end  for ( vim in 1:length(a) ) {
    }  ##  end  for ( i in 1:length(pou) ) {
    res <- list(index = pou, xa = xa)

  }  ##  end  if ( econ ) {

  res
}
