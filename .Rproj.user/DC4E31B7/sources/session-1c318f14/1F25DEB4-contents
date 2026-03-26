\name{Imputation of missing values using the Jensen-Shannon based k--NN algorithm}
\alias{alfa.knnimp}
\title{
Imputation of missing values using the Jensen-Shannon based k--NN algorithm
}

\description{
Imputation of missing values using the Jensen-Shannon based k--NN algorithm.
}

\usage{
alfa.knnimp(x, k = 2:10, a = seq(0.1, 1, by = 0.1), econ = TRUE)
}

\arguments{
\item{x}{
A matrix with the compositional data, zeros are allowed.
}
\item{k}{
The number of nearest neighbours to use. This can be a single number or a numerical vector.
}
\item{a}{
A vector (or a single number) with the values of \eqn{\alpha} to consider.
}
\item{econ}{
This argument is for memory saving purposes. If this is TRUE, the function will return
only the data with the imputed missing values, otherwise it will return the initial compositional data with the imputed values.
Evidently, the second option requires more memory.
}
}

\details{
Missing value imputation is performed using the k--NN algorithm with the (squared) Jensen-Shanon divergence.
}

\value{
A list including:
\item{index}{
A vector with the rows (indexes) that contain missing values.
}
\item{xa}{
A list with as many elements as the number of the values of the argument "a" (\eqn{\alpha} values),
and inside of each element there is a sub-list with as many elements as the number of neighbours (length of k basically).
If the argument "econ" is TRUE, then each element of the list containts a matrix of the imputed compositions only.
If the argument "econ" is FALSE, then each element contains the initial compositional data and the imputed missing
values.
}
}

\references{
Tsagris M., Alenazi A. and Stewart C. (2026). A Jensen--Shannon divergence based k--NN algorithm
for missing value imputation in compositional data.
Journal of Applied Statistics (Accepted for publication).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{cv.knnimp}
}
}

\examples{
x <- matrix( rgamma(100 * 7, runif(7, 4,6), 1), ncol = 7, byrow = TRUE )
x <- w <- x / Rfast::rowsums(x)  ## Dirichlet simulated values
dm <- dim(x)    ;    n <- dm[1]    ;    p <- dm[2]
xam <- matrix(nrow = 10, ncol =  3)
nu <- sample(1:n, 10)
for ( vim in 1:10)  xam[vim, ] <- sample(1:p, 3)
for ( vim in 1:10)  x[ nu[vim], xam[vim, ] ] <- NA
mod <- knnimp( x, k = 2:4, a = c(0.1, 0.2) )
}

