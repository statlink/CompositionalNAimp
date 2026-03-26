\name{Cross-validation for the Jensen-Shannon based a--k--NN algorithm}
\alias{alfaknnimp.tune}
\title{
Cross-validation for the Jensen-Shannon based \eqn{\alpha}--k--NN algorithm
}

\description{
Cross-validation for the Jensen-Shannon based \eqn{\alpha}--k--NN algorithm.
}

\usage{
alfaknnimp.tune(x, k = 2:10, a = seq(-1, 1, by = 0.1), type = "kl",
R = 50, graph = FALSE)
}
\arguments{
\item{x}{
A matrix with the compositional data, zeros are allowed.
}
\item{k}{
A vector with the number of nearest neighbours to explore.
}
\item{a}{
A vector with the values of \eqn{\alpha} to consider.
}
\item{type}{
The type of distance to use in order to see the difference between the true and the imputed data.
If there are no zeros, The Kullback-Leibler divergence ("kl", default value) may be used, otherwise the (squared) Jensen-Shannon
divergence ("js") may be used.
}
\item{R}{
The number of repetitions to perform.
}
\item{graph}{
Do you want to produce a heatmap with the performance of each combination of \eqn{\alpha} and k?
}
}

\details{
The algorithm performs a repeated cross-validation for the Jensen-Shanon divergence based
\eqn{\alpha}--k--NN algorithm. It selects the optimal pair of \eqn{\alpha} and k.
}

\value{
A list including:
\item{perf}{
A matrix with the average performance (Kullback-Leibler or Jensen-Shannon divergence) for each value
of \eqn{\alpha} and k.
}
\item{performance}{
The value of the optimal performance.
}
\item{best_a}{
The optimal value of \eqn{\alpha}.
}
\item{best_k}{
The optimal value of k.
}
\item{runtime}{
The runtime of the cross-validation procedure.
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
\code{\link{alfa.knnimp}, \link{cv.knnimp}
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
mod <- alfaknnimp.tune( x, k = 2:4, a = c(0.5, 0.75, 1), R = 20 )
}

