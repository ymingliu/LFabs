#' Generate a sparse adjacency matrix
#'
#' This function generates a sparse adjacency matrix for predictors.
#' The elements of the matrix indicate whether pairs of vertices are adjacent or not in the graph.
#' And each elements is based on the Pearson's correlation coefficients and proper cut-off value.
#' To save the storage space, this function only retains the information of nonzero ones
#' @param x The n x p design matrix.
#'
#' @return A list
#' \itemize{
#'   \item edge - The nonzero elements of the adjacency matrix.
#'   \item edgerow - The row index of each nonzero elements, which has the same length as \code{edge}.
#'   \item edgecol - The column index of each nonzero elements, which has the same length as \code{edge}.
#'   \item gamma - A \code{p} vector. Each elements is the sum of the column of the adjacency matrix.
#'   \item ledge - The number of the nonzero elements.
#' }
#' @references Jian Huang. Shuangge Ma. Hongzhe Li and Cun-Hui Zhang. (2011)
#' \emph{The sparse Laplacian shrinkage estimator for high-dimensional regression}
#' @examples
#' x = matrix(rnorm(200), 10, 20)
#' A <- adjmat(x)

adjmat<-function(x)
{
  n = nrow(x)
  p = ncol(x)
  r = cor(x)
  if (n<=p) {
    z = 0.5 * log((1 + r[upper.tri(r)])/(1 - r[upper.tri(r)]))
    c0 = mean(sqrt(n - 3)*z) + 2*sd(sqrt(n - 3)*z)
    cutoff = (exp(2*c0/sqrt(n - 3)) - 1)/(exp(2*c0/sqrt(n - 3)) + 1)
  } else {
    cutoff = 0
  }
  A = (r)^5 * (abs(r) > cutoff)
  diag(A) = 0

  gamma = colSums(abs(A))
  adj = c(A[upper.tri(A)])
  id = which(adj!=0)
  adj = adj[id]
  adjRow = c(c(sequence(sequence(p-1)))-1)[id]
  adjCol = c(rep(1:(p-1), 1:(p-1)))[id]

  val = list(edge    = adj,
             edgerow = adjRow,
             edgecol = adjCol,
             gamma   = gamma,
             ledge   = length(adj))

  return(val)
}














