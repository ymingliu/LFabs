# LFabs
A Forward and Backward Stagewise algorithm for High-dimensional smoothed partial rank estimation with sparse laplacian shrinkage.


LFabs uses coordinate descent with a fixed step size which consists of both forward and backward steps. At each step, the first-order Taylor's expansion is used to reflect the main part of the increment. For given tau , a turning parameter before the quadratic term, the laplacian penalty can be incorporated into the loss function.

# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("Amie-Liu/LFabs")

# Usage

- [x] [LFabs-manual](https://github.com/Amie-Liu/LFabs/inst/LFabs-manual.pdf) ------------ Details of the usage of the package.

# Example

    library(LFabs)
    library(Matrix)
    library(mvtnorm)
    n = 400
    p = 500
    d = 15
    g = 5
    sig = c(0.5, 1.5)
    rho = 0.9
    error = "contaminate"
    tran = "log"
    censor.rate = 0.1
    block = "Auto"

    dat = generator(n, p, d, g, sig, rho, error, tran, censor.rate, block)
    x = dat$x
    y = dat$y
    status = dat$status

    sigma = 1/sqrt(n)
    w = abs(1/drop(cor(x, y, method = "pearson")))
    w[which(w=="NaN"|w=="Inf")] = max(w[which(w!="NaN"&w!="Inf")])
    model = "spr"
    fit <- cv.LFabs(x, y, status, sigma, w, model)



# References

High-dimensional Smoothed Partial Rank Estimation with Sparse Laplacian Shrinkage for Nonparametric Transformation Survival Model. Manuscript.

# Development
The R-package is developed by Yiming Liu (liuyimingsufe@163.com), Xiao Zhang and Xingjie Shi.




