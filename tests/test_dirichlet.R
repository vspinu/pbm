context("Dirichlet conjugate priors")

test_that("CatDirichlet",{
    
    N <- 1000
    a <- (1:3)/6
    apr <- c(1, 1, 1)
    D <- dirichlet_rng(N, a)
    Y <- apply(D, 1, function(p) sample(1:3, 1, T, prob = p))

    M <- pbm("CatDirichlet",
             DATA = defBC("dc",
                 mixin = pdCat,
                 st = Y, N = 3, 
                 D = defP("pd(conj)(Dirich)",
                     varsize = 3, cix = rep.int(1, N), cix_dim = 1, 
                     hpD = defP("hc",
                         cix = 1, var = apr))))

    M$D$do.mc_st <- T
    update(M, nr_iter = 1000)
    matplot(M$D.$mc_st, type = "l")

    post <- unname(colMeans(M$D.$mc_st[-(1:100), ]))
    actual <- tabulate(Y)
    actual <- actual/sum(actual)
    rbind(actual, post)
    expect_close(post, actual, .003)
})

## ## JAGS
## library(rjags)

## js <- "model {
##     p ~ ddirch(apr)
##     for(ii in 1:N){
##         Y[ii] ~ dcat(p)
##     }
## }"

## jm <- jags.model(textConnection(js))
## update(jm, 1000)
## system.time(out <- coda.samples(jm, "p", 10000))
## jpost <- colMeans(out[[1]])




test_that("CatDirichlet (multivar)",{
    
    N <- 1000
    apr <- c(1, 1, 1, 1)

    a1 <- (1:3)/6
    D1 <- dirichlet_rng(N, a)
    Y1 <- apply(D1, 1, function(p) sample(1:3, 1, T, prob = p))

    a2 <- (3:1)/6
    D2 <- dirichlet_rng(N, a2)
    Y2 <- apply(D2, 1, function(p) sample(1:3, 1, T, prob = p))

    a3 <- c(1, 1, 1, 1)/4
    D3 <- dirichlet_rng(N, a3)
    Y3 <- apply(D3, 1, function(p) sample(1:4, 1, T, prob = p))

    M <- pbm("CatDirichlet",
             DATA = defBC("dc",
                 mixin = pdCat, varsize = 3, 
                 st = cbind(Y1, Y2, Y3), N = c(3, 3, 4), 
                 D = defP("pd(conj)(Dirich)",
                     cix = 1:3, cix_dim = 2,
                     size = 3, varsize = 4,
                     hpD = defP("hc",
                         cix = 1, var = apr))))

    M$D$do.mc_st <- T
    update(M, nr_iter = 10000)
    ## matplot(M$D.$mc_st, type = "l")

    (post <- colMeans(M$D.[["mc_st"]][-(1:100),, ]))
    actual <- rbind(tabulate(Y1, nbins = 4),
                    tabulate(Y2, nbins = 4),
                    tabulate(Y3, nbins = 4))
    (actual <- actual/rowSums(actual))

    expect_close(unname(post), actual, .003)
})
