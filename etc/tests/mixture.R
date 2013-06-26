## adapted from http://doingbayesiandataanalysis.blogspot.nl/2012/06/mixture-of-normal-distributions.html
trueM <- c(100, 145)
trueSD <- 15
N <- 1000
K <- 2
trueCL <- sample(1:K, N, T, c(.2, .8))
Y <- rnorm(N, trueM[trueCL] , trueSD)
burnin <- 1:500
hpDir = rep.int(.0001, K)

M <- pbm("NMix",
         D = defBC("dc.", distr = pdNormal, 
             st = Y, 
             mu = defP("dnorm",
                 cix = quote(pst("clust")),
                 scale = 1, size = 2, 
                 var = mean(Y), 
                 hp_mu = defP("hc",
                     var = c(mean = 0, tau = .00005))),
             tau = defP("dlnorm",
                 cix = 1, scale = .2, 
                 hp_tau = defP("hc",
                     var = c(meanlog = 1, taulog = .00005))),
             clust = defP("dcat",
                 cix = 1:N,
                 N = K, size = N, 
                 pr_clust = defP("ddirich",
                     cix = 1, cix_dim = 1, 
                     varsize = K, 
                     hp_clust = defP("hc", var = hpDir)))))


## ADAPTATION
## M$mu$mixin(adHST)
## M$mu$do.mc_scale <- T
## M$mu$do.ad.group_scale <- T
## M$tau$mixin(adHST)
## M$tau$do.mc_scale <- T
## update(M, 300)

## par(mfrow = c(2, 2))
## matplot(drop(M$mu.$mc_st)[, ], type = "l")
## matplot(drop(M$mu.$mc_st)[-burnin, ], type = "l")
## matplot(drop(M$tau.$mc_st)[-burnin], type = "l")
## matplot(drop(M$pr_clust.$mc_st)[], type = "l")


test_that("normal mixture works", {
    M$clust$st <- sample(1:2, N, T)
    update(M, 1000)

    par(mfrow = c(2, 2))
    matplot(drop(M$mu.$mc_st)[, ], type = "l")
    matplot(drop(M$mu.$mc_st)[-burnin, ], type = "l")
    matplot(drop(M$tau.$mc_st)[-burnin], type = "l")
    matplot(drop(M$pr_clust.$mc_st)[], type = "l")

    est <- sort(colMeans(drop(M$pr_clust$mc_st)[-burnin, ]))
    actual <- sort(table(trueCL)/N)
    rbind(est, actual)
    expect_close(est[[1]], actual[[1]], .01)
    
    est <- sort(colMeans(drop(M$mu.$mc_st)[-burnin, ]))
    actual <- sort(tapply(Y, trueCL, mean))
    rbind(est, actual)
    expect_close(est[[1]], actual[[1]], 1)
    expect_close(est[[2]], actual[[2]], 1)

    expect_close(mean(M$tau.$mc_st[-burnin]), 1/trueSD^2, .0005)
})

### jags
## ## Must have at least one data point with fixed assignment 
## ## to each cluster, otherwise some clusters will end up empty:
## ## clust <- rep(NA, N)
## ## clust[c(which.min(Y), which.max(Y))] <- 1:2

## library(rjags)
## js <- "model {
##     ## Likelihood:
##     for( i in 1 : N ) {
##         Y[i] ~ dnorm( mu[i] , tau ) 
##         mu[i] <- MU[ clust[i] ]
##         clust[i] ~ dcat( P )
##     }
##     ## Prior:
##     tau ~ dlnorm( 1 , .005 )
##     for ( k in 1:K ) {
##         MU[k] ~ dnorm( 0 , .005 )
##     }
##     P ~ ddirch( hpDir )
## }"

## jm <- jags.model(textConnection(js))
## update(jm, 1000)
## out <- coda.samples(jm, c("P", "MU"), 2000)
## plot(out)
## colMeans(out[[1]])
## table(trueCL)/N
