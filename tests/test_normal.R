### SIMPLE INFERENCE
### 
N <- 1000
sims <- 2000
mu <- 1
tau <- .25
Y <- rnorm(N, mu, sqrt(1/tau))
burnin <- 1:500



## KNOWN TAU
M <- pbm("NORMAL",
         DATA = defBC("dc.", mixin = pdNorm, 
             st = Y, 
             mu = defP("pd(Norm)",
                 cix = 1, size = 1, scale = .1, 
                 hp_mu = defP("hc",
                     var = c(mean = 0, tau = .00001))),
             tau = defP("hc",
                 cix = 1, st = tau)))

test_that("NORMAL with know TAU gives correct estimates",{
    update(M, nr_iter = sims)
    plot(M$mu$mc_st[-burnin,, ], type = "l")
    abline(h = mean(Y), col = "red", lwd = 2)
    expect_close(mean(Y), mean(M$mu.$mc_st[-burnin,, ]))
})


## KNOWN MU
M <- pbm("NORMAL",
         DATA = defBC("dc.", mixin = pdNorm, st = Y,
             mu =
             defP("hc", cix = 1, st = mu),
             tau =
             defP("pd(LogNorm)",
                  cix = 1, size = 1, scale = .1, 
                  hp_tau = defP("hc",
                      var = c(meanlog = 0, taulog = .005)))))

test_that("NORMAL with know MU gives correct estimates",{
    update(M, nr_iter = sims)
    plot(M$tau$mc_st[-burnin,, ], type = "l")
    abline(h = 1/var(Y), col = "red", lwd = 2)
    expect_close(mean(M$tau$mc_st[-burnin,, ]), 1/var(Y), .001)
})


## UNKNOWN MU, TAU
M <- pbm("NORMAL",
         DATA = defBC("dc.", mixin = pdNorm, 
             st = Y, 
             mu = defP("pd(Norm)",
                 cix = 1, size = 1, scale = .1, 
                 hp_mu = defP("hc",
                     var = c(mean = 0, tau = .00001))),
             tau =
             defP("pd(LogNorm)",
                  cix = 1, size = 1, scale = .1, 
                  hp_tau = defP("hc",
                      var = c(meanlog = 0, taulog = .005)))))


test_that("NORMAL with unknown MU and TAU gives correct estimates",{
    update(M, nr_iter = sims*2)

    par(mfrow = c(1, 2))
    plot(M$mu.$mc_st[-burnin,, ], type = "l")
    abline(h = mean(Y), col = "red", lwd = 2)
    expect_close(mean(Y), mean(M$mu.$mc_st[-burnin,, ]))

    plot(M$tau.$mc_st[-burnin,, ], type = "l")
    abline(h = 1/var(Y), col = "red", lwd = 2)
    expect_close(mean(M$tau.$mc_st[-burnin,, ]), 1/var(Y))
})

