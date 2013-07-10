sims <- 1500
N <- 2000
mulog <- 1
taulog <- .25
Y <- rlnorm(N, mulog, sqrt(1/taulog))
burnin <- 1:500
ustart <- 500
par(mfrow = c(1, 1))


## KNOWN taulog
M <- pbm("LogNormal",
         DATA = defBC("dc.", mixin = pdLogNormal, 
             st = Y, 
             mulog = defP("dnorm",
                 cix = 1, size = 1, scale = .5, 
                 hp_mulog = defP("hc",
                     var = c(mean = 0, tau = .00001))),
             taulog = defP("hc",
                 cix = 1, st = taulog)))

test_that("LogNormal with known TAULOG works", {
    update(M, nr_iter = sims)
    plot(M$mulog$mc_st[-burnin], type = "l")
    expect_close(mean(M$mulog$mc_st[-burnin]), mean(log(Y)), .01)
})


## KNOWN mulog
M1 <- pbm("LogNormal",
         DATA = defBC("dc", mixin = pdLogNormal,
             st = Y,
             mulog = defP("hc", cix = 1, st = mulog),
             taulog = defP("dlnorm", do.mc_ll = T, 
                 cix = 1, size = 1, scale = .1, 
                 hp_taulog = defP("hc",
                     var = c(meanlog = 0, taulog = .005)))))

test_that("LogNormal with known MULOG works", {
    update(M1, nr_iter = sims)
    plot(M1$taulog$mc_st[-burnin], type = "l")
    mn1 <<- mean(M1$taulog$mc_st[-burnin])
    expect_close(mean(M1$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})


## KNOWN mulog, with tranform cell
M2 <- pbm("LogNormal",
         DATA = defBC("dc.", mixin = pdLogNormal,
             st = Y,
             mulog = defP("hc", cix = 1, st = mulog),
             taulog = defP("tr", tr = tExp, do.mc_ll = T, 
                 tau = defP("dnorm", do.mc_ll = T, 
                     cix = 1, size = 1, scale = .1, 
                     hp_tau = defP("hc",
                         var = c(mean = 0, tau = .005))))))

test_that("LogNormal with known MULOG works", {
    update(M2, nr_iter = sims)
    plot(M2$taulog$mc_st[-burnin], type = "l")
    plot(M2$tau$mc_st[-burnin], type = "l")
    mn2 <<- mean(M2$taulog$mc_st[-burnin])
    expect_close(mean(M2$taulog$mc_st[-burnin]), 1/var(log(Y)), .001)
})

test_that("M1 (lnorm) and M2 (exp tr) give similar results", {
    layout(matrix(c(1, 3, 2, 3), 2))
    hist(log(st1 <- M1$taulog$mc_st))
    hist(log(st2 <- M2$taulog$mc_st))
    matplot(cbind(st1, st2)[-burnin, ], type = "l", lty = 1)
    expect_close(mn1, mn2, .001)
})

## M2$tau$do.debug <- T
## update(M2)

## c(M1$taulog$rejects, M2$tau$rejects)
## plot(window(M2$tau$mc_st, ustart))

## hist(log(st1 <- M1$taulog$mc_st))
## hist(log(st2 <- M2$taulog$mc_st))
## matplot(cbind(st1, st2)[-burnin, ], type = "l")
## qqplot(st1, st2)
## stdif <- (st1 - st2)[-burnin, ]
## mean(stdif)
## hist(stdif)
## ll1 <- M1$taulog$mc_ll[-burnin, ]
## ll2 <- (M2$tau$mc_ll - M2$tau$mc_st)[-burnin, ]
## lldif <- ll1-ll2
## hist(lldif)
## mean(lldif)
## plot(ts(lldif))



## UNKNOWN mulog, taulog
M <- pbm("LogNormal",
         DATA = defBC("dc.", mixin = pdLogNormal, 
             st = Y, 
             mulog = defP("dnorm",
                 cix = 1, size = 1, scale = .1, 
                 hp_mulog = defP("hc",
                     var = c(mean = 0, tau = .00001))),
             taulog =
             defP("dlnorm",
                  cix = 1, size = 1, scale = .1, 
                  hp_taulog = defP("hc",
                      var = c(meanlog = 0, taulog = .005)))))

test_that("LogNormal with unknown MULOG & TAULOG works", {

    update(M, nr_iter = 3*sims)

    par(mfrow = c(1, 2))
    plot(M$mulog.$mc_st[-burnin], type = "l")
    abline(h = mean(log(Y)), col = "red", lwd = 2)
    expect_close(mean(M$mulog.$mc_st[-burnin]), mean(log(Y)), .01)

    plot(M$taulog.$mc_st[-burnin], type = "l")
    abline(h = 1/var(log(Y)), col = "red", lwd = 2)
    expect_close(mean(M$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})
