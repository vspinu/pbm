N <- 1000
mulog <- 1
taulog <- .25
Y <- rlnorm(N, mulog, sqrt(1/taulog))
burnin <- 1:500


## KNOWN taulog
M <- pbm("LogNormal",
         DATA = defBC("dc.", distr = pdLogNormal, 
             st = Y, 
             mulog = defP("dnorm",
                 cix = 1, size = 1, scale = .5, 
                 hp_mulog = defP("hc",
                     var = c(mean = 0, tau = .00001))),
             taulog = defP("hc",
                 cix = 1, st = taulog)))

test_that("LogNormal with known TAULOG works", {
    update(M, nr_iter = 1000)
    plot(M$mulog$mc_st[-burnin], type = "l")
    expect_close(mean(M$mulog$mc_st[-burnin]), mean(log(Y)), .01)
})


## KNOWN mulog
M1 <- pbm("LogNormal",
         DATA = defBC("dc", distr = pdLogNormal,
             st = Y,
             mulog = defP("hc", cix = 1, st = mulog),
             taulog = defP("dlnorm",
                 cix = 1, size = 1, scale = .1, 
                 hp_taulog = defP("hc",
                     var = c(meanlog = 0, taulog = .005)))))

test_that("LogNormal with known MULOG works", {
    update(M1, nr_iter = 1000)
    plot(M1$taulog$mc_st[-burnin], type = "l")
    mn1 <<- mean(M1$taulog$mc_st[-burnin])
    expect_close(mean(M1$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})


## KNOWN mulog, with tranform
M2 <- pbm("LogNormal",
         DATA = defBC("dc.", distr = pdLogNormal,
             st = Y,
             mulog = defP("hc", cix = 1, st = mulog),
             taulog = defP("tr", tr = tExp, 
                 tau = defP("dnorm",
                     cix = 1, size = 1, scale = .1, 
                     hp_tau = defP("hc",
                         var = c(mean = 0, tau = .005))))))

test_that("LogNormal with known MULOG works", {
    update(M2, nr_iter = 1000)
    plot(M2$taulog$mc_st[-burnin], type = "l")
    mn2 <<- mean(M2$taulog$mc_st[-burnin])
    expect_close(mean(M2$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})


out <- replicate(50, {
    update(M1, 1000, reinit = T)
    mn1 <- mean(M1$taulog$mc_st[-burnin])
    update(M2, 1000, reinit = T)
    mn2 <- mean(M2$taulog$mc_st[-burnin])
    mn1-mn2
})


## UNKNOWN mulog, taulog
M <- pbm("LogNormal",
         DATA = defBC("dc.", distr = pdLogNormal, 
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
    update(M, nr_iter = 3000)

    par(mfrow = c(1, 2))
    plot(M$mulog.$mc_st[-burnin], type = "l")
    abline(h = mean(log(Y)), col = "red", lwd = 2)
    expect_close(mean(M$mulog.$mc_st[-burnin]), mean(log(Y)), .01)

    plot(M$taulog.$mc_st[-burnin], type = "l")
    abline(h = 1/var(log(Y)), col = "red", lwd = 2)
    expect_close(mean(M$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})
