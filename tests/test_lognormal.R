context("LogNorm")

sims <- 2000
N <- 1000
mulog <- 1
taulog <- .25
Y <- rlnorm(N, mulog, sqrt(1/taulog))
burnin <- 1:250
ustart <- 250
par(mfrow = c(1, 1))


test_that("LogNorm, known TAULOG", {

    M <- pbm("LogNorm",
             DATA = defBC("dc.", mixin = pdLogNorm, 
                 st = Y, 
                 mulog = defP("pd(Norm)",
                     cix = 1, size = 1, scale = .5, 
                     hp_mulog = defP("hc",
                         var = c(mean = 0, tau = .005))),
                 taulog = defP("hc",
                     cix = 1, st = taulog)))

    update(M, nr_iter = sims)
    plot(M$mulog$mc_st[-burnin], type = "l")
    expect_close(mean(M$mulog$mc_st[-burnin]), mean(log(Y)), .01)
})


test_that("LogNorm, known MULOG", {

    M <- M1 <- pbm("LogNorm",
                   DATA = defBC("dc", mixin = pdLogNorm,
                       st = Y,
                       mulog = defP("hc", cix = 1, st = mulog),
                       taulog = defP("pd(LogNorm)", do.mc_ll = T, 
                           cix = 1, size = 1, scale = .1, 
                           hp_taulog = defP("hc",
                               var = c(meanlog = 0, taulog = .005)))))

    update(M1, nr_iter = sims)
    plot(M1$taulog$mc_st[-burnin], type = "l")
    mn1 <<- mean(M1$taulog$mc_st[-burnin])
    expect_close(mean(M1$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})


## Known MULOG
M <- M1 <- pbm("LogNorm",
               DATA = defBC("dc", mixin = pdLogNorm,
                   st = Y,
                   mulog = defP("hc", cix = 1, st = mulog),
                   taulog = defP("norm.mhrw",
                       mixin = pdLogNorm,
                       mh_tr = tExp, 
                       do.mc_ll = T, 
                       cix = 1, size = 1, scale = .1, 
                       hp_taulog = defP("hc",
                           var = c(meanlog = 0, taulog = .005)))))

test_that("LogNorm with known MULOG works (dirrect specification)", {
    update(M1, nr_iter = sims)
    plot(M1$taulog$mc_st[-burnin], type = "l")
    mn1 <<- mean(M1$taulog$mc_st[-burnin])
    expect_close(mean(M1$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})


## Known MULOG (transform cell)
M <- M2 <-
    pbm("LogNorm",
        DATA = defBC("dc.", mixin = pdLogNorm,
            st = Y,
            mulog = defP("hc", cix = 1, st = mulog),
            taulog = defP("tr",
                st_tr = tExp, do.mc_ll = T, 
                tau = defP("pd(Norm)", do.mc_ll = T, 
                    cix = 1, size = 1, scale = .1, 
                    hp_tau = defP("hc",
                        var = c(mean = 0, tau = .005))))))

test_that("LogNorm with known MULOG works (with transform cell)", {
    update(M2, nr_iter = sims)
    plot(M2$taulog$mc_st[-burnin], type = "l")
    plot(M2$tau$mc_st[-burnin], type = "l")
    mn2 <<- mean(M2$taulog$mc_st[-burnin])
    expect_close(mean(M2$taulog$mc_st[-burnin]), 1/var(log(Y)), .001)
})

test_that("M1 (lnorm) and M2 (tExp) give similar results", {
    layout(matrix(c(1, 3, 2, 3), 2))
    hist(log(st1 <- M1$taulog$mc_st))
    hist(log(st2 <- M2$taulog$mc_st))
    matplot(cbind(st1, st2)[-burnin, ], type = "l", lty = 1)
    expect_close(mn1, mn2, .001)
    expect_close(M1$taulog$rejects, M2$tau$rejects, 50)
})



## Known MULOG, (1-to-1 internal transform)
M3 <- pbm("LogNorm",
          DATA = defBC("dc", mixin = pdLogNorm,
              st = Y,
              mulog = defP("hc", cix = 1, st = mulog),
              taulog = defP("norm.mhrw",
                  mixin = list(pdNorm, trExp), 
                  do.mc_ll = T, 
                  cix = 1, size = 1, scale = .1, 
                  hp_taulog = defP("hc",
                      var = c(mean = 0, tau = .005)))))

test_that("LogNorm with known MULOG (internal 1-to-1 transform)", {
    update(M3, nr_iter = sims)
    plot(M3$taulog$mc_st[-burnin], type = "l")
    mn3 <<- mean(M3$taulog$mc_st[-burnin])
    expect_close(mn1, mn3, .001)
    expect_close(M1$taulog$rejects, M3$taulog$rejects, 50)
    expect_close(mean(M3$taulog$mc_st[-burnin]), 1/var(log(Y)), .001)
})



test_that("LogNorm with unknown MULOG & TAULOG works", {

    M <- pbm("LogNorm",
             DATA = defBC("dc.", mixin = pdLogNorm, 
                 st = Y, 
                 mulog = defP("pd(Norm)",
                     cix = 1, size = 1, scale = .1, 
                     hp_mulog = defP("hc",
                         var = c(mean = 0, tau = .00001))),
                 taulog =
                 defP("pd(LogNorm)",
                      cix = 1, size = 1, scale = .1, 
                      hp_taulog = defP("hc",
                          var = c(meanlog = 0, taulog = .005)))))


    update(M, nr_iter = 3*sims)

    par(mfrow = c(1, 2))
    plot(M$mulog.$mc_st[-burnin], type = "l")
    abline(h = mean(log(Y)), col = "red", lwd = 2)
    expect_close(mean(M$mulog.$mc_st[-burnin]), mean(log(Y)), .01)

    plot(M$taulog.$mc_st[-burnin], type = "l")
    abline(h = 1/var(log(Y)), col = "red", lwd = 2)
    expect_close(mean(M$taulog.$mc_st[-burnin]), 1/var(log(Y)), .001)
})
