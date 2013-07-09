### SIMPLE INFERENCE
### Dirichlet conjugate priors
N <- 1000
a <- (1:3)/6
apr <- c(1, 1, 1)
D <- dirichlet_rng(N, a)
Y <- apply(D, 1, function(p) sample(1:3, 1, T, prob = p))

M <- pbm("CatDirichlet",
         DATA = defBC("dc.", distr = pdCategorical,
           st = Y, N = 3, 
           D = defP("ddirich",
             varsize = 3, cix = rep.int(1, N), cix_dim = 1, 
             hpD = defP("hc",
               cix = 1, var = apr))))

M$D$do.mc_st <- T
update(M, nr_iter = 1000)
matplot(M$D.$mc_st, type = "l")

test_that("CatDirichlet model returns correct posteriors",{
          post <- unname(colMeans(M$D.$mc_st[-(1:100), ]))
          actual <- tabulate(Y)
          actual <- actual/sum(actual)
          rbind(actual, post)
          expect_close(post, actual, .003)
      })

#### JAGS
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


