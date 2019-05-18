## VIRTUALS

pdMultivar <- mixin(
    initForms = list(
        init.R.validate.lldim = form(
            if(ncol(ll) > 1L)
            stop.pbm(sprintf("ncol of ll (%s) in multivariate cells should be 1", ncol(ll))))
        ),
    expr = quote(.lldim <- 1L))


### Dirichlet
dirichlet_log <- function(x, alpha) {
    ## x and alpha should be matrixes with same ncols and alpha is replicated
    ## over nrows
    rowSums((alpha - 1) * log(x)) + lgamma(rowSums(alpha)) - rowSums(lgamma(alpha))
}

dirichlet_rng <- function(n, alpha){
    if (is.matrix(alpha)) {
        k <- dim(alpha)[2]
        alpha <- c(t(alpha))
    } else
        k <- length(alpha)
    rand <- matrix(rgamma(n * k, shape = alpha), n, k, byrow = TRUE)
    rand/rowSums(rand)
}

## mt <- matrix(c(1, 2, 3, 10, 5, 1), 2)
## rng <- replicate(1000, dirichlet_rng(5, mt))
## rowMeans(rng, dims = 2)
## mt/rowSums(mt)

pdDirich <- mixin(
    setFields = list(multiparnames = "theta"),
    setForms = list(
        set.ll = form(
            ll[] <- dirichlet_log(st, PV("theta"))),
        ## do not delete! test the speed first
        ## set.ll = form({
        ##     for(i in seq_along(levels(ixs[[1]][[1L]]))){
        ##         ._curix <- ixs[[1L]] == i
        ##         ll[._curix] <- dirichlet_log(st[._curix], pv("theta")[i, ])
        ##     }}),
        set.rand.st = form(
            st[] <- dirichlet_rng(size, PV("theta"))
            )),
    initForms = list(
        init.C.normalize_st = form({
            st <- st/rowSums(st)
        }),
        init.C.validate.child_range = form(), ## todo
        init.R.validate.parent_varsize = form(
            if(pcell("theta")[["varsize"]] != varsize){
                stop.pbm(sprintf("varsize of the parent (%s) should be equal to varsize (%s) in Dirichlet cell",
                                 parents[[1]][["varsize"]], varsize))
            })),
    parentMixins = pdMultivar,
    subtype = "Dirich")

PBM$initCells(defBC(type = "pd", mixin = pdDirich,
                    prototype = "mhrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- dirichlet_rng(size, rep.int(1, varsize))),
                        set.alphas = form({
                            alphas <- exp(ll_all  - `_ll_all`)
                        }))))

pdDirichlet_conj <- mixin(
    setForms = list(
        set.st = form(
            ## todo: this is conjugate set.st for categorical child only. Split out!
            for(i in seq_len(size)){
                ._curix <- child[["ixs"]][[1]] == i
                ## fixme: limitation of matrix ST! Multivar categorical child
                ## might have different sizes N, which in turn requires
                ## different sizes of dirichlet prior. Potential solution: Use
                ## list, in matrix form. That to say: no multipar states. All
                ## multipars should be stored in lists?
                ._counts <- tabulate(child$st[._curix], nbins = varsize) ## child[["N"]][i] ??
                st[i, ] <- dirichlet_rng(1, ._counts + PV("theta")[i, ])
            }), 
        init.M.validate.child_type = form(
            if(!protoIs(child, "Cat"))
            stop.pbm("Child of '", .type, "' should be categorical (aka of subtype '(Cat)')"))),
    parentMixins = pdDirich,
    subtype = "conj")

PBM$initCells(defBC(type = "pd", prototype = "conj.uc",
                    mixin = pdDirichlet_conj))


### Categorical
pdCat <- mixin(
    ## TODO: make N stochastic parameter
    setFields = list(multiparnames = "P"),
    setForms = list(
        set.ll = form({
            ## tothink: wouldn't be more convenient to parametrize in terms of log?
            ll[] <- log(pv("P")[cbind(pix("P"), c(st))]) # matrix indexing
        }),
        set.rand.st = form(
            for(i in seq_along(levels(ixs[[1L]][[1L]]))){
                st[ixs[[1]] == i] <-
                    ## fixme: size?
                    sample(1:N[i], size, replace = TRUE, prob = pv("P")[i, ])
            })),
    initForms = list(
        ## FIXME: TODO: add setFormMaybe to handle cases when parent forms are
        ## set in mixin. Or (better?) make setForms install forms only if
        ## already installed, with no error!
        init.R.build_rwMin_rwMax = form({
            rwMin <- rep_len(1L, length(st))
            rwMax <- rep_len(N, length(st))
        }),
        init.C.build.integer_st = form({
            st[] <- round(st)
            storage.mode(st) <- "integer"
        }),
        init.C.validate.N = form({
            ## fixme: N is a vector.
            if(max(st) > max(N)) stop.pbm(sprintf("maximum in ST (%s) excedes the declared categorical dimmension N (%s)", max(st), max(N)))
            if(min(st) < 1) stop.pbm(sprintf("minimum in ST (%s) is less than 1. Categorical varialbe takes values in 1:N.", min(st)))
        })),
    initFields = list(N = 1),
    subtype = "Cat")

PBM$initCells(defBC(type = "pd",
                    prototype = "discr.unif.mhrw.acrej.uc",
                    mixin = pdCat))


### Norm
pdNorm <- mixin(
    setForms = list(
        set.ll = form(
            ll[] <- dnorm(st,
                          mean = PV("mean"),
                          sd=1/sqrt(PV("tau")),
                          log=TRUE)),
        set.rand.st = form(
            st[] <- rnorm(size,
                          mean = PV("mean"), 
                          sd = 1/sqrt(PV("tau"))))),
    setFields = list(
        parnames = c("mean", "tau")),
    subtype = "Norm")

PBM$initCells(defBC(type = "pd", mixin = pdNorm, 
                    prototype = "norm.mhrw.acrej.uc"))



###_ Gamma
pdGamma <- mixin(
    setForms = list(
        set.ll = form(
            ll[] <- dgamma(st,
                           shape = PV("shape"), 
                           rate = PV("rate"), 
                           log = TRUE)),
        set.rand.st = form(
            st[] <- rgamma(length(st),
                           shape = PV("shape"), 
                           rate = PV("rate")))), 
    setFields = list(
        parnames = c("shape", "rate")),
    subtype = "Gamma")

PBM$initCells(defBC(type = "pd", mixin = pdGamma,
                    prototype = "norm.mhrw.acrej.uc",
                    mh_tr = tExp))


###_ LogNorm
pdLogNorm <- mixin(
    setForms = list(
        set.ll = form(
            ll[] <- dlnorm(st,
                           meanlog = PV("meanlog"),
                           sdlog = 1/sqrt(PV("taulog")),
                           log = TRUE)),
        set.rand.st = quote(
            st[] <- rlnorm(length(st),
                           meanlog = PV("meanlog"),
                           sdlog = 1/sqrt(PV("taulog"))))),
    setFields = list(
        parnames = c("meanlog", "taulog")),
    subtype = "LogNorm")

PBM$initCells(defBC(type = "pd", mixin = pdLogNorm,
                    prototype = "norm.mhrw.acrej.uc",
                    mh_tr = tExp))


### Discrete Uniform
pdDiscrUnif <- mixin(
    setFields = list(
        parnames = c("min", "max")), 
    setForms = list(
        set.ll = form(
            ll[] <- log(ddunif(PV("min"), PV("max")))), 
        set.rand.st = form(
            st[] <- rdunif(PV("min"), PV("max")))
    ),
    subtype = "DiscrUnif")

PBM$initCells(defBC(type = "pd", mixin = pdDiscrUnif,
                    prototype = "mhrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- rdunif(PV("min"), PV("max"))), 
                        set.alphas = form(
                            alphas <- exp(ll_all  - `_ll_all`)))))


### GaNorm conjugate
pdGaNorm <- mixin(
    setFields = list(
        var = c(mu=0, tau=1),
        parnames = c("mu0", "n0", "alpha0", "beta0"),
        lldim = 1, 
        ## fixme: thisprotocol is old style
        protocol = list(
            parents = list(
                list(
                    varsize = 4L,
                    varnames = c("mu0", "n0", "alpha0", "beta0"))))),
    setForms = list(
        set.ll = form({
            for(i in seq_len(size)){
                ll[[i]] <-
                    dnorm(st[[i, 1L]], mean = pv("mu0")[i, ],
                          sd = 1/sqrt(st[[i, 2L]]*pv("n0")[i, ]),
                          log=T) + 
                              dgamma(st[[i, 2L]],
                                     shape = pv("alpha0")[i, ],
                                     rate = pv("beta0")[i, ],
                                     log=T)
            }
        }),
        set.rand.st = form(
            stop("GaNorm randset is not implemented (see rganorm functon)")
            )),
    initForms = list(
        init.R.build.nr_c_grs = form({## nr elements in each group as given by child[["ixs"]][[1L]]
            ## uses children! should be in very late stage! after M.build is done!
            nr_c_grs <- c(tapply(child$st, child[["ixs"]][[1L]], length))
            names(nr_c_grs) <- NULL}),
        init.M.build.posterior = form(
            posterior <- parents[[1]]$st), 
        set.st = form({
            ._trcst <- chITR(c(child$st))
            sum_c_grs <- rowsum(._trcst, group=child[["ixs"]][[1L]])
            mean_c_grs <- c(sum_c_grs/nr_c_grs)
            posterior[, "mu0"] <- (sum_c_grs + pv("mu0") * pv("n0"))/ (nr_c_grs + pv("n0"))
            posterior[, "n0"] <- nr_c_grs + pv("n0")
            posterior[, "alpha0"] <- pv("alpha0") + nr_c_grs/2
            posterior[, "beta0"] <- pv("beta0") +
                (rowsum.default((._trcst -
                                 mean_c_grs[child[["ixs"]][[1L]]])^2, child[["ixs"]][[1L]]) +
                 (nr_c_grs * pv("n0") * (mean_c_grs - pv("mu0"))^2)/(nr_c_grs + pv("n0")))/2
            for(i in seq_len(size)){
                st[[i, 2L]] <- rgamma(1L,
                                      shape=posterior[[i, "alpha0"]],
                                      rate=posterior[[i, "beta0"]])
                st[[i, 1L]] <- rnorm(1L,
                                     mean=posterior[[i, "mu0"]],
                                     sd=1/sqrt((st[[i, 2L]]*posterior[[i, "n0"]])))
            }}), 
        init.R.validate.child = form(
            if(length(children) != 1L){
                stop.pbm("conjugate GaNorm cell accepts only one child; supplied ", length(children))
            }
            ## , 
            ## if(!protoIs(children[[1L]], "dnorm")){
            ##     stop.pbm("child of conjugate GaNorm cell must be of type dnorm")
            ## }
            )),
    expr = expression(chITR <- identity), 
    subtype = "GaNorm")

PBM$initCells(defBC(type = "pd", prototype = "conj.uc", mixin = pdGaNorm))

rganorm <- function(n, mu0, n0, alpha0, beta0){
    T <- rgamma(n, shape = alpha0, rate = beta0)
    X <- rnorm(n, mean = mu0, sd = 1/sqrt(n0*T))
    cbind(X = X, T = T)
}

## .meanLNorm <- function(x){
##     exp(x[, "meanlog"] + 1/(x[, "taulog"]*2))
## }

## summary_ganorm <- function(N = 100000, mu0 = -5, n0 = 1, alpha0 = .01, beta0 = .01) {
##   x <- rganorm(N, mu0 = mu0, n0 = n0, alpha0 = alpha0, beta0 = beta0)
##   colnames(x) <- c("meanlog", "taulog")
##   colMeans(x)
##   tt <- .meanLNorm(x)
##   hist(tt[tt <10], 100)
##   str(list(median = median(tt, na.rm = T), 
##            mean = mean(tt, na.rm = T), 
##            sd = sd(tt, na.rm = T), 
##            median_fin = median(tt[is.finite(tt)], na.rm = T), 
##            mean_fin = mean(tt[is.finite(tt)], na.rm = T),
##            sd_fin = sd(tt[is.finite(tt)], na.rm = T), 
##            finite = as.list(table(is.finite(tt)))))
## }

## summary_ganorm(100000, mu0 = -1, n0 = 10, alpha0 = 5, beta0 = 10)
## summary_ganorm(100000, mu0 = -1, n0 = 10, alpha0 = 100, beta0 = 200)
## summary_ganorm(100000, mu0 = -5, n0 = 1, alpha0 = .01, beta0 = .01)
## summary_ganorm(100000, mu0=-5, n0=.1, alpha0=.01, beta0=.01)



### LogGaNorm conjugate
## fixme: this should rely on mixin inheritance rather than cell inheritance
pdLogGaNorm <- mixin(setFields = list(
                           var = c(meanlog=0, taulog=1),
                           lldim = 1L), 
                       expr = expression({
                           chITR <- log
                           ## chITR <- identity
                           protocol$rv <-list(dim = c(2L),
                                              dimnames=list(c("taulog", "meanlog")))  #1L, 2L
                           protocol$st <- list(dim = c(NA, 2L),
                                               dimnames = list(NULL, c("taulog", "meanlog")))
                       }), 
                       subtype = "LogGaNorm")

PBM$initCells(defBC(type = "pd", prototype = "pd(GaNorm)", mixin = pdLogGaNorm))



## Unif
## PBM$initCells(defBC(type = "Unif", prototype = "unif.mhrw.acrej.uc",
##                     setForms = list(
##                         set.ll = quote(
##                             ll[] <- dunif(st,
##                                           min=parents[[1]]$st[, "min"][ix0],
##                                           max=parents[[1]]$st[, "max"][ix0],
##                                           log=TRUE)),
##                         set.rand.st = quote(
##                             st[] <- runif(length(st),
##                                           min=parents[[1]]$st[, "min"][ix0],
##                                           max=parents[[1]]$st[, "max"][ix0]))),
##                     setFields = list(
##                         protocol = list(
##                             parents = list(
##                                 list(
##                                     dim=list(NULL, 2L),
##                                     dimnames=list(NULL, c("min", "max"))))
##                             ))))

