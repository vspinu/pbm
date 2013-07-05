## VIRTUALS

pdMultivar <- mixin(
    initForms = list(
        init.R.validate.lldim = form(
            if(ncol(ll) > 1L)
            stop.pbm(sprintf("ncol of ll (%s) in multivariate cells should be 1", ncol(ll))))
        ),
    expr = quote(.lldim <- 1L), 
    subtype = "multi")


### Dirichlet
dirichlet_log <- function(x, alpha)
    ## x and alpha should be matrixes with same ncols and alpha is replicated
    ## over nrows
    rowSums((alpha - 1) * log(x)) + lgamma(rowSums(alpha)) - rowSums(lgamma(alpha))

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

pdDirichlet <- mixin(
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
        setrand.st = form(
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
    parent_mixins = pdMultivar,
    subtype = "Dirichlet")

pdDirichlet_conj <- mixin(
    setForms = list(
        set.st = form(
            ## todo: this is conjugate set.st for categorical child only. Split out!
            for(i in seq_len(size)){
                ._curix <- child[["ixs"]][[1]] == i
                ._counts <- tabulate(child[["st"]][._curix], nbins = child[["N"]][i])
                st[i, ] <- dirichlet_rng(1, ._counts + PV("theta")[i, ])
            }), 
        init.M.validate.child_type = form(
            if(!protoIs(child, "Categorical"))
            stop.pbm("Child of '", .type, "' should be of type 'Categorical'"))),
    parent_mixins = pdDirichlet,
    subtype = "conj")

PBM$initCells(defBC(type = "ddirich",
                    prototype = "conj.uc",
                    distr = pdDirichlet_conj))


### Categorical
pdCategorical <- mixin(
    ## todo: make N also a stochastic parameter.
    setFields = list(multiparnames = "P"),
    setForms = list(
        set.ll = form({
            ## tothink: wouldn't be more convenient to parametrize in terms of log?
            ll[] <- log(pv("P")[cbind(pix("P"), c(st))]) # matrix indexing
        }),
        setrand.st = form(
            for(i in seq_along(levels(ixs[[1L]][[1L]]))){
                st[ixs[[1]] == i] <-
                    sample(1:N[i], size, replace = TRUE,
                           prob = pv("P")[i, ])
            })),
    initForms = list(
        ## FIXME: TODO: add setFormMaybe to handle cases when parent forms are
        ## set in mixins. Or (better?) make setForms install forms only if
        ## already installed, with no error!
        init.R.build_min_max = form({
            min <- rep_len(1L, length(st))
            max <- rep_len(N, length(st))
        }),
        init.C.build.integer_st = form({
            st[] <- round(st)
            storage.mode(st) <- "integer"
        }),
        init.C.validate.N = form({
            if(max(st) > N) stop.pbm(sprintf("maximum in ST (%s) excedes the declared categorical dimmension N (%s)", max(st), N))
            if(min(st) < 1) stop.pbm(sprintf("minimum in ST (%s) is less than 1. Categorical varialbe takes values in 1:N.", min(st)))
        })),
    initFields = list(N = 1),
    subtype = "Categorical")

PBM$initCells(defBC(type = "dcat",
                    prototype = "discr.unif",
                    distr = pdCategorical))


### Normal
pdNormal <- mixin(
    setForms = list(
        set.ll = form(
            ll[] <- dnorm(st,
                          mean = PV("mean"),
                          sd=1/sqrt(PV("tau")),
                          log=TRUE)),
        setrand.st = form(
            st[] <- rnorm(size,
                          mean = PV("mean"), 
                          sd = 1/sqrt(PV("tau"))))),
    setFields = list(
        parnames = c("mean", "tau")),
    subtype = "Normal")

PBM$initCells(defBC(type = "dnorm",
                    prototype = "norm.MHrw.like.uc",
                    distr = pdNormal))



###_ Gamma
pdGamma <- mixin(
    setForms = list(
        set.ll = form(
            ll[] <- dgamma(st,
                           shape = PV("shape"), 
                           rate = PV("rate"), 
                           log = TRUE)),
        setrand.st = form(
            st[] <- rgamma(length(st),
                           shape = PV("shape"), 
                           rate = PV("rate")))), 
    setFields = list(
        parnames = c("shape", "rate")),
    subtype = "Gamma")

PBM$initCells(defBC(type = "dgamma",
                    prototype = "lnorm.MHrw.like.uc",
                    distr = pdGamma))


###_ LogNorm
pdLogNormal <- mixin(
    setForms = list(
        set.ll = form(
            ll[] <- dlnorm(st,
                           meanlog = PV("meanlog"),
                           sdlog=1/sqrt(PV("taulog")),
                           log=TRUE)),
        setrand.st = quote(
            st[] <- rlnorm(length(st),
                           meanlog = PV("meanlog"),
                           sdlog = 1/sqrt(PV("taulog"))))),
    setFields = list(
        parnames = c("meanlog", "taulog")),
    subtype = "LogNormal")

PBM$initCells(defBC(type = "dlnorm",
                    prototype = "lnorm.MHrw.like.uc",
                    distr = pdLogNormal))



### GaNorm conjugate
pdGaNormal <- mixin(
    initFields = list(TR = function(x) x),
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
        set.ll = expression({
            for(i in seq_len(size))
                ll[[i]] <- dgamma(st[[i, 2L]],
                                  shape = pv("alpha0")[i, ],
                                  rate = pv("beta0")[i, ], log=T) +
                                      dnorm(st[[i, 1L]], mean = pv("mu0")[i, ],
                                            sd = 1/sqrt(st[[i, 2L]]*pv("n0")[i, ]), log=T)
        }),
        setrand.st = expression(stop("GaNorm randset is not implemented"))),
    initForms = list(
        init.R.build.nr_c_grs = expression({## nr elements in each group as given by child[["ixs"]][[1L]]
            ## uses children! should be in very late stage! after M.build is done!
            nr_c_grs <- c(tapply(get("st", child), child[["ixs"]][[1L]], length))
            names(nr_c_grs) <- NULL}),
        init.M.build.posterior = form(
            posterior <- parents[[1]][["st"]]), 
        set.st0 = expression({
            sum_c_grs <- rowsum(TR(c(child[["st"]])), group=child[["ixs"]][[1L]])
            mean_c_grs <- c(sum_c_grs/nr_c_grs)
            posterior[, "mu0"] <- (sum_c_grs + pv("mu0") * pv("n0"))/ (nr_c_grs + pv("n0"))
            posterior[, "n0"] <- nr_c_grs + pv("n0")
            posterior[, "alpha0"] <- pv("alpha0") + nr_c_grs/2
            posterior[, "beta0"] <- pv("beta0") +
                (rowsum.default((TR(c(get("st", child))) -
                                 mean_c_grs[child[["ixs"]][[1L]]])^2, child[["ixs"]][[1L]]) +
                 (nr_c_grs * pv("n0") * (mean_c_grs - pv("mu0"))^2)/(nr_c_grs + pv("n0")))/2
            for(i in seq_len(size)){
                st[[i, 2L]] <- rgamma(1L,
                                      shape=posterior[[i, "alpha0"]],
                                      rate=posterior[[i, "beta0"]])
                st[[i, 1L]] <- rnorm(1L,
                                     mean=posterior[[i, "mu0"]],
                                     sd=1/sqrt((st[[i, 2L]]*posterior[[i, "n0"]])))
            }})),
    subtype = "GaNormal")

PBM$initCells(defBC(type = "dganorm",
                    prototype = "conj.uc",
                    distr = pdGaNormal))



### LogGaNorm conjugate
## fixme: this should rely on mixin inheritance rather than cell inheritance
pdLogGaNormal <- mixin(setFields = list(
                           var = c( meanlog=0, taulog=1),
                           lldim = 1L, 
                           TR=function(x) log(x)),
                       expr = expression({
                           protocol$rv <-list(dim = c(2L),
                                              dimnames=list(c("taulog", "meanlog")))  #1L, 2L
                           protocol$st <- list(dim = c(NA, 2L),
                                               dimnames = list(NULL, c("taulog", "meanlog")))
                       }), 
                       subtype = "LogGaNormal")

PBM$initCells(defBC(type = "log",
                    prototype = "dganorm",
                    distr = pdLogGaNormal))



## ###_    * Unif
## PBM$initCells(defBC(type = "Unif", prototype = "unif.MHrw.like.uc",
##                     setForms = list(
##                         set.ll = quote(
##                             ll[] <- dunif(st,
##                                           min=parents[[1]][["st"]][, "min"][ix0],
##                                           max=parents[[1]][["st"]][, "max"][ix0],
##                                           log=TRUE)),
##                         setrand.st = quote(
##                             st[] <- runif(length(st),
##                                           min=parents[[1]][["st"]][, "min"][ix0],
##                                           max=parents[[1]][["st"]][, "max"][ix0]))),
##                     setFields = list(
##                         protocol = list(
##                             parents = list(
##                                 list(
##                                     dim=list(NULL, 2L),
##                                     dimnames=list(NULL, c("min", "max"))))
##                             ))))
