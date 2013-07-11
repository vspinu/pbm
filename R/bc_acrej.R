### Accept reject likelihood cells

PBM$initCells(defBC(type = "acrej", prototype="uc",
                    initFields = list(
                        adapt = protoField(.field.adapt),
                        scale = protoField(.field.scale),
                        rejects = protoReadOnlyField(".rej")),
                    expr = expression({
                        .rej <- rep.int(0, size)
                        scale <- 1
                        .basic_foldable_objects <- c(.basic_foldable_objects, "scale")
                    }), ## default
                    initForms = list(
                        resync = form(
                            ## tothink: move up the hierarchy?
                            doc = "Resync forms assume that objects 'rejects' and 'accepts are present and have been set.", 
                            ll = form({
                                ll[rejects] <- `_ll`[rejects]
                                `_ll`[accepts] <- ll[accepts]
                            }),
                            st = form(
                                st[rejects] <- `_st`[rejects], 
                                `_st`[accepts] <- st[accepts]
                                )),
                        init.M.build = form(
                            rejects_accepts = form({
                                rejects <- accepts <- vector(length=size)
                                .rej <- vector("integer", length=size)
                            }),
                            scale = form(.self$scale <- scale) ## replicate scale
                            ))))

### MHRW (random walk)

PBM$initCells(defBC(type = "mhrw",
                    prototype = "acrej.uc",
                    initFields = list(
                        mhtr = tIdentity), 
                    setForms = list(
                        set.st = form({
                            e(update.ll);	 .update.ll(children)
                            `_ll` <- ll
                            `_ll_all` <-
                                fastRS(ll) + fastRS(.get.ll(children, .type))
                            e(set.st_proposal)
                            e(set.ll_is_old); 	.set.ll_is_old(children)
                            e(update.ll); 	.update.ll(children)
                            ll_all <-
                                fastRS(ll) + fastRS(.get.ll(children, .type))
                            e(set.alphas)
                            u <- runif(n=size, min=0, max=1)
                            rejects <- (u > alphas) | is.na(alphas)
                            accepts <- !rejects
                            .rej <- .rej+rejects
                            e(resync.st)
                            ## e(update.ll) set it to F, reset
                            e(set.ll_is_old); 	.set.ll_is_old(children)
                        }),
                        UPDATE.main.ll = form({
                            ## Take advantage of rejects. Do not confuse with update.ll
                            e(resync.ll); 	.resync.ll(children, rejects, .type)
                            e(set.ll_is_updated); .set.ll_is_updated(children)
                        })),
                    initForms = list(
                        set.st_proposal = form(
                            doc = "API: should set ST with appropriate given the old state found in _st.", 
                            stop.pbm(" set.st_proposal is not specified.")),
                        set.alphas = form(
                            doc = "API: should set the MH alphas", 
                            stop.pbm(" set.alphas is not specified."))
                        )))

PBM$initCells(defBC(type = "norm", prototype="mhrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- mhtr@TR(mhtr@ITR(`_st`) + c(rnorm(length(st)) * scale))),
                        set.alphas = form({
                            alphas <- exp((ll_all - mhtr@LL(st))  -
                                          (`_ll_all` - mhtr@LL(`_st`)))
                        })),
                    expr = expression(TR <- identity, ITR <- identity)))

## PBM$initCells(defBC(type = "lnorm", prototype="mhrw.acrej.uc",
##                     setForms = list(
##                         set.st_proposal = form(
##                             st[] <- exp(log(`_st`) + c(rnorm(length(st)) * scale))),
##                         set.alphas = form({
##                             instr_lratios <- rowSums(log(st)-log(`_st`))
##                             alphas <- exp(ll_all  - `_ll_all` + instr_lratios)
##                         }))))

###_   + unif
PBM$initCells(defBC(type = "unif", prototype="mhrw.acrej.uc",
                    initFields = list(min = 0, max = 1), ## these are not distr parameters, but sampler parameters
                    initForms = list(
                        init.R.build_min_max = form({
                            min <- rep_len(min, length(st))
                            max <- rep_len(max, length(st))
                        })),
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- runif(size, min = min, max = max)),
                        set.alphas = form({
                            alphas <- exp(ll_all  - `_ll_all`)
                        }))))

###_   + unif discrete
PBM$initCells(defBC(type = "discr", prototype="unif.mhrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            for(i in 1:varsize){
                                        #fixme: min is a vector like ix?
                                st[, i] <- sample(min[i]:max[i], size, replace = TRUE)
                            }))))
