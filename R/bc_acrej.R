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
                        mh_tr = tIdentity), 
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
                            st[] <- mh_tr@TR(mh_tr@ITR(`_st`) + c(rnorm(length(st)) * scale))),
                        set.alphas = form({
                            ## reminder: transforms are of the proposal not of the ll(st)!
                            alphas <- exp((ll_all - fastRS(mh_tr@LL(st)))  -
                                          (`_ll_all` - fastRS(mh_tr@LL(`_st`))))
                        }))))

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
                    ## these are sampler parameters, not distr parameters 
                    initFields = list(rwMin = 0, rwMax = 1), 
                    initForms = list(
                        init.R.build_rwMin_rwMax = form({
                            rwMin <- rep_len(rwMin, length(st))
                            rwMax <- rep_len(rwMax, length(st))
                        })),
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- runif(size, min = rwMin, max = rwMax)),
                        set.alphas = form({
                            alphas <- exp(ll_all  - `_ll_all`)
                        }))))

###_   + unif discrete
PBM$initCells(defBC(type = "discr", prototype="unif.mhrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            ## fixme: rwMin, rwMax are vectors like ix
                            ## because sample() is not vectorized we are stuck
                            st[] <- rdunif(rwMin, rwMax)
                        ))))
