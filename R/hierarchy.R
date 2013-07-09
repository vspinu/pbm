## Sources("/Homes/vitoshka/Dropbox/works/R_dev/Classes/protoClasses.R")
## source("/home/vitoshka/Dropbox/works/R_dev/pbm/funcs.R")
## library(coda)
## library(Rgraphviz)
## load_all("~/works/protoClasses/")
source("~/works/pbm/R/utils.R")
source("~/works/pbm/R/fields.R")
source("~/works/pbm/R/sync.R")
source("~/works/pbm/R/valibuild.R")
source("~/works/pbm/R/S4.R")
source("~/works/pbm/R/mirror.R")
source("~/works/pbm/R/predict.R")
source("~/works/pbm/R/update.R")
source("~/works/pbm/R/transforms.R")
source("~/works/pbm/R/pbm.R")
source("~/works/pbm/R/mcmc.R")           

### HC
PBM$initCells(defBC(type = "hc", prototype="*",
                    ## do.mc_st=FALSE,
                    ## do.st=FALSE,
                    setFields = list(
                        do.update = FALSE, ## fixme: doesn't take effect, so need to list all of them
                        do.st = FALSE, 
                        do.ll = FALSE,
                        do.mc_ll = FALSE,
                        do.mc_ll = FALSE,
                        do.pc_st = FALSE,
                        do.pc_ll = FALSE),
                    initForms = list(
                        init.C.build.mc_ll = NULL,
                        init.C.build.mc_st = NULL,
                        init.C.build._ll   = NULL,
                        init.M.validate = NULL,
                        init.R = NULL),
                    expr = expression({
                        ## .basic_foldable_objects <- c("st", "ll")
                        .basic_foldable_inxs <- c()
                        init.M.build <- init.M.build[c("st", "pos_in_C")]
                    })))

###_* DC
PBM$initCells(defBC(type = "dc", prototype = "*",
                    setFields = list(
                        do.mc_st = FALSE,
                        do.st = FALSE,
                        do.pc_st = FALSE,
                        do.pc_ll = TRUE),
                    expr = expression({
                        init.M.build[c("pos_in_C")] <- NULL
                    })))

###_* UC
PBM$initCells(defBC(type = "uc", prototype = "*"))

###_ + acrej likelihood cells
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

PBM$initCells(defBC(type = "MHrw",
                    prototype = "acrej.uc",
                    setForms = list(
                        set.st = form({
                            e(update.ll)
                            .update.ll(children)
                            `_ll` <- ll
                            `_ll_all` <-
                                fastRS(ll) + fastRS(.get.ll(children, .type))
                            e(set.st_proposal) #!
                            e(set.ll_is_old)
                            .set.ll_is_old(children)
                            e(update.ll)
                            .update.ll(children)
                            ll_all <- fastRS(ll) + fastRS(.get.ll(children, .type))
                            e(set.alphas)
                            u <- runif(n=size, min=0, max=1)
                            rejects <- (u > alphas) | is.na(alphas)
                            accepts <- !rejects
                            .rej <- .rej+rejects
                            e(resync.st) ## UPDATE.main marks ..ll_is_old.. and ..st_is_old.. 
                        }),
                        UPDATE.main.ll = form(
                            ## Overwrite the updating routine to take advantage
                            ## of rejects. Do not confuse with update.ll
                            if(do.ll){
                                e(resync.ll)
                                e(set.ll_is_updated)
                                .resync.ll(children, rejects, .type)
                                .set.ll_is_updated(children)
                                })),
                    initForms = list(
                        set.st_proposal = form(
                            doc = "API: should set ST with appropriate given the old state found in _st.", 
                            stop.pbm(" set.st_proposal is not specified.")),
                        set.alphas = form(
                            doc = "API: should set the MH alphas", 
                            stop.pbm(" set.alphas is not specified."))
                        )))

PBM$initCells(defBC(type = "norm", prototype="MHrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- TR(ITR(`_st`) + c(rnorm(length(st)) * scale))),
                        set.alphas = form({
                            ## instr_lratios <- rowSums(log(st)-log(`_st`))
                            alphas <- exp(ll_all  - `_ll_all`)
                        })),
                    expr = expression(TR <- identity, ITR <- identity)))

PBM$initCells(defBC(type = "lnorm", prototype="MHrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- exp(log(`_st`) + c(rnorm(length(st)) * scale))),
                        set.alphas = form({
                            instr_lratios <- rowSums(log(st)-log(`_st`))
                            alphas <- exp(ll_all  - `_ll_all` + instr_lratios)
                        }))))

###_   + unif
PBM$initCells(defBC(type = "unif", prototype="MHrw.acrej.uc",
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
PBM$initCells(defBC(type = "discr", prototype="unif.MHrw.acrej.uc",
                    setForms = list(
                        set.st_proposal = form(
                            for(i in 1:varsize){
                                        #fixme: min is a vector like ix?
                                st[, i] <- sample(min[i]:max[i], size, replace = TRUE)
                            }))))

###_ + CONJ  (conjugate cells)
PBM$initCells(defBC(type = "conj", prototype = "uc",
                    initForms = list(
                        init.R.build_chITR = form(
                            chITR <- children[["i"]]$ITR),
                        init.M.validate.child_type = form())))

PBM$initCells(defBC(type = "tr", prototype = "*",
                    ## API: NOTE: tr cell, passes likelihood computation to next
                    ## cell. Examples: .update.ll, .set.ll_is_old, .get.ll
                    
                    ## simple one to one transformation
                    setFields = list(
                        ## tr cells are updated by parent cells in UPDATE.tr_children
                        do.update = FALSE), 
                    initFields = list(
                        tr = protoField(.field.tr,
                            doc = "API: this is an object of protoTransform class describing the state and likelihood transformation.")),
                    initForms = list(
                        init.R.build.st = form(e(init.M.build.st))), 
                    setForms = list(
                        set.st = form(st[] <- tr@TR(pst(1))),
                        set.ll = form(ll[] <- tr@LL(st)),
                        set.rand.st = form(e(set.st)),
                        init.M.build.st = form({
                            ## preclude the transform from droping dimensions 
                            .fields$st <- pst(1)
                            st[] <- tr@TR(st)
                        }),
                        update.ll = form(
                            if(..ll_is_old..){
                                e(update.st) # does nothing if already set
                                e(set.ll)
                                e(set.ll_is_updated)
                            })),
                    setMethods = list(
                        .resync.ll = ## fixme: st resync is meaningful here as well
                        function(rejects, ptype, token = runif(1L)){
                            ## next method; when implemented remove this:
                            if(!identical(..ctoken.., token)){
                                ix_rejects <- rejects[pix(, ptype)]
                                ix_accepts <- !ix_rejects
                                .self[["ll"]][ix_rejects] <- .self[["_ll"]][ix_rejects]
                                .self[["_ll"]][ix_accepts] <- .self[["ll"]][ix_accepts]
                                assign("..ctoken..", token, envir = .self)
                            }
                            .resync.ll(children, rejects, .type, token)
                        },
                        .get.ll =
                        function(ptype, token = runif(1L)){
                            ## return a matrix of this ll and parents LLs
                            structure(cbind(fastRS(ll), .get.ll(children, .type, token)),
                                      dimnames = list(NULL, c(.type, names(children))))
                        },
                        .update.ll = function() {
                            e(update.ll)
                            .update.ll(children)
                        },
                        .ll_is_old = function(){
                            e(set.ll_is_old)
                            .set.ll_is_old(children)
                        } ,
                        .set.ll_is_updated = function(){
                            e(set.ll_is_updated)
                            .set.ll_is_updated(children)
                        }), 
                    expr = expression(
                        tr <- tIdentity)))


source("~/works/pbm/R/distributions.R")
source("~/works/pbm/R/adapt.R")


