## Sources("/Homes/vitoshka/Dropbox/works/R_dev/Classes/protoClasses.R")
## source("/home/vitoshka/Dropbox/works/R_dev/pbm/funcs.R")
## library(coda)
## library(Rgraphviz)
library(protoClasses)
source("~/works/pbm/R/fields.R")
source("~/works/pbm/R/sync.R")
source("~/works/pbm/R/valibuild.R")
source("~/works/pbm/R/S4.R")
source("~/works/pbm/R/mirror.R")
source("~/works/pbm/R/predict.R")
source("~/works/pbm/R/update.R")

setContextClass("PBM", cellClass="BC")
    setMethod("initialize", "PBM",
          function(.Object, type = "--", prototype = NULL, ...){
              if(!is.null(prototype)){
                  cur_mir <- prototype$mirror
                  ## inherit .cells from base not from mirrors
                  prototype$mirror <- "base"
                  on.exit(prototype$mirror <- cur_mir)
              }

              pbm <- callNextMethod(.Object, type = type, prototype = prototype, ...)

              pbm$evalq({
                  current_mirror <- "base"
                  folds <- NULL
                  folds_names <- NULL
                  .cells <- as(.cells, "bCellContainer")
                  mirrors <- list(base = .cells)
              })
              pbm
          })

## don't delete
## setAs("PBM","graphNEL", .as_PBM_graphNEL)
## str(as(M, "graphNEL"))

pbm <- function(type = "--", rootParentEnv = NULL,
                initCells = list(), cells = list(), ...){
    initCells <- c(initCells, list(...))
    out <- new("PBM", type = type, rootParentEnv = rootParentEnv,
               initCells = initCells, cells = cells,
               expr = expression({}))
}

###_* PBM context
PBM <- getClassDef("PBM")@defaultContext
PBM$initCells("*")
PBM$initFields(mirror = protoField(.field.current_mirror),
               mirrors = protoReadOnlyField("mirrors"),
               ## mirrors_train = protoReadOnlyField("mirrors_train"),
               ## mirrors_test = protoReadOnlyField("mirrors_test"),
               folds_names = protoReadOnlyField("folds_names"))

PBM$initMethods(resetPC =
                function(){
                    map_over_model_cells(.cells, function(cell){
                        assign("._pc_done", FALSE, cell)
                    })},
                initTestMirror =
                function(test_folds, train_folds = NULL, switch = T,
                         test_mname = NULL, train_mname = NULL){
                    folds <- get("folds", .self, inherits = F)
                    if(is.numeric(test_folds))
                        test_folds <- folds[test_folds]
                    if(is.null(train_folds))
                        train_folds <- setdiff(folds, test_folds)
                    if(!all(what <- test_folds %in% folds))
                        stop("The following test folds are not declared: ", paste(test_folds[!what], sep = ", "))
                    if(!all(what <- train_folds %in% folds))
                        stop("The following train folds are not declared: ", paste(train_folds[!what], sep = ", "))
                    if(is.null(test_mname))
                        test_mname <- paste(test_folds, collapse = ":")
                    if(is.null(train_mname))
                        train_mname <- paste(train_folds, collapse = ":")
                    .create_folded_mirror(.self, train_folds, train_mname)
                    .create_folded_predict_mirror(.self, test_folds, test_mname, train_mname)
                    ## assign("mirrors_train", union(.self[["mirrors_train"]], test_mname), .self)
                    ## assign("mirrors_test", union(.self[["mirrors_test"]], test_mname), .self)
                    if(switch)
                        .self$mirror <- test_mname
                    invisible(NULL)})

troot <- root <- PBM$`*`

### OBJECTS root,
## todo: add init value argument to protoField
root$initFields(subtype = "", 
                size = protoField(.field.size),
                names = protoField(.field.names),
                var = protoField(.field.var),
                varnames = protoField(.field.varnames),
                ## should be readonly, but cannot initialize in setFields then
                parnames = character(),
                multiparnames = character(),
                varsize = protoField(.field.varsize),
                st = protoField(.field.st),
                ll = protoField(.field.ll),
                lldim = protoField(.field.lldim),
                sims =  protoField(.field.sims), 
                ## for debugging purpose only
                .pix_v = protoReadOnlyField(".pix_v"), 
                .pix_p = protoReadOnlyField(".pix_p"), 
                .pst = protoReadOnlyField(".pst"),
                .PST = protoReadOnlyField(".PST"),
                .pv = protoReadOnlyField(".pv"),
                .PV = protoReadOnlyField(".PV"))

root$evalq({
    is.initialised.M <- FALSE
    is.initialised.R <- FALSE
    size <- 1L
    ._size_init_done <- FALSE
    ._var_init_done <- FALSE
    folds <- NULL
    ## folds_fact <- NULL
    .folds_exclusive <- NULL ## T if folds can be converted to a factor
    st <- matrix(.8)
    varsize <- 1L
    ll <- matrix(-Inf)
    ._pc_done <- FALSE
    .lldim <- NULL ## if set by user, ll is more carefully initialised
    .basic_foldable_objects <- c("st", "ll")
    .basic_foldable_inxs <- c("ixs")
    do.mc_ll <- FALSE
    do.pc_ll <- FALSE
    do.mc_ll <- FALSE
    pc_ll_groups <- "sum"
    mc_ll_groups <- "sum"
})

root$initFields(do.update = TRUE,
                do.expand.q_UPDATE = FALSE,
                do.st = TRUE,
                do.timers = TRUE,
                do.mc_st = TRUE,
                do.mc_ll = FALSE,
                mc_ll_groups = protoField(.field.mc_ll_groups),
                do.first = TRUE,
                do.last = TRUE,
                do.debug = FALSE,
                do.ll = TRUE,
                do.pc_st = TRUE,
                do.pc_ll = FALSE,
                pc_ll_groups = protoField(.field.pc_ll_groups),
                do.children.ll = TRUE,  ## every cell with children must do!!
                ..ll_is_old.. = TRUE)

root$initFields(children = protoField(.field.children),
                parents = protoField(.field.parents),
                ixs = list(),
                ixs_dim = list(),
                protocol = list(),
                foldable_objects = character(),
                ## fix = integer(),
                folds = protoField(.field.folds),
                folds_factor = protoField(.field.folds_factor),
                folds_names = protoField(.field.folds_names),
                ## folds_fact = protoField(.field.folds_fact),
                pc_st = array(double()),
                mc_st = array(double()),
                pc_ll = array(double()),
                mc_ll = array(double()))
## fld_train = protoField(.field.train),
## fld_predict = protoField(.field.fld_predict))

root$initMethods(getMC =
                 function(rv_subset=NULL, ll_subset=NULL, start=1, ll=FALSE){
                     getMC(.self, rv_subset=rv_subset, ll_subset=ll_subset, start=start, ll=ll)
                 },
                 predict =
                 function(){
                     if(!isPredict(.self))
                         stop(.getType(.self), "is not a prediction cell")
                     .resetPC <- function(.cell){
                         assign("._pc_done", FALSE, .cell)
                         lapply(get("parents", .cell), .resetPC)
                     }
                     .resetPC(.self)
                     .set_pc_st_ll(.self)
                 })

## SETTERS root
root$initForms(set = form(
                   st = form(
                       stop.pbm("set.st is not specified.")),
                   ll = form(
                       stop.pbm("set.ll is not specified.")),
                   pc_st_ll = form(.set_pc_st(.self)),
                   children = form(
                       ll = form(
                           for(ch in children){
                               evalq(e(set.ll), ch)
                           }),
                       ll_is_old = form(
                           for(ch in children){
                               ch[["..ll_is_old.."]] <- TRUE
                           }),
                       ll_is_updated  = form(
                           for(ch in children){
                               ch[["..ll_is_old.."]] <- FALSE
                           })
                       )),
               setrand.st = quote(stop("set rand is not implemented for cell type ", .getType(.self))),
               set_if_old = form(
                   ll = form(
                       if(..ll_is_old..){
                           e(set.ll)
                       }),
                   children = form(
                       ll = form(
                           for(ch in children){
                               evalq(e(set_if_old.ll), ch)
                           }))))

root$initForms(init = form(
                   C = form(
                       build  = form(),
                       validate = form(
                           size = form(.validate_protocol_size(protocol$size, .self)),
                           vardim = form(.validate_protocol_vardim(protocol$vardim, .self)),
                           varnames = form(.validate_protocol_varnames(protocol$varnames, .self))
                           ## lldim = as.form(d_validate.lldim)
                           ## dimnames = as.form(d_validate.dimnames)
                           )),
                   M = form(
                       inform = form(cat(" *initializing :", .getType(.self), "\n")),
                       build = form(
                           ## assumes children and parents are already set
                           st = form(st <- st), ## copy for the convenience of foo[["st"]]
                           ll = form(ll[] <- -Inf), 
                           `_ll` = form(`_ll` <- ll), ## _ll is needed in all cells
                           `_st` = form(`_st` <- st), ## _st is not needed in all cells but placed for consistency
                           ixs = as.form(d_build.ixs),
                           child = form({child <- children[[1L]]}),
                           pos_in_C = form(
                               pos_in_C <- sapply(names(children),
                                                  function(nm) .get_pos_bc(children[[nm]][["parents"]], .self))),
                           pos_in_P = form(
                               pos_in_P <- sapply(names(parents), function(nm) .get_pos_bc(parents[[nm]][["children"]], .self))),
                           ## API: assign pst, pix and pv only if not set by user or
                           ## from cpars argument of defP
                           pix = form(for(pname in names(parents)) .assign_pix_maybe(.self, pname)),
                           pv = form({
                               for(pname in names(parents))
                                   ## assign all the actual parent variable names, if any
                                   .assign_pv_maybe(.self, parents[[pname]]$varnames, pname)
                               ## associate and assign this cells parameters
                               ._anas <- unique(c(setdiff(c(parnames, multiparnames), names(.pv)),
                                                  setdiff(c(parnames, multiparnames), names(.PV))))
                               .assign_pv_maybe(.self, ._anas) # consider all parents
                           })),
                       validate = form(
                           parents = form(.validate_protocol_parents(protocol$parents, .self)),
                           children = as.form(d_validate.children),
                           ixs = as.form(d_validate.ixs)
                           )),
                   R = form(
                       inform = form(cat(" *initializing (each_run):", .getType(.self), "\n")),
                       build = form(
                           timer = form(.t <- 0),
                           mc_ll = form({
                               if(do.mc_ll){
                                   if(length(mc_ll) == 0L){
                                       mc_ll <- .initial_xx_ll("mc", .self, 0L)
                                   }
                                   .prev_size <- nrow(mc_ll)
                                   mc_ll <- abind(mc_ll, array(NA_real_, dim = c(.nr_iter, dim(mc_ll)[ -1])), along = 1)
                                   attr(mc_ll, "prev_size") <- .prev_size
                               }}),
                           update_mc_ll = form(
                               ## depends on mc_ll_groups settings
                               update.main.mc_ll <- .gen_form_update.mc_ll(.self)),
                           mc_st = form({
                               if(do.mc_st){
                                   if(length(mc_st) == 0L){
                                       mc_st <- array(dim=c(0, dim(st)), dimnames=c(list(mc=NULL), dimnames(st)))
                                   }
                                   ._prev_size <- nrow(mc_st)
                                   mc_st <- abind(mc_st, array(NA_real_, dim=c(.nr_iter, dim(mc_st)[ -1])), along=1)
                                   attr(mc_st, "prev_size") <- ._prev_size
                               }})),
                     validate = form(
                         C = e(init.C.validate),
                         M = e(init.M.validate),
                         pixs = as.form(d_validate.pix),
                         ## pvs = as.form(d_validate.pv),
                         PVS = as.form(d_validate.PV)))),
               update = form(
                   pre = form(
                       timer = form(
                           if((!do.debug) && do.timers) .t_i <- proc.time()[[3]]
                           ),
                       debug = form(
                           if(do.debug) {
                               cat("\n * updating : \"")
                               cat(.getType(.self), "\"\n")
                               browser()
                           })
                       ),
                   main = form(
                       st = form(
                           if( do.st) e(set.st)),
                       ll = form(
                           if(do.ll){
                               e(set_if_old.ll)
                               ..ll_is_old.. <- FALSE
                           }),
                       children = form(
                           ll_is_old = form(
                               if(do.children.ll)  e(set.children.ll_is_old))),
                       mc_st = form(
                           if(do.mc_st) mc_st[.I + attr(mc_st, 'prev_size'),, ] <- st),
                       mc_ll = form(
                           ## auto-generated based on mc_ll_groups
                           if(do.mc_ll) stop.pbm("update.main.mc_ll is not defined.")),
                   post = form(
                       timer = form(
                           if((!do.debug) && do.timers) .t <- .t + proc.time()[[3]]-.t_i)
                       ))))

PBM$initCells(defBC(type = "hc", prototype="*",
                    ## do.mc_st=FALSE,
                    ## do.st=FALSE,
                    setFields = list(do.update=FALSE,
                        do.pc_st = FALSE,
                        do.pc_ll = FALSE),
                    initForms = list(
                        init.C.build.mc_ll = form(),
                        init.C.build.mc_st = form(),
                        init.C.build._ll   = form()),
                    setForms = list(
                        init.M.validate = form(),
                        init.R = form()),
                    expr = expression({
                        ## .basic_foldable_objects <- c("st", "ll")
                        .basic_foldable_inxs <- c()
                        init.M.build <- init.M.build[c("st", "child", "pos_in_C")]
                    })))

###_* DC
PBM$initCells(defBC(type = "dc", prototype = "*",
                    setFields = list(
                        do.mc_st = FALSE,
                        do.st = FALSE,
                        do.pc_st = FALSE,
                        do.pc_ll = TRUE),
                    expr = expression({
                        ## not coping ST in mirrored DCels, could be very big
                        ## all cells using child[["st"]] sould better use get("st", child),
                        ## otherwise it is safe.
                        ## not sure how to solve this more consi
                        init.M.build[c("child", "pos_in_C")] <- NULL
                    })))

###_* UC
PBM$initCells(defBC(type = "uc", prototype = "*"))


###_ + like likelihood cells
PBM$initCells(defBC(type = "like", prototype="uc",
                    initFields = list(
                        adapt = protoField(.field.adapt),
                        scale = protoField(.field.scale),
                        TR = function(x) x
                        ),
                    expr = expression({
                        scale <- 1
                        .basic_foldable_objects <- c(.basic_foldable_objects, "scale")
                    }), ## default
                    initForms = list(
### SETTERS (like)
                        set.ll_children = form({
                            for(nm in names(children)){
                                ## todo: considerable speed up could be achieved here
                                ## TODO: if ch[["ix"]] is 1:nr then not to do next step!!
                                ll_children[, nm] <-
                                    rowsum.default(c(children[[nm]][["ll"]]),
                                                   children[[nm]][["pix"]](, .type))
                            }}),
### SYNCS (like)
                        resync = form(
                            ll = form({
                                ll[rejects] <- `_ll`[rejects]
                                `_ll`[accepts] <- ll[accepts]} #_ll is set each update
                                        #_ll_all and ll_all are recalculated at each step any how
                                ),
                            st = expression()
                            ),
                        resync.children.ll = form({
                            for(nm in names(children)){
                                ## fixme: what is wrong with ix_regects <- !ix_accepts??
                                ix_accepts <- accepts[children[[nm]][["pix"]](, .type)]
                                ix_rejects <- rejects[children[[nm]][["pix"]](, .type)]
                                children[[nm]][["ll"]][ix_rejects] <- children[[nm]][["_ll"]][ix_rejects]
                                children[[nm]][["_ll"]][ix_accepts] <- children[[nm]][["ll"]][ix_accepts]
                            }}),
### INIT.M
                        init.M.build = form(
                            ll_all = form(ll_all <- array(-Inf, dim = size)),
                            rejects_accepts = form({
                                rejects <- accepts <- vector(length=size)
                                rej <- vector("integer", length=size)
                            }),
                            `resync,st` = form({
                                rep_commas <- paste(rep(",", length(dim(st))-1), collapse=" ")
                                assign("resync.st", parse(text=paste("{ st[rejects", rep_commas, " ] <- `_st`[rejects", rep_commas, " ]\n",
                                                              "`_st`[accepts", rep_commas, " ] <- st[accepts", rep_commas, "]}")))
                                remove(rep_commas)}
                                ),
                            ll_children = form(
                                ll_children <- array(-Inf, dim=c(size, length(children)), dimnames=list(NULL, names(children)))
                                ),
                            is_nrLlDim1 = form(
                                is_nrLlDim1 <- length(dim(ll))==1L || (length(dim(ll)==2L)&&dim(ll)[[2]]==1L)
                                ),
                            is_nrChildren1 = form( ## it's just length of children's vector?  waf did I do here?
                                is_nrChildren1 <- length(dim(ll_children))==1L || (length(dim(ll_children)==2L)&&dim(ll_children)[[2]]==1L)
                                ),
                            scale = form(.self$scale <- scale) ## make sure scale is ok
                            ))))


## .infoForms(PBM$like)
## .summaryForms(PBM$like)

###_  * MHrw (random walk)
PBM$initCells(defBC(type = "MHrw",
                    prototype = "like.uc",
                    setForms = list(
                        set.st = form({
                            e(set_if_old.children.ll)
                            e(set.ll_children)
                            e(set_if_old.ll)
                            `_ll` <- ll
                            `_ll_all` <- (if(is_nrLlDim1) ll else rowSums(ll)) +
                                (if(is_nrChildren1) c(ll_children) else rowSums(ll_children))
                            e(set.st_proposal) #!
                            e(set.children.ll)
                            e(set.ll_children)
                            e(set.ll)
                            ll_all <- (if(is_nrLlDim1) ll else rowSums(ll)) +
                                (if(is_nrChildren1) c(ll_children) else rowSums(ll_children))
                            e(set.alphas)
                            u <- runif(n=size, min=0, max=1)
                            rejects <- (u > alphas) | is.na(alphas)
                            accepts <- !rejects
                            rej <- rej+rejects
                            e(resync.st)}),
                        ## todo: check, why is this overwriting here?
                        update.main.ll = form(e(resync.ll))),
                    initForms = list(
                        init.M.build.pars.proposal = expression(),
                        set.st_proposal = expression(stop.pbm(" set.st_proposal is not specified.")),
                        set.alphas = expression(stop.pbm(" set.alphas is not specified.")),
                        ## ..ll_is_old.. <- FALSE (done outumaticaly in UPDATE)
                        ## may be not a good idea...if set manually like in
                        ## update.children.ll provides more uniformity.
                        update.main.children.ll = expression(
                            e(resync.children.ll),
                            e(set.children.ll_is_updated)
                            ))))

PBM$initCells(defBC(type = "norm", prototype="MHrw.like.uc",
                    setForms = list(
                      set.st_proposal = form(
                        st[] <- `_st` + c(rnorm(length(st)) * scale)),
                      set.alphas = form({
                          ## instr_lratios <- rowSums(log(st)-log(`_st`))
                          alphas <- exp(ll_all  - `_ll_all`)
                      }))))

PBM$initCells(defBC(type = "lnorm", prototype="MHrw.like.uc",
                    setFields = list(TR = log),
                    setForms = list(
                        set.st_proposal = expression(
                            st[] <- exp(TR(`_st`) + c(rnorm(length(st)) * scale))),
                        set.alphas = expression({
                            instr_lratios <- rowSums(log(st)-log(`_st`))
                            alphas <- exp(ll_all  - `_ll_all` + instr_lratios)
                        }))))

###_   + unif
PBM$initCells(defBC(type = "unif", prototype="MHrw.like.uc",
                    initFields = list(min = 0, max = 1), ## these are not distr parameters, but sampler parameters
                    initForms = list(
                        init.R.build_min_max = form({
                            min <- rep_len(min, length(st))
                            max <- rep_len(max, length(st))
                        })),
                    setForms = list(
                        set.st_proposal = form(
                            st[] <- runif(size, min = min, max = max)),
                        set.alphas = quote({
                            alphas <- exp(ll_all  - `_ll_all`)
                        }))))

###_   + unif discrete
PBM$initCells(defBC(type = "discr", prototype="unif.MHrw.like.uc",
                    setForms = list(
                        set.st_proposal = form(
                            for(i in 1:varsize){
                                #fixme: min is a vector like ix?
                                st[, i] <- sample(min[i]:max[i], size, replace = TRUE)
                            }))))

###_ + CONJ  (conjugate cells)
PBM$initCells(defBC(type = "conj", prototype = "uc",
                    initForms = list(
                        `init.M.build.set,st` = form({
                            ## fixme: this check is probably not needed
                            if(length(children) > 1) stop.pbm("multiple children are not allowed in conj cell.")
                            if(exists("set.st0", .self))
                                set.st <- get("set.st0")
                        }),
                        ## type should be set by mixins, otherwise this is useless
                        init.M.validate.child_type = form())))

## PBM$GaNorm
## infoForms(PBM$GaNorm, all.names = T)
## PBM$conj$set


source("~/works/pbm/R/distributions.R")
