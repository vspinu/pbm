### ROOT BCELL
PBM$initCells("*")
root <- PBM$`*`

## todo: add init value argument to protoField

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
                ll_tr = matrix(0), 
                lldim = protoField(.field.lldim),
                sims =  protoField(.field.sims),
                ## for debugging purpose only
                do.st = TRUE,
                do.timers = TRUE,
                do.mc_st = TRUE,
                do.mc_ll = FALSE,
                mc_ll_groups = protoField(.field.mc_ll_groups),
                do.first = TRUE,
                do.last = TRUE,
                do.debug = FALSE,
                do.ll = TRUE,
                do.ll_tr = TRUE,
                do.pc_st = TRUE,
                do.pc_ll = FALSE,
                pc_ll_groups = protoField(.field.pc_ll_groups),
                do.update = TRUE, 
                ..ll_is_old.. = TRUE,
                ..st_is_old.. = TRUE,
                ..ctoken..= -1, 
                pc_st = array(double()), 
                mc_st = array(double()), 
                pc_ll = array(double()), 
                mc_ll = array(double()),
                children = protoField(.field.children),
                parents = protoField(.field.parents),
                ixs = list(),
                ixs_dim = list(), ## fixme: readonly?
                protocol = list(),
                foldable_objects = character(),
                ## fix = integer(),
                folds = protoField(.field.folds),
                folds_factor = protoField(.field.folds_factor),
                folds_names = protoField(.field.folds_names),
                ## folds_fact = protoField(.field.folds_fact),
                pc_st = .make_field_mc("pc_st"),
                mc_st = .make_field_mc("mc_st"),
                pc_ll = .make_field_mc("pc_ll"), 
                mc_ll = .make_field_mc("mc_ll"))

root$initFields(.pix_v = protoReadOnlyField(".pix_v"), 
                .pix_p = protoReadOnlyField(".pix_p"), 
                .pst = protoReadOnlyField(".pst"),
                .PST = protoReadOnlyField(".PST"),
                .pv = protoReadOnlyField(".pv"),
                .PV = protoReadOnlyField(".PV"))

root$evalq({
    ..M.init.done.. <- FALSE
    ..R.init.done.. <- FALSE
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

## ROOT INITS
root$initForms(init.C = form(
                   build  = form(),
                   validate = form(
                       size = form(.validate_protocol_size(protocol$size, .self)),
                       vardim = form(.validate_protocol_vardim(protocol$vardim, .self)),
                       varnames = form(.validate_protocol_varnames(protocol$varnames, .self))
                       ## lldim = as.form(d_validate.lldim)
                       ## dimnames = as.form(d_validate.dimnames)
                       )))

root$initForms(init.M = form(
                   inform = form(cat(" *initializing :", .getType(.self), "\n")),
                   build = form(
                       ## assumes children and parents are already set
                       tr_children = form(
                           tr_children <- children[unlist(sapply(children, protoIs, "tr.*"))]),
                       st = form(.fields$st <- st), ## set up ll, size, varnames etc (needed only when st isnot explisitely supplied)
                       `_ll` = form(`_ll` <- ll), ## _ll is needed in all cells
                       `_st` = form(`_st` <- st), ## _st is not needed in all cells but placed for consistency
                       ixs = as.form(d_build.ixs),
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
                       }),
                       pos_in_P = form(
                           pos_in_P <- sapply(names(parents),
                                              function(nm) .get_pos_bc(parents[[nm]][["children"]], .self))),
                       pos_in_C = form(
                           pos_in_C <- sapply(names(children),
                                              function(nm) .get_pos_bc(children[[nm]][["parents"]], .self)))),
                   validate = form(
                       parents = form(.validate_protocol_parents(protocol$parents, .self)),
                       children = as.form(d_validate.children),
                       ixs = as.form(d_validate.ixs)
                       )))

root$initForms(init.R = form(
                   inform = form(cat(" *initializing (each_run):", .getType(.self), "\n")),
                   build = form(
                       timer = form(.t <- 0),
                       mc_ll = form({
                           if(do.mc_ll){
                               if(length(mc_ll) == 0L){
                                   ## takes into account mc_ll_groups 
                                   mc_ll <- .initial_xx_ll("mc", .self, 0L)
                               }
                               .prev_size <- nrow(mc_ll)
                               mc_ll <- abind(mc_ll, along = 1, 
                                              array(NA_real_,
                                                    dim = c(.nr_iter, dim(mc_ll)[ -1])))
                               attr(mc_ll, "prev_size") <- .prev_size
                           }}),
                       update_mc_ll = form(
                           ## takes into account mc_ll_groups 
                           UPDATE.main.mc_ll <- .gen_form_update.mc_ll(.self)),
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
                       PVS = as.form(d_validate.PV))))


## ROOT DEV INTERFACE
## DEV API: All communication from outside should be done with internal methods
## or at least through side effect internal functions like .mark_ll_is_old below

## INTERNAL dev interface
root$initForms(set = form(
                   ll = form(
                       stop.pbm("set.ll form is not defined"),
                       doc = "Set the ll object to the loglikelihood of the current state (st object)."), 
                   st = form(
                       stop.pbm("set.st form is not defined"),
                       doc = "Set the current state of the BCell during the update. This form is set by MC algorithms. "), 
                   rand.st = form(
                       stop.pbm("set.rand.st is not defined"),
                       doc = "Set st to randomly generated values given the parent cells parameters. Not to be confused with set.st."), 
                   pc_st_ll = form(
                       .set_pc_st(.self),
                       doc = "Set pc_st and pc_ll objects, stand for Prediction Chain."),
                   ll_is_old = form(..ll_is_old.. <- TRUE),
                   st_is_old = form(..st_is_old.. <- TRUE),
                   ll_is_updated = form(..ll_is_old.. <- FALSE), 
                   st_is_updated = form(..st_is_old.. <- FALSE)))

## EXTERNAL dev interface
.resync_object <- function(.self, rejects, obj, old_obj, ptype, token, force = FALSE){
    env <- as.environment(.self)
    if(force || !identical(env[["..ctoken.."]], token)){
        ix_rejects <- rejects[env[["pix"]](, ptype)]
        ix_accepts <- !ix_rejects
        env[[obj]][ix_rejects] <- env[[old_obj]][ix_rejects]
        env[[old_obj]][ix_accepts] <- env[[obj]][ix_accepts]
        assign("..ctoken..", token, envir = env)
        TRUE
    }else FALSE
}

## API: With some exceptions updating/settin/modiffying forms don't touch other
## states.

## API: external access api (LL only, no need to manipulate ST)
root$initMethods(.resync.ll =
                 function(rejects, ptype, token){
                     .resync_object(.self, rejects, "ll", "_ll", ptype, token)
                 },
                 .update.ll = function() e(update.ll),
                 .set.ll_is_old = function() e(set.ll_is_old),
                 .set.ll_is_updated = function() e(set.ll_is_updated),
                 .get.ll =
                 function(ptype, token, ...){
                     ## for some cells like TR should return a matrix
                     ## todo: considerable speed up could be achieved here
                     ## TODO: if ch[["ix"]] is 1:nr then not to do next step!!
                     if(!identical(..ctoken.., token)){
                         assign("..ctoken..", token, envir = .self)
                         rowsum.default(c(ll), pix(, ptype))
                     }else 0
                 })

root$evalq({
    ## EXTERNAL USEFUL SIDE EFFECTS
    ## fixme:? same names are a bit dificult, for example resync.ll is form,
    ## protoMethod and function. Confusing?
    .set.ll_is_old <- function(cells){
        for(c in cells) c$.set.ll_is_old()
    }
    .set.ll_is_updated <- function(cells){
        for(c in cells) c$.set.ll_is_updated()
    }
    .update.ll <- function(cells){
        for(c in cells) c$.update.ll()
    }
    .resync.ll <- function(children, rejects, ptype, token = runif(1L))({
        ## fixme: with transform cells might happen that child is recynced
        ## multiple times
        for(ch in children)
            ch$.resync.ll(rejects, ptype, token = token)
    })
    .get.ll <- function(cells, ptype, token = runif(1L), ...){
        ## Return a matrix of LLs of the forward blanket. All ll .get_ll should
        ## have the same length or dim[[1]] (equal to psize)
        out <- lapply(cells, function(c) c$.get.ll(ptype, token = token, ...))
        do.call(cbind, out)
    }
})

## utils
root$evalq({
    is_env_member <- function(env, env_set){
        for(el in env_set)
            if(identical(env, el)) return(TRUE)
        FALSE
    }
    fastRS <- function(x){
        dm <- dim(x)
        if(is.null(dm) || dm[[2]] == 1L) x else rowSums(x)
    }
})


## ROOT UPDATE
root$initForms(update = form(
                   doc = "Set objects if old and mark all relevant state variables cell and its children. This forms should be used in programs.",
                   ## API: update.st form sets children's ll_is_old *AFTER*
                   ## executing set.st. Algorithms can adjust children's ll and
                   ## should marke them with .set.ll_is_updated method in
                   ## update.ll or UPDATE.main.ll.
                   st = form({
                       if(..st_is_old..){
                           e(set.st)
                           e(set.st_is_updated)
                           e(set.ll_is_old);	.set.ll_is_old(children)
                       }}),
                   ll = form({
                       if(..ll_is_old..){
                           e(set.ll)
                           e(set.ll_is_updated)
                       }})), 
               UPDATE = form(
                   pre = form(
                       timer = form(
                           if((!do.debug) && do.timers) .t_i <- proc.time()[[3]]
                           ),
                       debug = form(
                           if(do.debug) {
                               cat("\n * updating : \"")
                               cat(.getType(.self), "\"\n")
                               browser()
                           })),
                   main = form(
                       doc =
                       "Main routine executed during the update.
By default consists of st, ll, mc_st and mc_ll forms.", 
                       st = form(
                           if(do.st){
                               e(update.st)
                           }),
                       ll = form(
                           if(do.ll){
                               e(update.ll)
                           }),
                       mc_st = form(
                           if(do.mc_st) mc_st[.I + attr(mc_st, 'prev_size'),, ] <- st),
                       mc_ll = form(
                           ## auto-generated based on mc_ll_groups
                           if(do.mc_ll) stop.pbm("update.main.mc_ll is not defined.")),
                       tr_children = form(
                           doc = "Update all transform cells",
                           for(ch in tr_children){
                               evalq(e(UPDATE), ch) ## recursively call on transformers
                           })),
                   post = form(
                       timer = form(
                           if((!do.debug) && do.timers) .t <- .t + proc.time()[[3]]-.t_i)
                       )))


