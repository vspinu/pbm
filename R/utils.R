### testthat enhancement
expect_close <- function (object, expected, tol = .005, ..., info = NULL,
                          label = NULL, expected.label = NULL)
{
    if (is.null(label)) {
        label <- testthat:::find_expr("object")
    }
    if (is.null(expected.label)) {
        expected.label <- testthat:::find_expr("expected")
    }
    expect_that(object,
                is_close_to(expected, tol = tol, label = expected.label),
                info = info, label = label)
}

is_close_to <- function (object, label = NULL, tol = .005, ...)
{
    if (is.null(label)) {
        label <- testthat:::find_expr("object")
    }else if (!is.character(label) || length(label) != 1) {
        label <- deparse(label)
    }
    function(actual) {
        same <- all.equal(actual, object, tolerance = tol, scale = 1, ...)
        expectation(identical(same, TRUE),
                    paste0("is not close to ", label, "\n", same))
    }
}



### .pv, .pst, .pixs_v, .pix_p etc construction

.assign_pst_maybe <- function(.self, pname, force = FALSE){
    env <- as.environment(.self)
    if(force || is.null(env$.pst[[pname]])){
        env$.pst[[pname]] <-
            bquote(parents[[.(pname)]][["st"]])
    }
    if(force || is.null(env$.PST[[pname]]))
        env$.PST[[pname]] <- 
            bquote(e(.pst[[.(pname)]])[e(.pix_p[[.(pname)]]), ])
}

.assign_pix_maybe <- function(.self, pname, force = FALSE){
    env <- as.environment(.self)
    if(force || is.null(env$.pix_p[[pname]]))
        env$.pix_p[[pname]] <- bquote(ixs[[.(pname)]])
}

.assign_pv_paired_maybe <- function(.self, pars, pname, force = FALSE, check = TRUE){
    ## pars is a named character c(cpar = ppar) where cpar is a name of the
    ## parameter in child and ppar is the name in parent.
    parent <- get("parents", .self, inherits = F)[[pname]]
    varnames <- parent$varnames
    env <- as.environment(.self)
    
    mapply(function(cpar, ppar){
        if(check && (!cpar %in% varnames))
            stop(sprintf("'%s' declared in cpars is not a valid variable name in parent object '%s'",
                         cpar, pname))

        ## this one should be here, because pname is dynamically computed in
        ## .assign_pv_maybe. fixme: this will likely cause trouble when cix is a
        ## form and this code is not executed :(
        if(force || is.null(env$.parent[[cpar]]))
            env$.parent[[cpar]] <- bquote(parents[[.(pname)]])
        
        if(force || is.null(env$.pix_v[[cpar]]))
            env$.pix_v[[cpar]] <- bquote(e(.pix_p[[.(pname)]]))
        
        if(force || is.null(env$.pv[[cpar]]))
            env$.pv[[cpar]] <- 
                bquote(e(.pst_p[[.(pname)]])[, .(ppar), drop = F])

        if(force || is.null(env$.PV[[cpar]]))
            env$.PV[[cpar]] <- 
                bquote(e(.pv[[.(cpar)]])[e(.pix_v[[.(cpar)]]),, drop = F])
        
    }, names(pars), c(pars))
}

.assign_pv_maybe <- function(.self, pars, pname = NULL, force = FALSE){
    ## Assign .pv[["xxx]],  .PV[["xxx]] and .pixv[["xxx]]. For each 'xxx' in pars:
    ##
    ## If pv[["xxx"]] is already installed in .self, do nothing. Try hard to
    ## associate varnames from parents to pv.mame forms. If xxx is found in
    ## parent's varnames, use it, otherwise, if no varnames found in parent,
    ## create those by side effect and give the first parent's variable name
    ## 'xxx'. Otherwise check for the first empty parent's varname "" and give
    ## it the 'xxx'. If there are not empty names in parents's varnames then
    ## just reparametrize (that is pv[["xxx"]] points to st[, reparam_ix]) with a
    ## warning.
    ##
    ## If parent is a multivariate, and 'xxx' is not found in parents varnames
    ## then pv[["xxx"]] will return a whole st.
    ##
    ## This function is called on parent per parent in new("BC") on cpar of
    ## parents. And then on all the uninitialised parameters and all the parents
    ## (pname = NULL) in init.R.build.

    parents <- get("parents", .self)
    if(!is.null(pname))
        parents <- parents[pname]

    ## vector all varnames as names, and coresponding parents' names as values
    var_p <-
        mapply(function(p, nm){
            nms <- p$varnames # always nonnull character vector
            structure(rep.int(nm, length(nms)), names = nms)
        }, parents, names(parents), SIMPLIFY = F)
    names(var_p) <- NULL
    var_p <- unlist(var_p)

    ## Only makes sense when all parents are involved
    if(is.null(pname)){
        if(length(pars) > (allsize <- sum(sapply(parents, function(x) x$varsize)))){
            stop.pbm(sprintf("Building .pv forms failed; there are more parameters (%s) than all parents can accept (%s)
Did you specify cpars in defP?", length(pars), allsize))
        }}

    env <- as.environment(.self)
    warned <- FALSE
    reparam_ix <- 1L
    opname <- pname
    for(par in pars[nzchar(pars)]){
        forms <- list()
        pname <- opname
        if(force || is.null(env$.pv[[par]])){
            if(par %in% .self$multiparnames){
                ## if multivariate, just return the whole parent st
                if(is.null(pname)){
                    pname <- var_p[[reparam_ix]]
                    reparam_ix <- reparam_ix + parents[[pname]]$varsize
                    warn.pbm(sprintf("multivariate parameter not found in parents' varnames,
   matching '%s:%s' ---> '%s:st'  ", .getType(.self, fullName = F), par, pname))
                }
                env$.pv[[par]] <- bquote(e(.pst[[.(pname)]]))
            }else{
                ## if univariate ...
                if(!(par %in% names(var_p))){
                    empty <- which(names(var_p) == "")
                    if(length(empty)){
                        ## ... assign by side effect varnames in parent st.
                        pname <- var_p[[empty[[1L]]]]
                        which_var <- which(pname == var_p)[[1]]
                        thispvars <- names(var_p)[which_var]
                        thispvars[[which(thispvars == "")[[1]]]] <- par
                        parents[[pname]]$varnames <- thispvars
                        if(warned)
                            cat(sprintf("   setting '%s:%s' ---> '%s:%s'\n",
                                        .getType(.self, fullName = F), par, pname, par))
                        else warn.pbm(sprintf("empty parent varnames detected,
   setting '%s:%s' ---> '%s:%s'  ", .getType(.self, fullName = F), par, pname, par))
                        warned <- TRUE
                        var_p <- var_p[-which_var]
                    }else{
                        ## ... associate
                        if(reparam_ix > length(var_p))
                            next ## fixme: happens when child has more parameters than this parent
                        pname <- var_p[[reparam_ix]]
                        new_par <- names(var_p)[[reparam_ix]]
                        warn.pbm(sprintf("parameter '%s' not found in parents' varnames,
   matching '%s:%s' ---> '%s:%s'  ", par, .getType(.self, fullName = F), par, pname, new_par))
                        par <- new_par
                        reparam_ix <- reparam_ix + 1L
                    }
                }else pname <- var_p[[which(par == names(var_p))[[1]]]]
                env$.pv[[par]] <- bquote(e(.pst[[.(pname)]])[, .(par), drop = F])
            }
        }

        if(!is.null(pname)){
            pair <- structure(par, names = par)
            .assign_pv_paired_maybe(.self, pair, pname, force = force, FALSE)
        }
    }
}



### ERROR HANDLING

stop.pbm <- function(..., env = parent.frame()){
    ## Stop with report of bcell in which the error occured
    while(!(found_self <- exists(".self", env))&&!identical(env, .GlobalEnv))
        env <- parent.env(env)
    if(found_self)
        bcn <- paste("in bcell '", .getType(get(".self", env), full = F), "' : \n", sep = "")
    else bcn <- "\n"
    stop("pbm -> ", bcn, ...,"\n", call. = FALSE)
}

warn.pbm <- function(..., env = parent.frame()){
    ## Warn with report of bcell in which the error occured
    while(!(found_self <- exists(".self", env)) &&! identical(env, .GlobalEnv))
        env <- parent.env(env)
    if(found_self)
        bcn <- paste("bcell '", .getType(get(".self", env), fullName = T), "' ", sep = "")
    else bcn <- "??? "
    invisible(cat("pbm: in ", bcn,  ..., fill=T, sep=""))
}



### MISC

.get_model_cells <- function(.cells){
    m_names <- leafNames(.cells)
    model_cells <- lapply(m_names, get, envir = .cells)
    names(model_cells) <- m_names
    model_cells
}

.get_pos_bc <- function(l, bc){
    not.identical <- TRUE
    i <- 0
    while(not.identical){
        i <- i+1
        if(i > length(l)) stop("BC was not found in the list of BCs.")
        not.identical <- !identical(as.environment(l[[i]]), as.environment(bc))
    }
    return(i)
}

getMC <- function(object, rv_subset = NULL, ll_subset = NULL, start=1,
                  ll=FALSE, transform_coda = TRUE){
    object <- as.environment(object)
    if(ll){
        stopifnot(start < NROW(object[["mc_ll"]]))
        ix <- start:NROW(object[["mc_ll"]])
        return(object[["mc_ll"]][ix])
    }
    stopifnot(start < NROW(object[["mc_st"]]))
    ix <- start:NROW(object$mc_st)
    st_subset <- rep(alist(a=), length(dim(object[["st"]])))
    if(!is.null(rv_subset)){
        if(is.atomic(rv_subset)) rv_subset <- list(rv_subset)
        if(length(rv_subset)!= length(object[["rv_dim"]]))
            stop.pbm(gettextf("length or 'rv_subset' (%s) is not the same as the nr of dimension of rv (%s)",
                              length(rv_subset), length(object$rv_dim)))
        st_subset[(length(dim(object[["st"]])) - length(object$rv_dim) + 1):length(dim(object[["st"]]))] <- rv_subset
    }
    if(!is.null(ll_subset)){
        if(is.atomic(ll_subset)) ll_subset <- list(ll_subset)
        if(length(ll_subset)!= length(dim(object[["ll"]])))
            stop.pbm(gettextf("length or 'll_subset' (%s) is not the same as the nr of dimension of ll (%s)",
                              length(ll_subset), length(object$ll_dim)))
        st_subset[1:length(object$ll_dim)] <- ll_subset
    }
    mc <- do.call("[", c(quote(object[["mc_st"]]), list(ix), st_subset, list(drop=FALSE)))
    if(transform_coda){
        dmn <- dimnames(mc)[-1L]
        null_dmn <- sapply(dmn, is.null)
        null_dm_len <- dim(mc)[-1L][null_dmn]
        new_names <- lapply(seq_len(sum(null_dmn)), function(i) if(null_dm_len[[i]] > 1L) seq_len(null_dm_len[[i]]) else "") #make.names(rep.int("", null_dm_len[[i]]), unique=TRUE))
        dmn[null_dmn] <- new_names
        join_names <- function(str1, str2){
            unlist(lapply(str2, function(s2) paste(str1, s2, sep="")))
        }
        flat_names <- Reduce(join_names, dmn)
        dim(mc) <- c(dim(mc)[[1L]], prod(dim(mc)[-1L]))
        dimnames(mc) <- c(list(NULL), list(flat_names))
        library(coda)
        mcmc(mc, start=start)
    }else{
        mc
    }
}

.strip_mcs <- function(cells){
    ## leaves only the last state of the chain in the mc objects
    ## fixme: todo: place this in init.M or even !!!
    lapply(cells, function(e){
        evalq(envir=e, expr={
            delayedAssign("mc_st", array(dim=c(0, dim(st)), dimnames=c(list(mc=NULL), dimnames(st))))
            delayedAssign("mc_ll", array(dim=c(0, length(ll)), dimnames=c(list(mc=NULL), list(names(ll)))))
        })
    })
    invisible(NULL)
}
