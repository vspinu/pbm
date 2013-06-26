.joint_mfold <- function(pred_cell){
    stopifnot(isPredict(pred_cell))
    train_cell <- get("._train_cell", pred_cell)
    get("._mirr_fold", pred_cell, inherits = F) & 
        get("._mirr_fold", train_cell, inherits = F)
}

.new_mfold <- function(pred_cell){
    ## logical vector giving elements that are to be predicted and not already
    ## set in train fold
    stopifnot(isPredict(pred_cell))
    train_cell <- get("._train_cell", pred_cell)
    get("._mirr_fold", pred_cell, inherits = F) & 
        !get("._mirr_fold", train_cell, inherits = F)
}

.subset <- function(obj, ix, dim, drop = T){
    nrdim <- length(dim(obj))
    sub <- rep.int(alist(a =), nrdim + 1L)
    names(sub) <- NULL
    sub[[1]] <- obj
    sub[[dim + 1L]] <- ix
    sub[["drop"]] <- drop
    do.call("[", sub)
}

.subset2 <- function(obj, ix, drop = T)
    .subset(obj, ix, 2, drop = drop)

.subset1 <- function(obj, ix, drop = T)
    .subset(obj, ix, 1, drop = drop)

.set <- function(obj, ix, dim, val){
    nrdim <- length(dim(obj))
    sub <- rep.int(alist(a =), nrdim + 1L)
    names(sub) <- NULL
    sub[[dim + 1L]] <- ix
    sub[[1]] <- obj
    sub[["value"]] <- val
    do.call("[<-", sub)
}

.gen_setter1 <- function(obj, ix_str, val_str){
    name <- as.character(as.name(substitute(obj)))
    nrdim <- length(dim(obj))
    (parse(text = paste0(name, "[", ix_str, paste0(rep.int(", ", nrdim-1), collapse = ""), "] <- ", val_str)))
}
.gen_setter2 <- function(obj_str, ix_str, nrdim, val_str){
    (parse(text = paste0(obj_str, "[, ", ix_str, paste0(rep.int(", ", nrdim-2), collapse = ""), "] <- ", val_str)))
}


## .set2 <- function(obj, ix, val)
##     .set(obj, ix, 2, val)

## .set1 <- function(obj, ix, val)
##     .set(obj, ix, 1, val)


.parents_set_ith_st_from_pc <- function(pEnv, i)
    lapply(pEnv$parents, function(p){
        p <- as.environment(p)
        if(get("do.pc_st", p)){
            p$st[] <- .subset1(p$pc_st, i) # don't change the dim
        }})

## set pc_st and pc_ll, taking into acount do.pc_ll that can be logical,
## character, or factor
.set_pc_st_ll <- function(pred_cell){
    pEnv <- as.environment(pred_cell)
    tEnv <- as.environment(pEnv$._train_cell)
    do.pc_st <- get("do.pc_st", pEnv)
    do.pc_ll <- get("do.pc_ll", pEnv)
    groups <- .get_ll_groups(do.pc_ll, pEnv)
    pc_ll_value <- .get_xx_ll_asigned_value(pEnv$pc_ll_groups, "pEnv$ll")
    
    ## can be called multiple times by child cells, don't redo
    if(!get("._pc_done", pEnv) &&
       (do.pc_st || do.pc_ll)){ # for HCels this is FALSE
        ## 0) first initialize the pc_st and pc_ll arrays if needed
        .sims <- get("mirrors", get(".homeContext", tEnv))[[tEnv$._mirror_name]]@sims
        nrSims <- max(.sims$end - .sims$start + 1L, 0L)
        ## dim is NULL if objects are not defined
        if(do.pc_st){
            if(is.null(dm <- dim(pEnv$st))) dm <- 1L
            pc_st <- array(, dim = c(nrSims, dm), dimnames = c(list(NULL), dimnames(pEnv$st)))
        }
        if(do.pc_ll){
            pc_ll <- .initial_xx_ll("pc", pEnv, nrSims)
        }
        ## 0.5) init all parents
        lapply(pEnv$parents, .set_pc_st_ll)

        ## Two cases:
        ## 1) values are not known, simulate from priors
        if(any(nix <- .new_mfold(pred_cell))){
            for(i in 1:nrSims){
                .parents_set_ith_st_from_pc(pEnv, i)
                if(do.pc_st){
                    evalq(e(setrand.st), pEnv)
                    eval(.gen_setter1(pc_st, "i", "pEnv$st"))
                    ## pc_st[] <- .set1(pc_st, i, pEnv$st)
                }
                if(do.pc_ll){
                    evalq(e(set.ll), pEnv)
                    eval(.gen_setter1(pc_ll, "i", pc_ll_value))
                }
            }}
        
        ## 2) values are known from posterior simulation, just use those
        if(any(jix <- .joint_mfold(pred_cell))){
            jsub <- jix[pEnv[["._mirr_fold"]]] ## subset of joined ix in pred_cell
            eval(.gen_setter2("pc_st", "jsub", length(dim(pc_st)), ".subset2(tEnv$mc_st, jix)"))
            ## pc_st <- .set2(pc_st, jsub, .subset2(tEnv$mc_st, jix))
            if(do.pc_ll){
                for(i in 1:nrSims){
                    pEnv$st[] <- .subset1(pc_st, i)
                    .parents_set_ith_st_from_pc(pEnv, i)
                    evalq(e(set.ll))
                    eval(.gen_setter1(pc_ll, "i", pc_ll_value))
                }}}
        if(do.pc_st)
            pEnv$pc_st <- pc_st
        if(do.pc_ll)
            pEnv$pc_ll <- pc_ll
    }
    pEnv$._pc_done <- TRUE
}

.set_folds_factor <- function(pEnv, error = T){
    pEnv <- as.environment(pEnv)
    local_folds <-
        if(isMirror(pEnv))
            get("folds", pEnv)[pEnv$._mirr_fold, pEnv$folds_names, drop = F]
        else
            get("folds", pEnv)
    if(pEnv$.folds_exclusive <- all(rowSums(local_folds) == 1L)){
        pEnv$folds_factor <- .ind2factor(local_folds)
    }else{
        mess <- "do.pc_ll is 'fold', but local folds are non-exclusive."
        if(error) stop.pbm(mess) else warn.pbm(mess)
    }}

.initial_xx_ll <- function(xx = "mc", object, nr){
    groups <- get(paste0(xx, "_ll_groups"), object)
    if(identical(groups, "sum")){
        dnames <- "LL"
        dm <- 1L
    }else if(identical(groups, "all") || identical(groups, TRUE)){
        ll <- get("ll", object)
        dnames <- dimnames(ll)
        if(is.null(dm <- dim(ll))) dm <- length(ll)
    } else {
        dnames <- levels(.get_ll_groups(groups, object)) # rowsum below preserves this
        dm <- length(dnames)
    }
    array(NA_real_, dim = c(nr, dm), dimnames = c(list(NULL), list(dnames)))
}

.get_ll_groups <- function(xx_ll_groups, object){
    objEnv <- as.environment(object)
    gr <- 
        if((identical(xx_ll_groups, "folds"))){
            if(is.null(objEnv$.folds_exclusive)){
                ## if non null, it have been already set before;
                ## also set .folds_exclusive
                .set_folds_factor(objEnv, error = TRUE)
            }
            ## .fold_exclusive is set by .set_folds_factor
                objEnv$folds_factor
        }else if (is.character(xx_ll_groups) && xx_ll_groups %in% names(objEnv$parents)){
            objEnv$ixs[[xx_ll_groups]]
        }else if(is.factor(xx_ll_groups) || is.integer(xx_ll_groups)){
            xx_ll_groups
        }else NULL
    as.factor(rep_len(gr, get("size", objEnv)))
}

.gen_form_update.mc_ll <- function(.self){
    with(.self, {
        if(do.mc_ll){
            value <-.get_xx_ll_asigned_value(mc_ll_groups)
            as.form(.gen_setter1(mc_ll, "attr(mc_ll, 'prev_size') + .I", value))
        }else form()
    })}

.get_xx_ll_asigned_value <- function(xx_ll_groups, ll_name = "ll"){
    name <- as.character(substitute(xx_ll_groups))
    if(identical(xx_ll_groups, "sum"))
        sprintf("sum(%s)", ll_name)
    else if(identical(xx_ll_groups, "all") || identical(xx_ll_groups, TRUE))
        ll_name
    else if(identical(xx_ll_groups, "rowSums"))
        sprintf("rowSums(%s)", ll_name)
    else sprintf("rowsum(c(%s), .get_ll_groups(%s, .self))", ll_name, name) # todo: store groups on first call
}
