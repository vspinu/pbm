map_over_model_cells <- function(.cells, func, ...){
    ## .cells <- as.environment(get(".cells", envir = pbm))
    root_env <- as.environment(get("*", .cells))
    m_names <- leafNames(.cells)
    model_cells <- lapply(m_names, get, envir = .cells)
    out <- lapply(model_cells, func, ...)
    names(out) <- m_names
    invisible(out)
}

map_from_current_cell <- function(cell, func, ..., .stamp = Sys.time()){
    ## map over children/parents till end, for side effects
    ## don't use it recursively except with explisite .stamp 
    cell <- as.environment(cell)
    if(is.null(cell$._map_stamp) || cell$._map_stamp != .stamp){
        func(cell$.self)
        cell$._map_stamp <- .stamp
        lapply(cell$parents, function(cl)
            map_from_current_cell(cl, func, ..., .stamp = .stamp))
        lapply(cell$children, function(cl)
               map_from_current_cell(cl, func, ..., .stamp = .stamp))
        TRUE
    }else FALSE
}

map_from_current_cell_with_old <- function(cell, old_cell = NULL, func, ..., .stamp = Sys.time()){
    ## map over children/parents till end, for side effects
    ## don't use it recursively except with explisite .stamp 
    cell <- as.environment(cell)
    if(is.null(cell$._map_stamp) || cell$._map_stamp != .stamp){
        func(cell$.self, old_cell)
        cell$._map_stamp <- .stamp
        lapply(cell$parents, function(cl)
               map_from_current_cell(cl, cell, func, ..., .stamp = .stamp))
        lapply(cell$children, function(cl)
               map_from_current_cell(cl, cell, func, ..., .stamp = .stamp))
        TRUE
    }else FALSE
}

## map_from_current_cell(M$LUR, function(cell) cat(.getType(cell), "\n" ))

isPredict <- function(obj)
    exists("._train_cell", obj, inherits = F) && !is.null(get("._train_cell", obj))

.gen_mirr_name <- function(mname)
    paste0("_", mname)

.create_mirror <- function(pbm, mname, init = T){
    pbm <- as.environment(pbm)
    if(mname == "base")
        stop( "cannot create a mirror with name 'base'")
    mirrors <- get("mirrors", pbm, inherits = F)
    .cells <- mirrors[["base"]]
    new_mirror_cont <- new("bCellContainer", is_mirror = TRUE)
    parent.env(new_mirror_cont) <- .cells
    full_mname <- .gen_mirr_name(mname)
    map_over_model_cells(.cells, function(cell){
        mirror <- new("BC", prototype = cell,
                      homeContext = pbm$.self, 
                      type = full_mname)
        assign("._mirror", T, mirror)
        assign("._mirror_name", mname, mirror)
        assign(.getType(cell), mirror, new_mirror_cont)
        assign(full_mname, mirror, cell)
    })

    map_over_model_cells(.cells, function(cell){
        parents <- get("parents", cell)
        children <- get("children", cell)
        mparents <- new("bcParents")
        mchildren <- new("bcChildren")
        for(nm in names(parents))
            mparents[[nm]] <- get(full_mname, parents[[nm]])
        for(nm in names(children))
            mchildren[[nm]] <- get(full_mname, children[[nm]])
        mirror <- get(full_mname, cell)
        assign("parents", mparents, mirror)
        assign("children", mchildren, mirror)
    })

    ## map_over_model_cells(.cells, function(cell){
    ##     mirror <- get(full_mname, cell)
    ##     .initialize_M(mirror)
    ## })
    pbm$mirrors[[mname]] <- new_mirror_cont
    ## assign("mirrors", mirrors, pbm)
    pbm$sims[[mname]] <- list(start = 1L, end = 0L)
    ## assign("sims", sims, pbm)
    if(get("current_mirror", pbm, inherits = F) == mname)
        assign(".cells", new_mirror_cont, pbm)
    if(init)
        .initialize_M(pbm$mirrors[[mname]])
}

.subfold <- function(object, mfold, oldsize){
    if(is.list(object))
        sapply(object, .subfold, mfold, oldsize, simplify = F)
    else if(is.array(object) && length(dim(object)) > 1) ## so fucked up :(
        ## adjust only  arrays that explicitely match the old size
        if(dim(object)[[1]] == oldsize) object[mfold,, drop = F]
        else object
    else if(length(object) %% oldsize == 0L) ## so fucked up :(
        object[rep(mfold, length.out = length(object))] ## this silently accepts shorter objects
    else object
}

.subfold_ixs <- function(object, parent, size){
    if(is.array(object) && length(dim(object)) > 1){
        ## adjust only arrays that explisitely match the old size
        if(dim(object)[[1]] == size)
            object[] <- get("._subfold_ix", parent)[object]
    }else object[] <- get("._subfold_ix", parent)[object]
    object
}

.create_folded_mirror <- function(pbm, FLDS, mname, init = T){
    pbm <- as.environment(pbm)
    if(is.numeric(FLDS))
        FLDS <- pbm$folds_names[FLDS] # recognize by name
    .create_mirror(pbm, mname, init = F)
    ## fixme: cleanup in case of the error?
    map_over_model_cells(pbm$mirrors[[mname]], function(cell){
        cell$initFields(folds = protoField(.field.mfolds),
                        folds_names = protoField(.field.mfolds_names))
        cell <- as.environment(cell)
        if(is.null(folds <- get("folds", cell)))
            stop("FOLDS is not defined for cell ", .getType(cell$.self))
        ## mirror's main fold consisting of all elementary folds given by FLDS
        mfold <- as.logical(rowSums(folds[, FLDS, drop = F]))
        oldsize <- get("size", cell)
        for(nm in c(get("foldable_objects", cell),
                    get(".basic_foldable_inxs", cell), 
                    get(".basic_foldable_objects", cell))){
            ## .initialize_M has been already performed
            obj <- .subfold(get(nm, cell), mfold, oldsize)
            assign(nm, obj, cell)}
        newsize <- sum(mfold)
        cell$size <- newsize
        subfix <- integer(oldsize)
        subfix[mfold] <- 1:newsize
        ## vector of 0s with 1:newsize inserted in places given by mfold
        cell$._subfold_ix <- subfix
        ## a boolean vector, giving the positions of this fold in the original vector
        cell$._mirr_fold <- mfold
        cell$folds_names <- FLDS
    })
    map_over_model_cells(pbm$mirrors[[mname]], function(cell){
        cell <- as.environment(cell)
        for(nm in get(".basic_foldable_inxs", cell)){
            obj <- get(nm, cell)
            if(is.list(obj)){
                obj <- mapply(.subfold_ixs, obj, get("parents", cell), cell$size, SIMPLIFY = F)
            }else{
                obj <- .subfold_ixs(obj, get("parent", cell), cell$size)
            }
            assign(nm, obj, cell)
        }})
    if(init){
        cat("Post index subfolding initialization:\n")
        invisible(.initialize_M(pbm$mirrors[[mname]]))
    }
}

.create_folded_predict_mirror <- function(pbm, FLDS, mname, train_mname){
    .create_folded_mirror(pbm, FLDS, mname, init = F)
    pbm <- as.environment(pbm)
    train_cont <- pbm$mirrors[[train_mname]]
    train_names <- ls(train_cont, all.names = T)
    pbm$mirrors[[mname]]@train_mirror_name <- train_mname
    pred_cont <- pbm$mirrors[[mname]]
    map_over_model_cells(pred_cont, function(cell){
        name <- .getType(cell, base_only = T)
        tr_name <- train_names[grepl(name, train_names, fixed = T)]
        tr_cell <- get(tr_name, train_cont, inherits = F)
        ## a bit of a cludge
        cell$evalq({init.M.build.mc_st = expression()
                   init.M.build.mc_ll = expression()})
        ## remove(mc_st, mc_ll, envir = cell)
        assign("._train_cell", tr_cell, cell)
        assign("._setter_redirect_to_", tr_cell, cell)
        ## kludge?
        ## replanting only now; FOLDS should be taken from base cells.
        ## not clear how to cleanly overcome this
        setPrototype(cell, tr_cell)
        assign("._train_mirror_name", train_mname, cell)
    })
    cat("Post index subfolding initialization:\n")
    invisible(.initialize_M(pbm$mirrors[[mname]]))
}
        

.set_folds <- function(cell, folds_fact){
    folds_fact <- as.factor(folds_fact)
    ## fails if in mirror
    ## re-initialize.M current cells (.cells) if not initialized (we need ixs_dim here)
    pbm <- get(".homeContext", cell)
    assign("folds_names", NULL, envir = pbm)
    assign("folds", NULL, envir = pbm)
    if(isMirror(cell))
        stop("cannot set folds in mirrored objects; reset the mirror to 'base' first")
    ## optionally initialize
    .initialize_M(.get_model_cells(get(".cells", pbm)), force = FALSE)
    cell <- as.environment(cell)
    if(length(folds_fact) != cell$size)
        stop("folds lenth (", length(folds_fact), ") is not the same as cells SIZE (", cell$size, ")")
    folds <- .class.ind(folds_fact)
    storage.mode(folds) <- "logical"
    cell$folds <- folds
    cell$folds_names <- levels(folds_fact)
    ## these two are valid only in some cells
    cell$.folds_exclusive <- TRUE
    cell$folds_factor <- folds_fact
    stamp <- Sys.time() # don't remove, this is needed in d_rec_set_folds
    eval(d_rec_set_folds)
    assign("folds", colnames(folds), envir = pbm)
}

d_rec_set_folds <- expression({
    cell$._map_stamp <- stamp
    for(i in seq_along(cell$pos_in_C))
        .set_folds_from_parent(cell$pos_in_C[[i]], cell$children[[i]], stamp)
    for(i in seq_along(cell$pos_in_P))
        .set_folds_from_child(cell$pos_in_P[[i]], i, cell$parents[[i]], stamp)
})

.set_folds_from_child <- function(cnr, pnr, cell, stamp){
    cell <- as.environment(cell)
    if(is.null(cell$._map_stamp) || cell$._map_stamp != stamp){
        child <- as.environment(cell$children[[cnr]])
        if(is.null(child))
            stop.pbm("folds are not set in '", cnr, "'")
        new_folds <- apply(child$folds, 2, function(fld){
            ## fld <- rep.int(fld, prod(c(1L, dim(ll)[-1]))) # make it size of inx0
            stop("fixme: ixs_dim are no longer there, ixs have meaning of ixs0")
            1:get("size", cell)%in% unique(child$ixs0[[pnr]][fld]) # recicle fld as needed
        })
        ## uff, apply is so stupid
        if(is.null(dim(new_folds)))
            dim(new_folds) <- c(1, length(new_folds))
        colnames(new_folds) <- colnames(child$folds)
        cell$folds <- new_folds
        ## reset 
        cell$.folds_exclusive <- NULL
        cell$folds_factor <- NULL
        eval(d_rec_set_folds)
    }}

.set_folds_from_parent <- function(pnr, cell, stamp){
    cell <- as.environment(cell)
    if(is.null(cell$._map_stamp) || cell$._map_stamp != stamp){
        parent <- as.environment(cell$parents[[pnr]])
        stop("fixme: ixs_dim are no longer there, ixs have meaning of ixs0")
        if(cell$ixs_dim[[pnr]] != 1L)
            stop.pbm("Folds cannot propagate from parents to children if ix_dim is not 1.
Not the case of ", .getType(cell, F), " <- ", .getType(parent, F))
        cell$folds <- parent$folds[cell$ixs[[pnr]],, drop = F]
        ## don't know if folds are exclusive, so reasign
        cell$.folds_exclusive <- NULL
        cell$folds_factor <- NULL
        eval(d_rec_set_folds)
    }}

.class.ind <- function(cl){
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)) )
    x[(1:n) + n*(unclass(cl)-1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
}

.ind2factor <- function(class_ind){
    storage.mode(class_ind) <- "logical"
    out <- character(NROW(class_ind))
    names <- colnames(class_ind)
    for(d in 1:length(names))
        out[class_ind[, d]] <- names[d]
    as.factor(out)
}
        


.field.current_mirror <- function(mname){
    if(missing(mname))
        current_mirror
    else{
        . <- as.environment(.self)
        if(is.null(mname))
            mname <- "base"
        else if(!mname %in% names(mirrors))
            stop("Mirror ", mname, " has not been defined.")
        .$mirrors[[.$current_mirror]] <- .$.cells
        .$current_mirror <-  mname
        .$.cells <-  mirrors[[mname]]
        mname
    }}

## .field.train_fold <- function(fld){
##     if(missing(fld))
##         ## todo: infer from predict_ix if not found
##         get("train_fold", .self)
##     else{
##         if(is.null(folds))
##             stop("FOLDS not defined")
##         if(!as.character(fld) %in% rownames(folds))
##             stop("FOLD must one of ", paste(rownames(folds), collapse = ", "))
##         assign("train_fold", fld, .self)
##         .create_folded_mirror(folds[, fld])
##         ## assign("predict_ix", !train_ix, .self)
##         .self
##     }}

