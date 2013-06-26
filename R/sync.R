.redirect_kins_toMe <- function(me, container){
    ## make children and parents of ALL cells in container
    ## to point to cell ME and add/correct a reversal link if needed (sync)
    ## i.e. make all the children and parents of type .getType(me) to really be the ME cell.
    ## invoked for side effect
    type <- .getType(me)
    if(exists(type, container, inher = FALSE)){
        lapply(ls(container, all.names = TRUE), function(nm){
            cellEnv <- as.environment(get(nm, envir = container, inherits = FALSE))
            if(length(children <- cellEnv[["children"]]))
                for(i in seq_along(children))
                    if(.getType(children[[i]]) == type){
                        cellEnv[["children"]][[i]] <- me
                        .synchronize_child(cellEnv[["children"]][[i]], me)
                    }
            if(length(parents <- cellEnv[["parents"]]))
                for(i in seq_along(parents))
                    if(.getType(parents[[i]]) == type){
                        cellEnv[["parents"]][[i]] <- me
                        .synchronize_parent(cellEnv[["parents"]][[i]], me)
                    }
        })
    }else{
        warn.pbm(gettextf("cell %s is not installed in the container of type %s; no parents / children synchronization has been performed.",
                          type, container@typeContainer))
    }
    invisible(NULL)
}

.synchronize_parent <- function(me, parent){
    ## checks if type of ME is present in children of PARENT,
    ## if not,  inserts ME in the "children" bclist of PARENT
    ## if exists,  checks if it is the same object,  if not,
    ## re-point the PARENT'S child to "right" child ME
    stopifnot(is(parent, "BC"))
    parentEnv <- as.environment(parent)
    me_type <- .getType(me)
    changed <- FALSE
    len <- length(children <- parentEnv[["children"]])
    i <- 1
    found_type <- FALSE
    while(i <= len && !found_type){
        if(.getType(children[[i]]) == me_type){
            found_type <- TRUE
            if(!identical(as.environment(children[[i]]), me)){
                ## replace if not identical
                parentEnv[["children"]][[i]] <- me
                changed <- TRUE
            }
        }
        i <- i + 1
    }
    if(!found_type){
        parentEnv[["children"]][[get(".type", envir = me)]] <- me
        changed <- TRUE
    }
    if(changed)
        assign("child", parentEnv[["children"]][[1L]], envir = parentEnv)
    invisible(changed)
}

.synchronize_child <- function(me, child){
    ## checks if type of ME is present in parents of CHILD,
    ## if not,  inserts ME in the "parents" bclist of CHILD
    ## if exists,  checks if it is the same object,  if not,
    ## re-point the CHILD'S parent to "right" parent ME
    stopifnot(is(child, "BC"))
    childEnv <- as.environment(child)
    me_type <- .getType(me)
    changed <- FALSE
    len <- length(parents <- childEnv[["parents"]])
    i <- 1
    found_type <- FALSE
    while(i <= len && !found_type){
        if(.getType(parents[[i]]) == me_type){
            found_type <- TRUE
            if(!identical(as.environment(parents[[i]]), me)){
                ## replace if not identical
                childEnv[["parents"]][[i]] <- me
                changed <- TRUE
            }
        }
        i <- i + 1
    }
    if(!found_type){
        childEnv[["parents"]][[get(".type", envir = me)]] <- me
        changed <- TRUE
    }
    if(changed)
        assign("parent", childEnv[["parents"]][[1L]], envir = childEnv)
    invisible(changed)
}


###_* BCLIST,  PARENTS,  CHILDREN
## .insert_in_bclist <- function(bclist, bc){
## ## does not insert if already in the bclist, even if with different name.
## ## overwrites the existing cell with the name bc[["type"]]
##     if(!is(bc, "BC")) stop("\npbm: 'bc' must be of class `BC`.")
##     bc_in_bclist <- unlist(lapply(bclist, function(el)
##                                   identical(as.environment(el),
##                                             as.environment(bc))))
##     if(is.null(bc_in_bclist)||!any(bc_in_bclist))
##         bclist[[bc[["type"]]]] <- bc
##     bclist
## }

## .syncronize_bclist <- function(bclist, info=TRUE){
##     bclist <- as.bclist(bclist)
##     changed <- rep.int(TRUE, length(bclist))
##     while(any(changed)){
##         for(i in seq_along(bclist)){
##             changed[i] <- any(.syncronize_parents(bclist[[i]]), .syncronize_children(bclist[[i]]))
##             if(info) cat(" * Synced ", i, bcType(bclist[[i]]),":", !changed[i], fill=TRUE)
##         }}
##     bclist
## }

## .syncronize_parents <- function(bcell){
##     stopifnot(is(bcell, "BC"))
##     changed <- FALSE
##     bc <- as.environment(bcell)
##     lapply(bc[["parents"]], function(p){
##         changed <<- .syncronize_parent(bcell, p) || changed
##     })
##     return(changed)
## }

## .syncronize_children <- function(bcell){
##     stopifnot(is(bcell, "BC"))
##     changed <- FALSE
##     bc <- as.environment(bcell)
##     lapply(bc[["children"]], function(p){
##         changed <<- .syncronize_child(bcell, p) || changed
##     })
##     return(changed)
## }

## .redirect_my_kins <- function(me, container){
##     ## redirect the kins of me to point to right type from the container
##     ## i.e. children[[i]] of type "aa.bb.cc" should be identical to container[["aa.bb.cc"]]
##     ## if type of a children is not present skip with a message.
##     ## used for side effect.
##     if(length(children <- get("children", envir = me,  inher = FALSE)))
##         for(nm_c in names(children))
##             if(exists(c_type <- .getType(children[[nm_c]]), envir = container,  inher = FALSE))
##                 me$children[[nm_c]] <- get(c_type, envri = container) ## fixme:  use of $operator here!!
##             else
##                 message(gettextf("child of type %s is present in the cell %s but not in the container %s "
##                                  c_type, .getType(me), container@typeContainer))
##     if(length(parents <- get("parents", envir = me,  inher = FALSE)))
##         for(nm_c in names(parents))
##             if(exists(c_type <- .getType(parents[[nm_c]]), envir = container,  inher = FALSE))
##                 me$parents[[nm_c]] <- get(c_type, envri = container)  ##fixme:  use of $operator here.
##             else
##                 message(gettextf("parent of type %s is present in the cell %s but not in the container %s "
##                                  c_type, .getType(me), container@typeContainer))
##     invisible(NULL)
## }
