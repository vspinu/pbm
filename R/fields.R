## API: size, varsize, st are primitives. They re-size each other and other
## state vars. varsize == dim(st)[[2]], size == dim(st)[[1]]
##
## API: ll, varnames, names are not primitives and throw errors in case of
## dimension mismatch
##
## API: varnames, names, var are abstract fiedls (not explicitly stored)
.reinit_state <- function(.self){
    for( nm in c("size", "varsize", "st", "ll")){
        if (exists(nm, .self, inherits = F))
            remove(list = nm, envir = .self)
    }}

.field.size <- function(arg){
    ## API: size is primitive
    ## API: size<- changes everything
    ## fixme: check: how it plays with tr
    if (missing(arg))
        return(size)
    else{
        if(!is.numeric(arg)) stop("supplied SIZE is not numeric")
        if(length(arg) > 1) stop("SIZE should be a number")
        assign("size", arg, .self)
        assign("._size_init_done", TRUE, .self)
        ## very silent R type recycling
        if(dim(ll)[[1]] != size)
            assign("ll", array(-Inf, c(size, ncol(ll))), .self)
        if(dim(st)[[1]] != size){
            if(length(dim(st)) == 1L) stop.pbm("Oups, st dimnum is 1. Should not have happened.")
            assign("st", st[rep(1:(dim(st)[[1]]), length.out = size),, drop = FALSE], .self)
        }
        .self
    }}

.field.varsize <- function(arg){
    ## API: varsize is primitive fields
    ## API: varsize can affect only ST (not SIZE, nor LL)
    if(missing(arg))
        return(varsize)
    else{
        assign("varsize", arg, .self)
        if(!identical(dim(st)[-1], varsize)){
            ## API: llcols is varsize if .lldim is not supplied (NULL)
            if(is.null(.lldim)) llsize <- varsize
            else{
                ## API: llcols is minimum of .lldim and varsize
                llsize <- min(.lldim, varsize)
                if(.lldim != length(llsize)){
                    assign(".lldim", length(llsize), .self)
                    warn.pbm("LLDIM was reset to ", .lldim)
                }}
            ## standard ST replication by columns
            assign("ll", array(-Inf, c(size, llsize)), .self)
            assign("st", array(st, c(size, varsize)), .self)
        }
        .self
    }}

## tt$initFields(v = protoField(
.field.var <- function(arg){
    ## API: var is abstract field derived from st[1, ]
    ## API: When var<- is always replicated across first dimension
    ## API: var<- can not change SIZE, but can change varsize and LL
    ## fixme: check: how it plays with tr
    if(missing(arg)){
        out <- st[1,, drop = F]
        if(size > 1) rownames(out) <- NULL
        out
    }else{
        if (is.null(arg)){
            .reinit_state(.self)
            return(.self)
        }
        rv <- as.array(arg)
        if(length(dim(rv)) == 1L)
            rv <- t(rv)
        varsize <- dim(rv)[-1]
        nr_dim <- length(varsize) + 1L
        st <- aperm(array(rv, dim = c(varsize, size)), c(nr_dim, 1:(nr_dim-1L)))
        if(!is.null(dimnames(rv)))
            dimnames(st) <- c(list(NULL), dimnames(rv)[-1])
        if(is.null(.lldim)) llsize <- varsize #uninitialised by the user or protoparent
        else{
            llsize <- varsize[1:min(.lldim, length(varsize))]
            if(.lldim != length(llsize)){
                assign(".lldim", length(llsize), .self)
                warn.pbm("LLDIM was reset to ", .lldim)
            }}
        assign("ll", array(-Inf, c(size, llsize)), .self)
        assign("varsize", varsize, .self)
        assign("st", st, envir=.self)
        .self
    }}

.field.st <- function(arg){
    ## API: st is concrete primitive fields
    if(missing(arg)){
        return(st)
    }else{
        arg <- as.array(arg)
        if(!(._size_init_done || ._var_init_done)){
            ## give power to ST for the first time only, useful for DataCells (DC)
            nm <- NULL
            if(length(dim(arg)) == 1L){
                nm <- names(arg)
                dim(arg) <- c(length(arg), 1L)
            }
            ## fixme: check: how it plays with tr
            assign("st", arg, .self)
            if(!is.null(nm))
                .self$names <- nm
            .self$size <- dim(st)[[1L]]
        }else{
            if(length(st) %% length(arg) != 0)
                stop.pbm("Supplied vector cannot be evenly recycled to fit ST object")
            st[] <- arg
            assign("st", st, .self)
        }
        if(length(rownames(arg)) == size)
            .self$names <- rownames(arg)
        if(identical(dim(arg), dim(st)) && !is.null(dimnames(arg)))
            .self$varnames <- dimnames(arg)[-1]
    }
}

.field.varnames <- function(arg){
    ## API: varnames is abstract fields
    if(missing(arg)){
        names <- dimnames(st)[-1]
        names <-
            if (length(names) == 1L) names[[1]]
            else names
        if(is.null(names))
            ## could never be null
            rep.int("", varsize)
        else names
    }else{
        if (is.null(arg) ){
            rownames <- rownames(st)
            dimnames(st) <- NULL
            rownames(st) <- rownames
        }else{
            if(is.character(arg)) arg <- list(arg)
            else if(!is.list(arg))
                stop.pbm("varnames must be a character vector or a list of character vectors.")
            if(length(arg) != length(dim))
                stop.pbm("varnames dimension should match the length of varsize (", length(varsize), ")")
            ## don't copy st localy on dimnames assignment
            eval(substitute(dimnames(st) <- c(list(dimnames(st)[[1]]), arg), list(arg = arg)),
                 .self)
            .self
        }}}

.field.ll <- function(arg){
    ## API: ll is concrete primitive field
    ## dimnames are not sticky! new ARG replaces everything (ie. dim and dimnames)
    if(missing(arg)){
        return(ll)
    }else{
        stop.pbm("Cannot assign LL. LL is initialised automatically to -Inf.
You can control its size with LLDIM field.")
    }}

.field.lldim <- function(arg)
    ## API: lldim if not set means lldim == varsize (fixme: name it llcols)
    ## fixme: why .lldim is assigned, cannot do without?
    ## because distributions like Multivar sets it to 1
    ## fixme: rename into llsize?
    ## fixme: doesn't take effect if precedes var field assignment
    ## fixme: lldim is not inherited (see log.dganorm)
    if(missing(arg)){
        if(length(out <- dim(ll)[-1]) > 0)
            out
        else 1L
    }else{
        arg <- as.integer(arg)
        if(length(arg) != 1L) stop.pbm("LLDIM must be a number")
        if(arg > length(varsize))
            stop.pbm(sprintf("LLDIM (%s) cannot be bigger than varsize (%s)", arg, varsize))
        assign(".lldim", arg, .self)
        assign("ll", array(-Inf, c(size, arg)), .self)
    }

.field.names <- function(arg){
    ## API: names is abstract fields derived from rownames(st)
    ## fixme: if called before size or st fields, give error, order fields initialization somehow
    if(missing(arg)){
        rownames(st)
    }else{
        eval(substitute({
            rownames(st) <- nm
            rownames(ll) <- nm
        }, list(nm = arg)), .self)
        .self
    }}

## todo, this is a commong case, make argument to protoField
.field.sims <- function(el){
    if(missing(el)){
        .cells@sims
    }else{
        stop("Object 'sims' is read-only")
    }}

.field.children <- function(children){
    if(missing(children)){
        get("children", envir=.self)
    }else{
        assign("children", as(children, "bcChildren"), envir=.self)
        for(ch in children)
            .syncronize_child(.self, ch)
        .self
    }
}

.field.parents <- function(parents){
    if(missing(parents)){
        get("parents", envir=.self)
    }else{
        assign("parents", as(parents, "bcParents"), envir=.self)
        for(pr in parents)
            .syncronize_parent(.self, pr)
        .self
    }
}

.field.pc_ll_groups <- function(gr)
    .field.xx_ll_groups(gr, "pc", .self)

.field.mc_ll_groups <- function(gr)
    .field.xx_ll_groups(gr, "mc", .self)

.field.xx_ll_groups <- function(gr, xx = "mc", .self){
    xx_ll_groups <- sprintf("%s_ll_groups", xx)
    xx_ll <- sprintf("%s_ll", xx)

    if(missing(gr)){
        get(xx_ll_groups, .self)
    }else{
        if(is.character(gr)){
            values <- c("all", "folds", "sum", "rowSums", names(get("parents", .self)))
            if(!(gr %in% values))
                stop(sprintf("character %s should be one of ", xx_ll_groups),
                     paste0("'", values, "'", collapse = ", "))
            if(identical(gr, "folds") && is.null(folds))
                stop.pbm(sprintf("folds are not set, cannot use as %s", xx_ll_groups))
        }else if (is.factor(gr) || is.integer(gr)) {
            if (length(gr) != size)
                stop(sprintf("if factor, %s must have the same length [%s] as cell size [%s]",
                             xx_ll_groups, length(gr), size))
            gr <- as.factor(gr)
        }else stop(sprintf("%s should be either character, integer or factor", xx_ll_groups))
        assign(xx_ll, array(double()), .self)
        warn.pbm("reseting mc_ll")
        assign(xx_ll_groups, gr, .self)
    }
}

.field.folds <- function(folds){
    if(missing(folds))
        get("folds", .self)
    else{
        ## stop("Cannot assign folds into model object, use FOLDS on a cell instead.")
        if(length(folds) != size)
            stop(sprintf("ll and folds length differ (%d vs %d)",
                         size, length(folds)))
        .set_folds(.self, folds)
        .self
    }}

.field.mfolds <- function(folds){
    ## derive at runtime from mfolds_factor
    if(missing(folds)){
        if(!exists("folds_factor", .self, inherits = F))
            .set_folds_factor(.self)
        out <- .class.ind(get("folds_factor", .self))
        storage.mode(out) <- "logical"
        out
    }else{
        stop("mirror 'folds' field is readonly.")
    }}

.field.folds_factor <- function(folds){
    if(missing(folds)){
        if(!exists("folds_factor", .self, inherits = F))
            .set_folds_factor(.self)
        get("folds_factor", .self)
    }else{
        stop("'fold_factor' is readonly.")
    }}

.field.folds_names <- function(obj){
    if(missing(obj)){
        if(!exists("folds"))
            stop.pbm("'folds' have not been set in this cell")
        colnames(folds)
    }else{
        stop("'fold_names' field is readonly.")
    }}

.field.mfolds_names <- function(obj){
    ## mfolds_names is primitive
    ## mfolds and mfold_factor are derived from it
    if(missing(obj)){
        get("mfolds_names", .self)
    }else{
        stop("mirror 'fold_names' field is readonly.")
    }}

.field.adapt <- function(name){
    if(missing(name)){
        if(!exists("adapt", .self, inherits = FALSE)){
            ## if does not exist install the default algo
            if(!exists(".ALGO", .self, inherits = FALSE))
                assign(".ALGO",
                       new("ALGO", type = "local_algorithms", rootParentEnv = .self),
                       envir = .self)
            .ALGO$initCells("adapt.*") # install the algorithm in the local context
            .ALGO$`adapt.*`$.host <- .self
            assign("adapt", .getCell("adapt.*", .ALGO), envir = .self)
        }
        get("adapt", .self, inherits = FALSE)
    }else{
        ## init .ALGO
        if(!exists(".ALGO" , envir=.self, inherits = FALSE))
            assign(".ALGO", new("ALGO", type = "local_algorithms", rootParentEnv = .self),
                   envir = .self)
        if(is.null(name)){
            remove(adapt, envir = .self) ## does not set do.adapt to FALSE (do.adapt can be inherited)
            return(invisible(NULL))
        }
        if(is.character(name)){
            if(length(name) == 0L || name == "none")
                name <- "adapt.*"
            type_long <- .getType(.getPartial(name, .ALGO[[".cells"]], object_name = "adaptation algorithm"))
            if(!grep("adapt", type_long, fixed = TRUE))
                stop("Cell ", type_long, " is not an 'adapt' algorithm")
            .ALGO$initCells(name) # install the algorithm in the local context
            .cell <- .getCell(name, .ALGO)
            .cell$.host <- .self
            assign("adapt", .cell, envir = .self)
            assign("do.adapt", TRUE, envir = .self)
            if(.sims$end > 1) ## needs reinitialization only if is M initialized
                assign("do.init_M_in_R", TRUE, .cell)
            ## only if is.character try to do the reinitialization
        }else if(is(name, "algo")){
            ##todo: check if it is the same algo from .ALGO container
            type_long <- .getType(name)
            if(!grep("adapt", type_long, fixed = TRUE))
                stop("Cell ", type_long, " is not an 'adapt' algorithm")
            if(.sims$end > 1) ## needs reinitialization only if is M initialized
                assign("do.init_M_in_R", TRUE, .cell)
            assign("adapt", name, envir = .self)
            assign("do.adapt", TRUE, envir = .self)
        }else{
            stop.pbm("Only objects of class 'character' or 'algo' could be assigned to 'adapt' fields;\n supplied an object of class \"", class(name), "\"")
        }
    }
}

.field.scale <- function(x){
    if(missing(x)){
        get("scale", envir = .self)
    }else{
        ## scale is for st not ll
        scale <- as.numeric(x)
        if(length(scale) == length(st))
            if(!identical(dim(scale), dim(st)))
                dim(scale) <- dim(st)
        else if(!(length(scale) == 1 || length(scale) == size))
            stop("SCALE parameters must be of length 1, ",
                 length(ll), "(length(ll)) or ", size, " (size)")
        assign("scale", as.array(scale), .self)
    }
}


.make_field_mc <- function(mcname){
    protoField(eval(substitute(
        function(x)
        if(missing(x)){
            pbmmc(get(mcname, .self))
        }else{
            stop(mcname, "is read only")
        })))
}


.field.st_tr <- function(x){
        if(missing(x))
            get("st_tr", .self)
        else{
            if(!is(x, "protoTransform"))
                stop("tr fields should be of class 'protoTransform'; supplied ", class(x))
            assign(".subtypes", c(x@name, .subtypes))
            assign("st_tr", x, .self)
            .self
    }
}
