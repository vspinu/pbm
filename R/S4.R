setCellClass("BC", contextClass = "PBM")
setClass("bclist", list(names = "character"), contains = "list")
setClass("bcParents", contains = "bclist")
setClass("bcChildren", contains = "bclist")
setClass("parentDefinition", contains = "namedList")
setClass("bCellContainer", representation(sims = "namedList", is_mirror = "logical", train_mirror_name = "character"),
         prototype = list(is_mirror = FALSE, sims = new("namedList", list(start = 1L, end = 0L))),
         contains = "cellContainer")

instantiate_parents <- function(parents, homeContext = NULL){
    ## parents is a list of parents defenitions
    ## RETURN: list of parents (cells!!), ixs and ixs_dim all with the same names

    if(any((not_P <- !(sapply(parents, class) %in% c("parentDefinition", "character")))))
        stop("parents arg should contain objects of class 'parentDefinition' or 'character' only")

    out <- list(parents = list(), cixs = list(), cixs_dim = list(), cpars = list())
    names <- allNames(parents)
    for(i in seq_along(parents)){
        ## generate the bcell
        .p <- parents[[i]]
        parent <- .p$parent
        pn <-
            if(names[[i]] != "")
                names[[i]] ## external names in defBC have precedence
            else if(is(.p$parent, "protoCell"))
                .getType(.p$parent)
            else if(is.character(parent))
                parent
            else stop("cannot infer parent's name")

        if(is.character(parent)){
            ## GET from homeContext
            if(is.null(homeContext)){
                stop("Cannot create bcell from string /'", parent, "' when homeContext is NULL")
            }else{
                ## previously declared "parent" are accessed here
                .cells <- get(".cells", envir = homeContext)
                typeLong <- protoClasses:::.complete_partial_name(parent, .cells)
                if(!is.character(typeLong))
                    stop("Cannot find bcell with (partial) name '", parent, "', or match is not unique.")
                if(exists(typeLong, envir = .cells, inherit = FALSE)){
                    ## GET existing one from current context
                    parent <- get(typeLong, envir = .cells)
                    pn <- .getType(parent, F)
                }else{
                    ## CLONE the existing one from inhereted context (rarely needed)
                    parent <- new("BC", type = typeLong, homeContext = homeContext) ## todo: this fails for sure
                    pn <- .getType(parent, F)
                }
            }
        }else if(is(parent, "protoCellDefinition")){
            ## GENERATE FROM INFO
            loc_type <- parent[["type"]]
            if(is.null(loc_type) || loc_type == "--" )
                parent[["type"]] <- pn
            parent <- cellFromDefinition(parent, homeContext = homeContext)
        }else if(is(parent, "protoCell")){
            ## SIMPLY RETURN
            if(identical(.getType(parent, F), "--"))
                protoClasses:::.setType(parent, pn)
        }else{
            stop("Parent in object of class '.p', should be of class character,  protoCellDefinition, or protoCell;",
                 "\nSomething went wrong,  please report")
        }
        out[["parents"]][[pn]] <- parent
        out[["cixs"]][[pn]] <- .p[["cix"]]
        out[["cixs_dim"]][[pn]] <- .p[["cix_dim"]]
        out[["cpars"]][[pn]] <- .p[["cpars"]]
    }
    out
}

instantiate_children <- function(children, homeContext = NULL){
    ## children is a character vector of already instantiated cells
    ## RETURN: list of children (cells!!), ixs and ixs_dim all with the same names
    out <- list()

    if(any((not_C <- sapply(children, class)) != "character"))
        stop("children arg should contain objects of class 'character',",
             "representing already instantiated cells\n",
             "Not true for:", paste(which(not_C), collapse = ", "), "/")

    for(nm in children){
        .cells <- get(".cells", envir = homeContext)
        typeLong <- protoClasses:::.complete_partial_name(nm, .cells)
        if(!is.character(typeLong))
            stop("Cannot find bcell with (partial) name '", nm, "'")
        if(exists(typeLong, envir = .cells, inherit = FALSE))
            ## GET existing one from current context
            out[[nm]] <- get(typeLong, envir = .cells)
        else stop(sprintf("child cell of type `%s` was not yet instantiated", typeLong))
    }
    out
}

setMethod("installBinding", "BC",
          function(bindDefinition, container, bindName, ... ){
              ## bindName has precedence:
              if(!missing(bindName) && !identical(bindName, ""))
                  protoClasses:::.setType(bindDefinition, bindName)
              else
                  bindName <- .getType(bindDefinition, F)
              where <- as.environment(container@host)[[".self"]]
              cell <- protoClasses:::.installCellInContainer(bindDefinition, container)
              cellEnv <- as.environment(cell)
              if(length(children <- cellEnv[["children"]])){
                  ## throws an error if some of the children do not exist
                  ## in the container, (parent - orientedness shows here)
                  if(any(not_exist <-
                         unlist(sapply(children, function(child){
                             !exists(.getType(child), envir = container, inherit = FALSE)
                         }))))
                      stop.pbm("Children must be already installed; not true for \n",
                               paste(names(children)[not_exist], sep = ", "))
                  for(i in seq_along(children)){
                      ## ensures that point to the object from container
                      cellEnv[["children"]][[i]] <- get(.getType(children[[i]]), envir = container, inherit = FALSE)
                      ## children's parent point to right object (cell)
                      .synchronize_child(cell, cellEnv[["children"]][[i]])
                  }
                  assign("child", cellEnv[["children"]][[1L]], envir = cellEnv)
              }
              if(length(parents <- cellEnv[["parents"]])){
                  ## replacing parents if needed
                  for (i in seq_along(parents)){
                      cellEnv[["parents"]][[i]] <- installBinding(parents[[i]], container = container)
                      .synchronize_parent(cell, cellEnv[["parents"]][[i]])
                  }
                  assign("parent", cellEnv[["parents"]][[1L]], envir = cellEnv)
              }
              cell
          })

setMethod("show", "BC",
          function(object){
              callNextMethod()
              cat("\nKins: \n")
              if(!is.null(child <- object[["child"]])) child <- format(as.environment(child))
              ## cat("\n* child: ", child, "\n* ")
              print(object[["children"]])
              if(!is.null(parent <- object[["parent"]])) parent <- format(as.environment(parent))
              ## cat("\n* parent: ", parent, "\n* ")
              print(object[["parents"]])
              cat("\nHierarchy: \n")
              print_bc_rec(object, 0)
              cat("\nBasic_objects: ")
              cat("\n st: \n")
              if(length(object$st) > 100){
                  str(object$st)
              }else{
                  print(object$st)
              }
              cat("\n ll: \n")
              if(length(object$ll) > 100)
                  str(object$ll)
              else
                  print(object$ll)
              cat("\n ixs: ")
              str(object$ixs)
              ## cat(" ixs_dim: ")
              ## str(object$ixs_dim)
          })


print_bc_rec <- function(bcell, nest.lev=0, c_nest.lev=c()){
    format.bc <- function(bcell) format(as.environment(bcell))
    .string_repr <- function(bclist){
        "Returns the vector of unique string representation of bcells in bclist."
        unlist(lapply(bclist, format.bc))
    }
    c_in_c_nest.lev <- names(c_nest.lev) %in% .string_repr(bcell$children)
    if(any(c_in_c_nest.lev)){
        c_nest <- c_nest.lev[c_in_c_nest.lev]
        c_nest_diff <- diff(c(0, c_nest))
        cat(rep.int(" :  ", max(c_nest_diff[1]-1, 0)), sep="")
        cat(" ^--")
        lapply(c_nest_diff[-1], function(i){
            cat(rep.int(" ---", max(i-1, 0)), sep="")
            cat(" ^--")
        })
    }
    cat(" ", .getType(bcell, fullName = F), "       ", sub("environment: ", "", format(as.environment(bcell))), "\n",  sep="")
                                        #  cat(rep.int(" -  ", nest.lev), " |", "\n", sep="")
    if(length(bcell$parents)){
        c_nest.lev[[format.bc(bcell)]] <- nest.lev+1
        lapply(bcell$parents, print_bc_rec, nest.lev=nest.lev+1, c_nest.lev=c_nest.lev)
    }
}

setMethod("print", "bclist",
          function(x, str=FALSE, fullNames=TRUE){
              bcl <- x
              ## recursive not needed, print.bc is recursive!
              if(length(bcl)==0) return(print("Empty bclist!!"))
              if(str) str(bcl)
              print_bc <- function(bcell, nest.lev=1){
                  bc_name <- ifelse(fullNames, .getType(bcell, fullName = F),
                                    .getType(bcell))
                  cat("", bc_name, " ",
                      sub("environment: ", "", format(as.environment(bcell))), "\n")
                  if(exists("children", envir=bcell, inherits=FALSE)&&length(bcell[["children"]])){
                      str(lapply(bcell[["children"]],
                                 function(bc) paste("C", sub("environment: ", "", format(as.environment(bc))))),
                          nest.lev=nest.lev,
                          no.list=TRUE,
                          give.head=F)
                  }
                  if(exists("parents", envir=bcell, inherits=FALSE)&&length(bcell[["parents"]])){
                      str(lapply(bcell[["parents"]],
                                 function(bc) paste("P", sub("environment: ", "", format(as.environment(bc))))),
                          nest.lev=nest.lev,
                          no.list=TRUE,
                          give.head=F)
                  }
              }
              lapply(bcl, print_bc)
              invisible(NULL)
          })

setMethod("show", "bclist",
          function(object) print(object))

setMethod("print", "bcParents",
          function(x, str=FALSE, fullNames=TRUE){
              cat("parents: \n")
              callNextMethod()})

setMethod("print", "bcChildren",
          function(x, str=FALSE, fullNames=TRUE){
              cat("children:\n")
              callNextMethod()})

setMethod("initialize", "bcParents",
          function(.Object, ...){
              .Object <- callNextMethod()
              if(length(.Object)){
                  if(dup <- anyDuplicated(names(.Object)))
                      stop("Error: duplicated names in children vector (", paste(names(out)[dup], collapse=", "), ")")
                  is_bc <- unlist(lapply(.Object, is, "BC"))
                  if(any(!is_bc))
                      stop("\npbm: some elements are not BCs, cannot create a bcParents object.")
              }
              .Object
          })

setMethod("initialize", "bcChildren",
          function(.Object, ...){
              .Object <- callNextMethod()
              if(length(.Object)){
                  if(dup <- anyDuplicated(names(.Object)))
                      stop("Error: duplicated names in children vector (", paste(names(out)[dup], collapse=", "), ")")
                  is_bc <- unlist(lapply(.Object, is, "BC"))
                  if(any(!is_bc))
                      stop("\npbm: some elements are not BCs, cannot create a bcChildren object.")
              }
              .Object
          })

setMethod("clone", "BC",
          function(x, ...){
              callNextMethod(x, exclude_names = c("child", "parent"), ...)
          })

setMethod("clone", "PBM",
          function(x, ...){
              y <- callNextMethod(x, exclude_names = c("mirrors"))
              evalq(mirrors <- list(base = .cells), y)
              y
          })

defP <- function(prototype = NULL, ...,
                 cpars = c(), cix=NULL, cix_dim=NULL, cell=NULL){
    ## API: cix can be an index, or form, in latter case it will be assigned as pix.pname
    ## Just like defBC but allows for direct specification of the cell (either
    ## as character or BC)
    if(is.null(prototype))
        if(is.null(cell))
            stop("prototype or cell arguments should be specified")
        else # ignoring ...
            parent <- cell
    else # ignoring cell
        parent <- defBC(prototype = prototype, ...)
    new("parentDefinition", list(parent = parent, cpars = cpars,
                                 cix=cix, cix_dim=cix_dim))
}


## CALL SEQUENCE:
## .initCells ->
##    (-> defBC ->  protoCellDefinition)
## installBinding_protoCellDefinition ->
##    (-> initialize "BC" -> instantiate_parents + instantiate_children)
## installBinding_BC ->
##    (-> installCellInContainer)

defBC <- function(prototype, ..., distr = NULL, var = NULL, varnames = NULL,
                  varsize = NULL, size = NULL, st = NULL, lldim = NULL, parents
                  = list(), children = list(), type = "--", initFields = list(),
                  initForms = list(), initMethods = list(), setFields = list(),
                  setForms = list(), setMethods = list(), expr = NULL){

    ## Parent can be supplied as a list of PARENTS,
    ## or as individual object of class parentDefinition
    ## Return a completed object of class protoCellDefinition
    dots <- list(...)
    names <- allNames(dots)
    parents <- as.list(parents)
    children <- as.list(children)

    for(i in seq_along(dots)){
        ## interpret all proto cells
        if(is(dots[[i]], "parentDefinition")){
            tp <- dots[[i]]
            if(length(names[[i]]))
                parents[[names[[i]]]] <- tp
            else
                parents <- c(parents, tp)
            names[[i]] <- "_to_remove_"
        }
    }
    dots[names == "_to_remove_"] <- NULL

    setFields <- c(list(var = var, varsize = varsize, varnames = varnames, size = size,
                        st = st, lldim = lldim) , dots, setFields)

    nulls <- sapply(setFields, is.null)
    setFields[nulls] <- NULL

    objDef <- new("protoCellDefinition", cellClass = "BC",
                  list(prototype = prototype,
                       type = type,
                       parents = parents,
                       children = children,
                       setFields = setFields,
                       initFields = initFields,
                       initMethods = initMethods,
                       setMethods = setMethods,
                       initForms = initForms,
                       setForms = setForms,
                       mixin = mixin,
                       expr = expr))
    objDef
}


## setCellClass("AC", contextClass="pbm", contains="BC")
## setCellClass("ZC", contextClass="pbm", contains="AC")
setMethod("initialize", "BC",
          function(.Object,
                   homeContext = NULL,
                   parents = list(), ## either character or parentDefinition objects
                   children = list(), ## only character objects
                   ## capture important names
                   mixin = list(),
                   initMethods = list(), initFields = list(), initForms = list(),
                   setMethods = list(), setFields = list(), setForms = list(),
                   expr = expression(),                    
                   INIT.C = TRUE, ...){

              ## Initialize new object and match and initialize  the parents
              ## Simmilarly to .initialize_protoCell, Don't install anything into the context
              ## initCells install the cells in the context.
              .Object <- callNextMethod(.Object, homeContext = homeContext, ...)
              if(INIT.C &&  exists("init.C", envir = .Object))
                  evalq(e(init.C), envir=.Object)
              .Object$evalq({
                  .pix_v <- .pix_p <- list()
                  for(nm in c(".pst", ".PST", ".pv", ".PV"))
                      assign(nm, list(), .self)

                  pix <-
                      eval(substitute(
                          function(v, p){
                              if(missing(p)) e(vname[[v]])
                              else e(pname[[p]])
                          }, list(vname = as.name(".pix_v"), 
                                  pname = as.name(".pix_p"))))
                  
                  pst <- function(p) e(.pst[[p]])
                  PST <- function(p) e(.PST[[p]])
                  pv <- function(v) e(.pv[[v]])
                  PV <- function(v) e(.PV[[v]])
                  ## fixme: .parent should be .pcells, pcell should be parent, parents should be pcells
                  pcell <- function(v) e(.parent[[v]])
              })
              if(length(parents)){
                  P <- instantiate_parents(parents, homeContext)
                  assign("parents", new("bcParents",
                                        P[["parents"]]), envir = .Object)
                                        #class(.Object[["parents"]]) <- c("parents", "list")
                  assign("ixs_dim", P[["cixs_dim"]], envir = .Object)
                  assign("ixs", P[["cixs"]], envir = .Object)

                  objenv <-  as.environment(.Object)
                  for(pname in names(P[["parents"]])){
                      ## if cixs are forms, assign to .fpix
                      if(is.language(P[["cixs"]][[pname]])){
                          objenv$.pix_p[[pname]] <- P[["cixs"]][[pname]]
                      }
                      .assign_pst_maybe(.Object, pname)
                      ## set all pix.xxx forms provided by user
                      .assign_pix_maybe(.Object, pname)
                  }

                  for(pname in names(P[["parents"]])){
                      ## create 'pv.xx', 'PV.xxx' and 'pixv.xxx' forms
                      pars <- P[["cpars"]][[pname]]
                      if(!is.null(pars)){
                          if(is.null(names(pars))){
                              ## API: if no names, then it is a list of cpars
                              ## declared in current parent
                              if(any(notin <- !pars %in% c(.Object$parnames, .Object$multiparnames)))
                                  stop(paste(pars[notin], collapse = ", "),
                                       " are not declared parameters in cell ", .Object$type)
                              .assign_pv_maybe(.Object, pars, pname)
                          }else{
                              ## API: if at least one name in cpars, all should be named
                              if(any(!nzchar(names(pars))))
                                  stop("if cpars is a named vector, all names must be suplied")
                              .assign_pv_paired_maybe(.Object, pars, pname)
                          }}}
              }else{
                  assign("parents", new("bcParents"), .Object)
              }

              if(length(children)){
                  ptcl <- get("protocol", get(".prototype", envir = .Object))$children
                  if((len <- length(ptcl)) > 0 &&
                     len != length(children)){
                      stop(sprintf("Supplied number of children /%s/ is not equial to defined number of children /%s/ in the protocol of object of type '%s' ",
                                   length(children), len, .getType(.Object)))
                  }
                  C <- instantiate_children(children, homeContext)
                  assign("children", new("bcChildren", C), envir = .Object)
              }else{
                  assign("children", new("bcChildren"), .Object)
              }

              ## only at the end
              .mixin(mixin, .Object,  initMethods = initMethods,
                     initFields = initFields, initForms = initForms,
                     setMethods = setMethods, setFields = setFields,
                     setForms = setForms, expr = expr)

              .Object
          })
