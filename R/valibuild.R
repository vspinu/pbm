### BASIC, size, varnames, vardim
.validate_protocol_size <- function(S, .self, msg){
    ## S is a function, or numeric of length 1 or 2
    if(!is.null(S))
        if(is.function(S))
            S(.self)
        else if (is.numeric(S))
            if (length(S) == 1L){
                if(S != .self$size)
                    stop.pbm("Protocol.size validation failed", msg, " expected, ", S, ", found ", .self$size)
            }else if (length(S) == 2){
                s <- .self$size
                if( !(S[[1]] <= s && s <= S[[2]]) )
                    stop.pbm(sprintf("Protocol.size validation failed %s: expecting %s <= size <= %s, got %s",
                                     msg, S[[1]], S[[2]], s))
                else warn.pbm("Don't know how to handle SIZE protocol. Skipping.")
            }}

.validate_protocol_varnames <- function(VN, .self, msg)
    if(!is.null(VN)){
        if(is.function(VN))
            VN(.self)
        else if(is.list(VN)){
            varnames <- .self$varnames
            for(i in seq_along(VN))
                if(!is.null(VN[[i]]) && !identical(varnames[[i]], VN[[i]]))
                    stop.pbm("Protocolo.varnames failed ", msg, ": dim[", i, "] is expected to be: (",
                             paste(VN[[i]], collapse = ", "), ")\n, but found:  (",
                             paste(varnames[[i]], collapse = ", "),")")
        }else
            warn.pbm("Don't know how to handle varnames protocol. Skipping.")
    }

    
.validate_protocol_vardim <- function(VD, .self, msg = "")
    if(!is.null(VD)){
        if(is.function(VD))
            VD(.self, msg)
        else if(is.numeric(VD)){
            vd <- .self$vardim
            if(length(VD) != length(vd))
                stop.pbm("vardim length protocol failed", msg, ": expecting length (", length(VD) ,
                         "), found (", length(vd), ")")
            for(i in seq_along(VD))
                if(!is.na(VD[[i]]) && VD[[i]] != vd[[i]])
                    stop.pbm(" vardim protocol failed ", msg, ": dim[",  i,
                             "] is expected to be: (", VD[[i]], ")\n found:  (", vd[[i]], ")")
        }}
        


### PARENTS
.validate_protocol_parents <- function(P, .self){
    parents <- .self$parents
    if(!exists("parents", .self, inherits=FALSE) || length(parents)==0)
        stop.pbm(" 'parents' is not found in this BCell.")
    ##0 parents exist (if one validates parents, it expects some parents?)
    if(length(P) > 0 &&
       length(parents)!=length(P)) ##probably is not needed (done in initialization); problem with "..."
        stop.pbm("Parents validation failed: number of 'parents' in protocol is ", length(protocol$parents), ", supplied ", length(parents))
    ##1 nr_parents == length(parents) [ will not work with ... tothink]
    ##2 by construction in bcell names are ok and parents are sorted in accordance to protocol specification.
    for(n_pr in names(parents)){
        .validate_protocol_size(.eval_maybe(P[[n_pr]]$size, .self),
                                parents[[n_pr]], sprintf("in parent %s", n_pr))
        .validate_protocol_vardim(.eval_maybe(P[[n_pr]]$vardim, .self),
                                  parents[[n_pr]], sprintf("in parent %s", n_pr))
        .validate_protocol_varnames(.eval_maybe(P[[n_pr]]$varnames, .self),
                                    parents[[n_pr]], sprintf("in parent %s", n_pr))
    }
}

.eval_maybe <- function(obj, .self){
    if(is.language(obj)) eval(obj, .self)
    else obj
}

.validate_basic_objects_exist <- function(.self)
    for( nm in c("st", "ll", "size"))
        if(!exists("st",.self, inherits = F))
            stop.pbm(sprintf("%s was not initialised"))

d_validate.lldim <- form(
  if(is.null(lldim)){
      if(length(dim(st)) < length(dim(ll)))
          stop.pbm("number of LL dimensions is bigger than number of ST dimmenstions")
  }else{
      if(lldim != length(dim(ll)) - 1L)
          stop.pbm(sprintf("LLDIM (%s) is not equal to number of dimensions of LL (%s)",
                           lldim, length(dim(ll))))
  })

d_validate.children <- form()



### IXS
d_build.ixs <- form({
    for(i in names(parents)){
        if(!is.null(ixs_dim[[i]]) && !ixs_dim[[i]] %in% 0:2)
            ## API: ixs_dim is in 0:2
            stop.pbm("ixs_dim[[\"",i,"\"]] should be in 0:2")

        if(is.null(ixs[[i]]))
            if(!is.null(ixs_dim[[i]]))
                ## API: when cix is NULL and cix_dim is given, build ix = 1:N
                ixs[[i]] <- # validation comes later
                    seq_len(switch(as.character(ixs_dim[[i]]),
                                   "0" = length(ll), "1" = size, "2" = varsize))
            else if(parents[[i]]$size == 1L)
                ## API: when nor cix, nor cix_dim are givven, but parent size is 1, replicate 1
                ixs[[i]] <- 1L
            else
                stop.pbm(sprintf("cix (%s) was not supplied, nor ix_dim", i))

        if(!is.language(ixs[[i]])){
            ## if ixs_dim is given replicate along oposite dimmension
            if(!is.null(ixs_dim[[i]])){
                ## derive ixs from ix_dim
                ixs[[i]] <-
                    switch(as.character(ixs_dim[[i]]),
                           "0", 
                           "1" = rep_len(ixs[[i]], length(ll)), 
                           "2" = c(matrix(ixs[[i]], nrow = size, ncol = ncol(ll), byrow = T)))
            }
            ixs[[i]] <- rep_len(as.integer(ixs[[i]]), length(ll))
        }
    }
    ix <- ixs[[i]]
})

d_validate.ixs <- form({
    if(length(ixs)!=length(parents)) stop.pbm(" length of 'ixs' and 'parents' differ")
    for(i in names(parents)){
        ## ixs
        ## - integer check
        ## - unique == 1:dim(ll)[ixs_dim[[i]]]
        ## - length ixs[[i]]== dim(parents_st)[[1]]
        if(!is.language(ixs[[i]])){
            if(!identical(as.integer(sort(unique(ixs[[i]]))), 1:parents[[i]]$size))
                stop.pbm("ixs[[", i, "]] must be in the range 1:", parents[[i]]$size, 
                         " to agree with SIZE in parent '", i, "'")
            else
                if(length(ixs[[i]]) != length(ll))
                    stop.pbm(" length of ixs[[", i, "]] (", length(ixs[[i]]),
                             ") must be equal to length of ll  (", length(ll), ")")
            }}})

d_validate.pix <- form({
    for(pname in names(parents)){
        ## API: ixs should be of the same length as ll
        if(length(pix(, pname)) != length(ll))
            stop.pbm(sprintf("length (%s) of index returned by pix('%s') does not match the length of ll (%s)",
                             length(pix(, pname)), pname, length(ll)))
    }})

## API: PV validate only names declared in parnames (forms derived from parent's cpars not checed)
d_validate.PV <- form({
    for(vname in parnames){
        if(length(PV(vname)) != length(ll))
            stop.pbm(sprintf("length (%s) of indexes returned by PV('%s')' does not match the length of ll (%s)",
                             length(PV(vname)), vname, length(ll)))
    }})
