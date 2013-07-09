library(coda)
### extension of coda mcmc for 2D variables
pbmmc <- function (data = NA, start = 1, end = NULL, thin = 1) 
{
    pbm_class <- "pbmmc"
    if (is.array(data)) {
        if (length(dim(data)) > 3L) {
            stop("arrays with more than 3 dimmensions are not supported")
        }else if (length(dim(data)) == 3){
            pbm_class <- "pbmmc2d"
            orig.dim <- dim(data)
            orig.dn <- dimnames(data)
            vnames <- orig.dn[[3]]
            if(!is.null(vnames))
                vnames <- rep(vnames, each = orig.dim[[2]])
            cnames <- orig.dn[[1]]
            if(is.null(cnames))
                cnames <- as.character(1:orig.dim[[2]])
            cnames <- rep.int(cnames, orig.dim[[3]])
            dim(data) <- c(orig.dim[1], prod(orig.dim[-1]))
            colnames(data) <- paste0(vnames, cnames)
            attr(data, "orig.dim") <- orig.dim
            attr(data, "orig.dn") <- orig.dn
            if(is.null(end))
                end <- orig.dim[[1]]
        } # else, guaranted to be a matrix
    }else{
        if(is.null(end))
            end <- length(data)
    }
    data <- coda::mcmc(data, start = start, end = end, thin = thin)
    class(data) <- c(pbm_class, class(data))
    data
}


`[.pbmmc2d` <- function (x, i, j, z, drop = FALSE){
    xstart <- start(x)
    xthin <- thin(x)
    xend <- end(x)
    if(!missing(z)){
        orig.dim <- attr(x, "orig.dim")
        if(is.null(orig.dim)) stop("not a valid pbmmc2d object, attribute orig.dim is missing")
        orig.dn <- attr(x, "orig.dn")
        if(is.null(orig.dn)) stop("not a valid pbmmc2d object, attribute orig.dn is missing")
        dim(x) <- orig.dim
        dimnames(x) <- orig.dn
        y <- unclass(x)[i, j, z, drop = drop]
    }else
        y <- unclass(x)[i, j, drop = drop]
    if (length(y) == 0 || is.null(y)) 
        return(y)
    if (missing(i)) 
        return(pbmmc(y, start = xstart, end = xend, thin = xthin))
    else return(y)
}
