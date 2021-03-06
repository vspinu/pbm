adRoot <- mixin(
    initFields = list(
        do.adapt = TRUE,
        adapt_sims = 1000L),
    initForms = list(
        UPDATE.main.adapt = form(
            if(do.adapt) e(adapt))))

adHST <- mixin(
    ## Algorithm proposed by Haario, Heikki, Eero Saksman, and Johanna
    ## Tamminen. 2001. An Adaptive Metropolis Algorithm. Bernoulli 7, no. 2 (April
    ## 1): 223-242. doi:10.2307/3318737 http://www.jstor.org/stable/3318737 and
    ## implemented as Algorithm 2 (p. 358) in Andrieu, Christophe, and Johannes
    ## Thoms. 2008. A tutorial on adaptive MCMC. Statistics and Computing 18, no. 4
    ## (12): 343-373. doi:10.1007/s11222-008-9110-y
    ## http://www.springerlink.com/content/979087678366r78v/.
    initFields = list(
        ## "scale" field is already defined
        abeta = 1, # todo: make a check of beta > 0
        aalpha = .5, # todo: make a check for alpha \in (0, 1])
        do.group_scale = TRUE, # whether to compute scale per parent ix
        do.mc_scale = FALSE,
        mc_scale = array(double()), 
        dilate_scale = 1), # multiplicative factor to adjust scale 
    initForms = list(
        init.M.build.ad = form(
            st = form({
                aix <- do.call(base::interaction,
                                 c(lapply(names(parents), function(p) pix(, p)),
                                   list(drop = TRUE)))
                anr_gr_ix <- c(table(as.integer(aix)))
                e(set.ad.st)
            }),
            mu = form({
                ## mu is never grouped
                amu <- array(0, dim = dim(st), dimnames = dimnames(st))
            }),
            sigma = form({
                alambda <- 2.38^2/ncol(st) # = varsize
                asigma <- rep.int(scale^2L, length(st))
                aoffset <- (0.0001/ncol(st))*asigma
            })),
        init.R.build.ad = form(
            rebuild_maybe = form(if(!exists("aix", inherits = F)) e(init.M.build.ad)), 
            mc = form({
                if(do.mc_scale){
                    if(length(mc_scale) == 0L){
                        mc_scale <- array(dim=c(0, dim(st)), dimnames=c(list(mc=NULL), dimnames(st)))
                    }
                    ._prev_size <- nrow(mc_scale)
                    mc_scale <- abind(mc_scale, array(NA_real_, dim=c(.nr_iter, dim(mc_scale)[ -1])), along=1)
                    attr(mc_scale, "prev_size") <- ._prev_size
                }})),
        set.ad = form(
            st = form(ast <- mh_tr@ITR(st)),
            gamma = form(agamma <- abeta/.N^aalpha),
            mu = form({
                amu <- amu + agamma*(ast - amu)
            }),
            sigma = form({
                asigma <- asigma + agamma*(((ast - amu)^2L) - asigma)
            }),
            scale = form({
                scale <-
                    if(do.group_scale)
                        sqrt(rowsum(c((dilate_scale*alambda)*asigma + aoffset), aix)/anr_gr_ix)[aix]
                    else
                        sqrt(c((dilate_scale*alambda)*asigma + aoffset))
            })),
        adapt = form({
            e(set.ad.st)
            e(set.ad.gamma)
            e(set.ad.sigma) ## sigma should be first,  depends on previous amu!!
            e(set.ad.scale)
            e(set.ad.mu)
            if(do.mc_scale){
                if(do.mc_scale) mc_scale[.I + attr(mc_scale, 'prev_size'),, ] <- c(scale)
            }
        })),
    parentMixins = adRoot,
    subtype = "HST")

