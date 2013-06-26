adRoot <- mixin(
    initFields = list(
        do.adapt = TRUE,
        adapt_sims = 1000L),
    initForms = list(
        update.adapt = form(
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
        ad_beta = 1, ## todo: make a check of beta > 0
        ad_alpha = .5, ## todo: make a check for alpha \in (0, 1])
        do.ad.group_scale = TRUE,
        do.mc_scale = FALSE,
        mc_scale = array(double())), 
    initForms = list(
        init.M.build.ad = form(
            st = form({
                ad_ix <- do.call(base::interaction,
                                 c(lapply(names(parents), function(p) pix(, p)),
                                   list(drop = TRUE)))
                ad_nr_gr_ix <- c(table(as.integer(ad_ix)))
                e(set.ad.st)
            }),
            mu = form({
                ## mu is never grouped
                ad_mu <- array(0, dim = dim(st), dimnames = dimnames(st))
            }),
            sigma = form({
                ad_lambda <- 2.38^2/ncol(st) # = varsize
                ad_sigma <- rep.int(scale^2L, length(st))
                ad_offset <- (0.001/ncol(st))*ad_sigma
            })),
        init.R.build.ad = form(
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
            st = form(ad_st <- TR(st)),
            gamma = form(ad_gamma <- ad_beta/.N^ad_alpha),
            mu = form({
                ad_mu <- ad_mu + ad_gamma*(ad_st - ad_mu)
            }),
            sigma = form({
                ad_sigma <- ad_sigma + ad_gamma*(((ad_st - ad_mu)^2L) - ad_sigma)
            }),
            scale = form({
                scale <-
                    if(do.ad.group_scale)
                        sqrt(rowsum(c(ad_lambda*ad_sigma + ad_offset), ad_ix)/ad_nr_gr_ix)[ad_ix]
                    else
                        sqrt(c(ad_lambda*ad_sigma + ad_offset))
            })),
        adapt = form({
            e(set.ad.st)
            e(set.ad.gamma)
            e(set.ad.sigma) ## sigma should be first,  depends on previous ad_mu!!
            e(set.ad.scale)
            e(set.ad.mu)
            if(do.mc_scale){
                if(do.mc_scale) mc_scale[.I + attr(mc_scale, 'prev_size'),, ] <- c(scale)
            }
        })),
    parent_mixins = adRoot,
    subtype = "HST")

