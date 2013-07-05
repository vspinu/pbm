###_* UPDATE

setMethod("update", "PBM",
          function(object, nr_iter=10, thin=1, append=T, reinit = F, do.pr_bar=TRUE, ...){
              . <- as.environment(object)
              library(abind)
              ## get current mirror
              .curcells <- get(".cells", envir = object)
              is_test <- .curcells@is_mirror && length(.curcells@train_mirror_name)
              ## change the mirror to train if the current is a test mirror
              .cells <-
                  if(is_test) .$mirrors[[.curcells@train_mirror_name]]
                  else .curcells
              root_env <- as.environment(get("*", .cells))
              model_cells <- .get_model_cells(.cells)
              to_iterate <- sapply(model_cells, function(bc){
                  bc$do.update
              })
              .UP <- model_cells[rev(which(to_iterate))]
              .cells@sims$thin <- as.integer(thin)
              nr_iter <- as.integer(nr_iter)
              if(append ||  .cells@sims$end == 0L){
              }else{
                  .strip_mcs(.UP)  ## define the func <- todo:
                  .cells@sims$start <- .cells@sims$end+1L
              }
              assign(".nr_iter", nr_iter, envir = root_env) ## nee for build.mc_st
              assign(".thin", as.integer(thin), envir = root_env)
              if(reinit || .cells@sims$end==0L){
                  root_env$.N <- 0L
                  .$.runEnds <- integer()
                  .cells@sims$end <- 0L
                  .cells@sims$start <- 1L
                  .initialize_M(model_cells)
                  ## .strip_mcs(.UP) init.M.build.mc_st and mc_ll do the job
              }
              .cells@sims$end <- .cells@sims$end + nr_iter ## available for R.builders!!
              ## so ugly
              if(is_test) .$mirrors[[.curcells@train_mirror_name]] <- .cells
              else .$.cells <- .cells
              .initialize_R(model_cells)
              iter_along <- seq_len(nr_iter)
              nr_BCs <- length(.$.MCs)
              iter_along_thin <- seq_len(thin-1)
              if (do.pr_bar){
                  .pb <- txtProgressBar(min = iter_along[[1L]],
                                        max = tail(iter_along, 1),
                                        initial = 1L, style = 3, width = 50,
                                        char = "*")
                  on.exit(close(.pb))
              }
              .N <- get(".N", root_env, inherits = F)
              itime <- proc.time()[[3]]
              for(.I in iter_along){
                  .Internal(assign(".I", .I, root_env, FALSE))
                  .Internal(assign(".T", 0L, root_env, FALSE))
                  .N <- .N + 1L
                  .Internal(assign(".N", .N, root_env, FALSE))
                  for(bc in .UP)
                      evalq(e(update), envir=bc)
                  for(.T in iter_along_thin){
                      .Internal(assign(".T", .T, root_env, FALSE))
                      .N <- .N + 1L
                      .Internal(assign(".N", .N, root_env, FALSE))
                      for(bc in .UP) evalq(e(update), envir=bc)
                  }
                  if(do.pr_bar) setTxtProgressBar(pb = .pb, value = .I)
              }
              ftime <- proc.time()[[3]]
              cat("\nTotal Time Elapsed: ",
                  round((total.time <- ftime-itime)/60, 2), " mins ", fill=TRUE)
              if(root_env$do.timers && !root_env$do.debug){
                  cat("\nPartial Times: \n")
                  tdata <- as.data.frame(do.call(rbind,
                                                 lapply(.UP, function(bc){
                                                     c(minutes=round((bc.t <- get(".t", envir = bc))/60, 3),
                                                       `%total`=round(bc.t/total.time*100, 1))
                                                 })))
                  tdata <- tdata[order(tdata$minutes, decreasing=TRUE), ]
                  print(do.call(rbind, list(tdata, TOTAL=colSums(tdata))))
              }
          })

.initialize_M <- function(mCells, force = TRUE){
    if(is(mCells, "BC"))
        mCells <- list(mCells)
    lapply(mCells, function(bc){
        if(force || !get("is.initialised.M", bc))
            eval(bc$init.M, envir=bc)
        bc$evalq(is.initialised.M <- TRUE)
    })
}

.initialize_R <- function(mCells){
    if(is(mCells, "BC"))
        mCells <- list(mCells)
    lapply(mCells, function(bc){
        eval(bc$init.R, envir=bc)
        bc$evalq(is.initialised.R <- TRUE)
        ## if(bc$do.expand.ue_UPDATE && bc$do.update && exists("ue_UPDATE", envir=bc))
        ##   bc[["ue_UPDATE"]] <- expand.e(bc$ue_UPDATE, bc)
    })
}

