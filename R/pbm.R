setContextClass("PBM", cellClass="BC")
setMethod("initialize", "PBM",
          function(.Object, type = "--", prototype = NULL, ...){
              if(!is.null(prototype)){
                  cur_mir <- prototype$mirror
                  ## inherit .cells from base not from mirrors
                  prototype$mirror <- "base"
                  on.exit(prototype$mirror <- cur_mir)
              }
              pbm <- callNextMethod(.Object, type = type, prototype = prototype, ...)

              pbm$evalq({
                  current_mirror <- "base"
                  folds <- NULL
                  folds_names <- NULL
                  .cells <- as(.cells, "bCellContainer")
                  mirrors <- list(base = .cells)
              })
              pbm
          })

## don't delete
## setAs("PBM","graphNEL", .as_PBM_graphNEL)
## str(as(M, "graphNEL"))

pbm <- function(type = "--", rootParentEnv = NULL,
                initCells = list(), cells = list(), ...){
    initCells <- c(initCells, list(...))
    out <- new("PBM", type = type, rootParentEnv = rootParentEnv,
               initCells = initCells, cells = cells,
               expr = expression({}))
}

###_* PBM context
PBM <- getClassDef("PBM")@defaultContext

PBM$initFields(mirror = protoField(.field.current_mirror),
               mirrors = protoReadOnlyField("mirrors"),
               ## mirrors_train = protoReadOnlyField("mirrors_train"),
               ## mirrors_test = protoReadOnlyField("mirrors_test"),
               folds_names = protoReadOnlyField("folds_names"))

PBM$initMethods(getLL =
                function(sum = FALSE, mc = FALSE){
                    if(mc){
                        LL <- map_over_model_cells(.cells, function(cell){
                            if(cell$do.mc_ll) rowSums(as.matrix(cell$mc_ll))
                        })
                        if(sum) rowSums(do.call(cbind, LL)) else LL
                    }else{
                        LL <- unlist(map_over_model_cells(.cells, function(cell){
                            if(cell$do.ll) sum(cell$ll)
                        }))
                        if(sum) sum(LL) else LL
                    }
                }, 
                resetPC =
                function(){
                    map_over_model_cells(.cells, function(cell){
                        assign("._pc_done", FALSE, cell)
                    })},
                initTestMirror =
                function(test_folds, train_folds = NULL, switch = T,
                         test_mname = NULL, train_mname = NULL){
                    folds <- get("folds", .self, inherits = F)
                    if(is.numeric(test_folds))
                        test_folds <- folds[test_folds]
                    if(is.null(train_folds))
                        train_folds <- setdiff(folds, test_folds)
                    if(!all(what <- test_folds %in% folds))
                        stop("The following test folds are not declared: ", paste(test_folds[!what], sep = ", "))
                    if(!all(what <- train_folds %in% folds))
                        stop("The following train folds are not declared: ", paste(train_folds[!what], sep = ", "))
                    if(is.null(test_mname))
                        test_mname <- paste(test_folds, collapse = ":")
                    if(is.null(train_mname))
                        train_mname <- paste(train_folds, collapse = ":")
                    .create_folded_mirror(.self, train_folds, train_mname)
                    .create_folded_predict_mirror(.self, test_folds, test_mname, train_mname)
                    ## assign("mirrors_train", union(.self[["mirrors_train"]], test_mname), .self)
                    ## assign("mirrors_test", union(.self[["mirrors_test"]], test_mname), .self)
                    if(switch)
                        .self$mirror <- test_mname
                    invisible(NULL)})
