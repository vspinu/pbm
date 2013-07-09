### HC
PBM$initCells(defBC(type = "hc", prototype="*",
                    ## do.mc_st=FALSE,
                    ## do.st=FALSE,
                    setFields = list(
                        do.update = FALSE, ## fixme: doesn't take effect, so need to list all of them
                        do.st = FALSE, 
                        do.ll = FALSE,
                        do.mc_ll = FALSE,
                        do.mc_ll = FALSE,
                        do.pc_st = FALSE,
                        do.pc_ll = FALSE),
                    initForms = list(
                        init.C.build.mc_ll = NULL,
                        init.C.build.mc_st = NULL,
                        init.C.build._ll   = NULL,
                        init.M.validate = NULL,
                        init.R = NULL),
                    expr = expression({
                        ## .basic_foldable_objects <- c("st", "ll")
                        .basic_foldable_inxs <- c()
                        init.M.build <- init.M.build[c("st", "pos_in_C")]
                    })))

###_* DC
PBM$initCells(defBC(type = "dc", prototype = "*",
                    setFields = list(
                        do.mc_st = FALSE,
                        do.st = FALSE,
                        do.pc_st = FALSE,
                        do.pc_ll = TRUE),
                    expr = expression({
                        init.M.build[c("pos_in_C")] <- NULL
                    })))

###_* UC
PBM$initCells(defBC(type = "uc", prototype = "*"))

###_ + CONJ  (conjugate cells)
PBM$initCells(defBC(type = "conj", prototype = "uc",
                    initForms = list(
                        init.R.build_chITR = form(
                            chITR <- children[["i"]]$ITR),
                        init.M.validate.child_type = form())))
