## TRANFROMER BCELLS

PBM$initCells(defBC(type = "tr", prototype = "*",
                    ## API: NOTE: tr cell, passes likelihood computation to next
                    ## cell. Examples: .update.ll, .set.ll_is_old, .get.ll

                    ## API: LL and ST are one thing in TR cells. All external ll
                    ## updating also involve st updating. See below.
                    
                    setFields = list(
                        ## tr cells are updated by parent cells in UPDATE.tr_children
                        do.update = FALSE), 
                    initFields = list(
                        tr = protoField(.field.tr,
                            doc = "API: this is an object of protoTransform class describing the state and likelihood transformation.")),
                    initForms = list(
                        init.R.build.st = form(e(init.M.build.st))), 
                    setForms = list(
                        set.st = form(st[] <- do.call(tr@TR, lapply(seq_along(parents), pst))),
                        set.ll = form(ll[] <- tr@LL(st)),
                        set.rand.st = form(e(set.st)),
                        init.M.build.st = form({
                            ## preclude the transform from droping dimensions 
                            .fields$st <- pst(1)
                            st[] <- tr@TR(st)
                        })),
                    setMethods = list(
                        .get.ll =
                        function(ptype, token = runif(1L), exclude_tr = TRUE){
                            if(exclude_tr)
                                .get.ll(children, .type, toke = token, exclude_tr = exclude_tr)
                            else
                                ## return a matrix of this ll and parents LLs
                                structure(cbind(rowsum.default(c(ll), pix(, ptype)),
                                                .get.ll(children, .type, token = token, exclude_tr = exclude_tr)),
                                          dimnames = list(NULL, c(.type, names(children))))
                        },
                        ## API: all external api calls that manipulate ll, also manipulate st:
                        .resync.ll = 
                        function(rejects, ptype, token = runif(1L)){
                            if(.resync_object(.self, rejects, "ll", "_ll", ptype, token))
                                .resync_object(.self, rejects, "st", "_st", ptype, token, force = TRUE)
                            .resync.ll(children, rejects, .type, token)
                        },
                        .update.ll = function() {
                            e(update.st) # does nothing if already set
                            e(update.ll)
                            .update.ll(children)
                        },
                        .set.ll_is_old = function(){
                            e(set.st_is_old)
                            e(set.ll_is_old)
                            .set.ll_is_old(children)
                        } ,
                        .set.ll_is_updated = function(){
                            e(set.st_is_updated)
                            e(set.ll_is_updated)
                            .set.ll_is_updated(children)
                        }), 
                    expr = expression(
                        tr <- tIdentity)))

