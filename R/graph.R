setAs("bCellContainer", "graphNEL",
      function(from){
          .gr <- as(as(from, "cellContainer"), "graphNEL")
          leaf_names <- leafNames(from)
          ## nodeData(.gr, leaf_names, "type") <- "model" ## change leaf to model
          for(nm in leaf_names){
              bc <- get(nm, from, inherits = F)
              bcname <- .getType(bc)
              if(length(parents <- bc$parents)){
                  parents_names <- sapply(parents, .getType)
                  .gr <- addEdge(parents_names, bcname, .gr)
                  edgeData(.gr, parents_names, bcname, "type") <- "leaf"
              }
          }
          .gr
      })

setMethod("plot", c("bCellContainer", "ANY"), 
          function(x, y, types = "leaf", ...){
              callNextMethod(x, types = types, ...)
          })
