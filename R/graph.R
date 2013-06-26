setAs("bCellContainer", "graphNEL",
      function(from){
          .gr <- as(as(from, "cellContainer"), "graphNEL")
          leaf_names <- leafNames(from)
          nodeData(.gr,leaf_names, "type") <- "model" ## change leaf to model
          for(nm in leaf_names){
              bc <- get(nm, from, inherits = F)
              if(length(parents <- bc$parents)){
                  parents_names <- sapply(parents, .getType)
                  .gr <- addEdge(parents_names,nm, .gr)
                  edgeData(.gr, parents_names, nm, "type") <- "model"
              }
          }
          .gr
      })

setMethod("plot", c("bCellContainer", "ANY"), 
          function(x, y, types = "model", ...){
              callNextMethod(x, types = types, ...)
          })
