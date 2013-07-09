
setClass("protoTransform",
         list(name = "character",
              TR = "function", 
              ITR = "function",
              LL = "function"))

protoTransform <- function(name = "identity",
                           TR = function(x) x,
                           ITR = function(y) y,
                           LL = function(y) 0){
    new("protoTransform", name = name,
        TR = TR, ITR = ITR, LL = LL)
}

tIdentity <- protoTransform()

tLog <- protoTransform("log",
                       TR = log,
                       ITR = exp,
                       LL = function(y) y) ## log(exp(y))

tExp <- protoTransform("exp",
                       TR = exp,
                       ITR = log,
                       LL = function(y) -log(y))

