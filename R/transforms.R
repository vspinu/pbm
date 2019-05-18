
setClass("protoTransform",
         list(name = "character",
              TR = "function", 
              ITR = "function",
              LL = "function"))

protoTransform <- function(name = "identity",
                           TR = function(x) x,
                           ITR = function(y) y,
                           LL = function(y) 0){
    new("protoTransform", name = name, TR = TR, ITR = ITR, LL = LL)
}

tIdentity <- protoTransform()

## fixme: @LL should be changed into @ILL 
tLog <- protoTransform("tLog",
                       TR = log,
                       ITR = exp,
                       LL = function(y) y) ## log(exp(y))

tExp <- protoTransform("tExp",
                       TR = exp,
                       ITR = log,
                       LL = function(y) -log(y))

trRoot <- mixin(
    initFields = list(
        st_tr = tIdentity, 
        st = protoField(
            doc = "Call next method to retrive ST and then apply the transfrom.",
            function(arg){
                if(missing(arg))
                    st_tr@TR(nextProtoField("st")())
                else
                    nextProtoField("st")(st_tr@ITR(arg))
            }),
        mc_st = protoField(
            doc = "Call next method to retrive MC_ST and then apply the transfrom.",
            function(arg){
                st_tr@TR(nextProtoField("mc_st")(arg)) ## mc_st is readonly
            }),
        st = protoField(
            doc = "Call next method to retrive ST and then apply the transfrom.",
            function(arg){
                if(missing(arg))
                    st_tr@TR(nextProtoField("st")())
                else
                    nextProtoField("st")(st_tr@ITR(arg))
            }),
        mc_st = protoField(
            doc = "Call next method to retrive MC_ST and then apply the transfrom.",
            function(arg){
                st_tr@TR(nextProtoField("mc_st")(arg))
            })),
    initForms = list(
        init.R.build.ll_tr = form(
            if(!exists("ll_tr", inherits = F) || !identical(dim(ll), dim(ll_tr))){
                ll_tr <- array(0, dim = dim(ll), dimnames = dimnames(ll))
            }), 
        update.ll.ll_tr = form(
            if(do.ll_tr){
                ## needed for DIC but not for MHRW fixme: how to cache nicely?
                ## Both st method and ll_tr form apply the transformation. May
                ## be set.st.hook?
                ll_tr[] <- st_tr@LL(st_tr@TR(st))
            })))

trLog <- mixin(
    setFields = list(st_tr = tLog),
    parentMixins = trRoot, 
    subtype = "trLog")

trExp <- mixin(
    setFields = list(st_tr = tExp),
    parentMixins = trRoot, 
    subtype = "trExp")
        
        
        


