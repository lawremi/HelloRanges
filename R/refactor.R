### =========================================================================
### refactoring utilities
### -------------------------------------------------------------------------
###

x <- quote({
    genome <- Seqinfo(genome = NA_character_)
    ga_a <- import(a, genome = genome)
    cov <- coverage(resize(granges(ga_a), 100L))
    ans <- GRanges(cov)
    ans <- subset(ans, score > 0)
    ans
})

extractFunction <- function(x, where=parent.frame()) {
    fun <- as.function(list(x))
    globals <- findGlobals(fun)
    argNames <- Filter(function(sym) !exists(sym, where), globals)
    formals(fun) <- setNames(rep(alist(x=), length(argNames)), argNames)
    fun
}
