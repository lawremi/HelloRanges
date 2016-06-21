bedtools_groupby <- function(cmd = "--help") {
    do_R_call(R_bedtools_groupby, BEDTOOLS_GROUPBY_DOC, cmd)
}

normG <- function(i, g) {
    groups <- as.integer(if (is.character(g)) strsplit(g, ",", fixed=TRUE)[[1L]]
                         else g)
    cn <- COLNAMES_FOR_FORMAT[[getFormat(i)]][groups]
    as.formula(paste("~", paste(cn, collapse="+")))
}

R_bedtools_groupby <- function(i, g = 1:3, c, o="sum", delim=",")
{
    stopifnot(isSingleString(i) || hasRanges(i),
              isSingleString(g) || (is.numeric(g) && !anyNA(g) && all(g > 0L)),
              isSingleString(c) || (is.numeric(c) && !anyNA(c) && all(c > 0L)),
              isSingleString(o),
              isSingleString(delim))

    co <- normCandO(i, c, o, delim)

    R(genome <- Seqinfo(genome=NA_character_))

    .gr_i <- importA(i, co$cn)
    if (isVcf(i)) {
        .gr_i <- .R(granges(expand(.gr_i)))
    }

    formula <- normG(i, g)
    vars <- all.vars(formula)
    if (setequal(setdiff(vars, "strand"), c("seqnames", "start", "end"))) {
        formula <- NULL
        if (!"strand" %in% vars) {
            .gr_i <- .R(unstrand(.gr_i))
        }
    }
    .agg <- .aggregateCall(.gr_i, formula, co$exprs)
    
    R(ans <- .agg)

    R(ans)
}

BEDTOOLS_GROUPBY_DOC <-
    "Usage:
       bedtools_groupby [options]
     Options:
       -i <FILE>  BAM/BED/GFF/VCF file.
       -g <COL1,...>  Specifies which column(s) (1-based) should be used to
                      group the input. Columns may be comma-separated with each
                      column must be explicitly listed. Or, ranges (e.g. 1-4)
                      are also allowed. [default: 1,2,3]
       -c <COL1,...>  Specify columns from the input file to operate upon
           (see -o option, below).
           Multiple columns can be specified in a comma-delimited list.
       -o <OP1,...>  Specify the operation that should be applied to -c.
           Valid operations:
             sum, min, max, absmin, absmax,
             mean, median, mode, antimode, sstdev,
             collapse (i.e., print a delimited list (duplicates allowed)),
             distinct (i.e., print a delimited list (NO duplicates allowed)),
             count
             count_distinct (i.e., a count of the unique values in the column),
             freq (table), first, last
           Default: sum
           Multiple operations can be specified in a comma-delimited list.
           If there is only column, but multiple operations, all ops will be
           applied on that column. Likewise, if there is only one operation, but
           multiple columns, that operation will be applied to all columns.
           Otherwise, the number of columns must match the the number of ops,
           and will be applied in respective order.

           E.g., -c 5,4,6 -o sum,mean,count will give the sum of column 5,
           the mean of column 4, and the count of column 6.
           The order of output cols matches the ordering given in the command.
       --delim <DELIM>  Specify a custom delimiter for the collapse operation
               Example: -delim \"|\"
               [default: ,]"

do_bedtools_groupby <- make_do(R_bedtools_groupby)
