bedtools_merge <- function(cmd = "--help") {
    do_R_call(R_bedtools_merge, BEDTOOLS_MERGE_DOC, cmd)
}

COLNAMES_FOR_FORMAT <- list(
    bed = c("seqnames", "start", "end", "name", "score", "strand",
            "start(thick)", "end(thick)", "itemRgb",
            "lengths(blocks)", "width(blocks)", "start(blocks)"),
    ## we make no attempt to handle the last column
    gff = c("seqnames", "source", "type", "start", "end", "score", "strand",
            "phase"),
    ## this subscripts into granges(vcf)
    vcf = c("seqnames", "start", "names", "REF", "ALT", "QUAL"),
    bam = c("names", NA, "seqnames", "start", "mapq", "cigar",
            "mrnm", "mpos", "isize", "seq", "qual")
)

EXPRS_FOR_OPS <- list(
    sum=quote(sum(X)),
    min=quote(min(X)),
    max=quote(max(X)),
    absmin=quote(min(abs(X))),
    absmax=quote(max(abs(X))),
    mean=quote(mean(X)),
    median=quote(median(X)),
    mode=quote(distmode(X)),
    antimode=quote(distmode(X, anti=TRUE)),
    collapse=quote(unstrsplit(X, delim)),
    distinct=quote(unstrsplit(unique(X), delim)),
    count=quote(lengths(X)),
    count_distinct=quote(lengths(unique(X))),
    sstdev=quote(sd(X)),
    freq=quote(table(X)),
    first=quote(drop(phead(X, 1L))),
    last=quote(drop(ptail(X, 1L)))
)

normCandO <- function(i, c, o, delim = ",") {
    cols <- as.integer(if (is.character(c))
                           strsplit(c, ",", fixed=TRUE)[[1L]]
                       else c)
    if (any(is.na(cols) | cols <= 0L)) {
        stop("'c' must be a comma-separated list of positive integers")
    }
    ops <- strsplit(o, ",", fixed=TRUE)[[1L]]
    len <- max(length(cols), length(ops))
    cols <- recycleIntegerArg(cols, "c", len)
    ops <- recycleCharacterArg(ops, "o", len)
    cn <- COLNAMES_FOR_FORMAT[[getFormat(i)]][cols]
    if (any(is.na(cn))) {
        stop("unsupported columns: ", paste(cols[is.na(cn)], collapse=", "))
    }
    expr <- EXPRS_FOR_OPS[ops]
    badops <- vapply(expr, is.null, logical(1L))
    if (any(badops)) {
        stop("unsupported ops: ", paste(ops[badops], collapse=", "))
    }
    op_exprs <- mapply(function(cni, expri) {
        X <- parse(text=cni)[[1L]] # cni might be an unparsed call
        eval(call("substitute", expri, list(X=X, delim=delim)))
    }, cn, expr, SIMPLIFY=FALSE)
    names(op_exprs) <- paste0(cn, ".", ops)
    list(cn=cn, exprs=op_exprs)
}

.aggregateCall <- function(.gr_i, grouping, exprs, drop=TRUE) {
    .agg <- as.call(c(list(quote(aggregate), .gr_i), grouping, exprs))
    if (!drop)
        .agg$drop <- FALSE
    .agg
}

R_bedtools_merge <- function(i, s=FALSE, S=c("any", "+", "-"),
                             d=0L, c=NULL, o="sum", delim = ",")
{
    S <- match.arg(S)
    stopifnot(isSingleString(i) || hasRanges(i),
              isTRUEorFALSE(s),
              isSingleNumber(d), d >= 0L,
              is.null(c) || is.numeric(c) || isSingleString(c),
              is.null(o) || isSingleString(o),
              isSingleString(delim))

    if (!is.null(c)) {
        co <- normCandO(i, c, o, delim)
        cn <- co$cn
    } else {
        cn <- NULL
    }

    R(genome <- Seqinfo(genome=NA_character_))

    .gr_i <- importA(i, cn)
    if (isVcf(i)) {
        .gr_i <- .R(granges(expand(.gr_i)))
    }
    
    .gr_i_o <- prepOverlapRanges(i, FALSE)

    ignore.strand <- !s && S == "any"

    if (S != "any") {
        .gr_i_o <- .R(subset(.gr_i_o, strand == S))
    }

    .reduce <- .R(reduce(.gr_i_o, ignore.strand=ignore.strand))

    if (!is.null(c)) {
        .reduce$with.revmap <- TRUE
    }

    if (d > 0L) {
        .reduce$min.gapwidth <- d + 1L
    }
    
    R(ans <- .reduce)
    
### NOTE: the grouping information is preserved, so this should cover
### bedtools_cluster.
    if (!is.null(c)) {
        .agg <- .aggregateCall(.gr_i, quote(mcols(ans)$revmap), co$exprs,
                               drop=FALSE)
        R(mcols(ans) <- .agg)
    }
    
    R(ans)
}

BEDTOOLS_MERGE_DOC <-
    "Usage:
       bedtools_merge [options]
     Options:
       -i <FILE>  BAM/BED/GFF/VCF file.
       -s  Force strandedness. That is, only merge features that are the
           same strand. By default, this is disabled.
       -S <STRAND>  Force merge for one specific strand only.
           Follow with + or - to force
           merge from only the forward or reverse strand, respectively.
           By default, merging is done without respect to strand.
       -d <DIST>  Maximum distance between features allowed for features
           to be merged. Default is 0. That is, overlapping and/or book-ended
           features are merged [default: 0].
       -c <COL1,...>  Specify columns from the input file to operate upon
           (see -o option, below).
           Multiple columns can be specified in a comma-delimited list.
       -o <OP1,...>  Specify the operation that should be applied to -c.
           Valid operations:
             sum, min, max, absmin, absmax,
             mean, median,
             collapse (i.e., print a delimited list (duplicates allowed)),
             distinct (i.e., print a delimited list (NO duplicates allowed)),
             count
             count_distinct (i.e., a count of the unique values in the column),
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

do_bedtools_merge <- make_do(R_bedtools_merge)
