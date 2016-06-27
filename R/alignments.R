### Alignments -----------------------------------------------------------------

get_alignment <- function(lst, target, markers) {
    stopifnot(identical(length(markers), length(target)))
    res <- chopper::mergeSeq(lst, output = tempdir(),
                             markers = markers)
    f <- attr(res, "aligned_files")

    res_files <- character(length(f))
    for (i in seq_along(f)) {
        alg <- ape::read.dna(file = f[i], format = "fasta", as.matrix = TRUE)
        lbl <- labels_from_extractions(dimnames(alg)[[1]])
        dimnames(alg)[[1]] <- make.unique(lbl)
        alg <- chopper::cleanSeqLabels(alg)
        ape::write.dna(alg, file = target[i], format = "fasta", colw = 1000, colsep = "")
        res_files[i] <- target[i]
    }
    res_files
}

get_all_loci <- function(lst, genus) {
    markers <- c("COI", "H3a", "ITS", "LSU", "18S")
    targets <- paste0("data/alignments/", genus, "-", markers, ".afa")
    get_alignment(lst, targets, markers)
}


### Converts to NEXUS ---------------------------------------------------------

all_loci_nexus <- function(lst_files) {
    lf <- lst_files
    res <- sapply(lf, chopper::alg2nex, format = "fasta")
    gsub("\\.[a-z]{2,}", "\\.nex", lf, ignore.case = TRUE)
}

### Remove empty sequences -----------------------------------------------------

remove_empty_sequences <- function(alg_file) {
    alg_noempty <- gsub("\\.afa$", "-noEmpty.afa", alg_file)
    chopper::removeEmptySeqs(file = alg_file, output = alg_noempty)
}

make_no_empty_sequences <- function(alg_files) {
    sapply(alg_files, remove_empty_sequences)
}



trim_alignment <- function(alg_file) {
    alg <- read.dna(file = alg_file, format = "fasta", as.character = TRUE, as.matrix = TRUE)
    is_gap <- apply(alg, 2, function(x) sum(x == "-"))
    first_no_gap <- min(which(is_gap == 0))
    last_no_gap <- length(is_gap) - min(which(rev(is_gap) == 0)) + 1
    ape::write.dna(alg[, first_no_gap:last_no_gap], file = gsub("noEmpty", "trimmed", alg_file),
                   format = "fasta", colsep = "", colw = 10000)
}

make_trim_alignment <- function() {
    alg_files <- list.files(path = "data/alignments", pattern = "noEmpty.afa",
                            full.names = TRUE)
    sapply(alg_files, trim_alignment)
}
