make_trait_ind_table <- function(alg_file, species) {

    alg <- ape::read.dna(file = alg_noempty, format = "fasta", as.matrix = TRUE)

    res <- array(0, dim = c(nrow(alg), length(species)))

    for (sp in species) {
        i <- grep(sp, dimnames(alg)[[1]])
        j <- match(sp, species)
        res[i, j] <- 1
    }

    test_sum <- apply(res, 1, sum)

    if (any(test_sum == 0)) {
        stop(paste(dimnames(alg)[[1]][which(test_sum == 0)], collapse = ", "),
             " sequence wasn't assigned any species.")
    }

    res <- rbind(c(species), res)
    res <- cbind(c("", dimnames(alg)[[1]]), res)

    trait_table <- paste0("data/traits/", basename(alg_file),
                          "_trait.table")
    write.table(res, file = trait_table, quote = FALSE,
                row.names = FALSE, col.names = FALSE, sep = ",")
    alg_noempty
}

make_trait_table <- function(alg_files, species) {
    sapply(alg_files, function(x) make_trait_ind_table(x, species))
}
