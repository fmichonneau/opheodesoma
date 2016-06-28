make_raxml_trees <- function() {
    lst_alg <- list.files(path = "data/alignments",
                          pattern = "-noEmpty.afa$",
                          full.names = TRUE)
    lapply(lst_alg, raxml_tree)
}


raxml_tree <- function(alg_file, raxml_path = "data/./raxmlHPC-PTHREADS-SSE3") {
    phy_file <- fas2phy(alg_file, overwrite = TRUE)
    phy_file <- names(phy_file)

    raxml_output <- "data/raxml"
    raxml_suffix <- gsub("-noEmpty.phy", "", basename(phy_file))

    if (!file.exists(raxml_output)) {
        dir.create(raxml_output)
    } else {
        lst_raxml_suffix <- list.files(path = raxml_output,
                                       pattern = paste0(raxml_suffix, "$"),
                                       full.names = TRUE)
        file.remove(lst_raxml_suffix)
    }

    raxml_phy <- file.path(raxml_output, basename(phy_file))
    file.rename(phy_file, raxml_phy)

    raxml_cmd <- paste(raxml_path, "-s", raxml_phy, "-m GTRGAMMA",
                      "-f a -p 10101 -x 10101 -# 500 -n ", raxml_suffix,
                      "-T8 -w", file.path(getwd(), raxml_output))
    system(raxml_cmd)
}
