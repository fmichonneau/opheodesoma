build_haplotype_network <- function(alg_file, genus_name = "Opheodesoma",
                                    unique_name = "ESU1|ESU2",
                                    palette) {
    alg <- ape::read.dna(alg_file, format = "fasta")
    alg <- alg[grepl(unique_name, dimnames(alg)[[1]]), ]
    hap <- haplotype(alg)
    net <- haploNet(hap)
    pie <- stack(setNames(attr(hap, "index"), rownames(hap)))
    reg_expr <- paste0("^(", genus_name, ")_(", unique_name, ")", "_(.+)")
    pie$sp <- gsub(reg_expr, "\\1_\\2", dimnames(alg)[[1]])
    tab <- table(pie$ind, pie$sp)
    cols <- wesanderson::wes_palette(palette, min(5, ncol(tab)))

    if (dim(net)[1] < 2) return(NULL)

    plot(net, size = sapply(attr(hap, "index"), length),  pie = tab, labels = FALSE,
         bg = cols, threshold = 0, fast = FALSE, scale.ratio = 2)
    legend("bottomleft", legend = colnames(tab), fill = cols)
}

all_networks <- function(genera = c("Euapta", "Opheodesoma"),
                         species = c("tahitiensis|godeffroyi|lappa",
                                     "ESU1|ESU2|giantRed|redsea")
                         ) {
    all_alg <- list.files(path = "data/alignments", pattern = "-trimmed.afa$")

    for (i in seq_along(genera)) {
        sub_alg <- grep(tolower(genera[i]), all_alg, value = TRUE)
        sapply(sub_alg, function(x) {
            message(x, appendLF = FALSE)
            if (grepl("euapta.+LSU|COI", x)) {message ("... skipped ..."); return(NULL)}
            svg(file = file.path("networks", gsub("-trimmed.afa", "-network.svg", x)),
                width = 10)
            on.exit(dev.off())
            build_haplotype_network(file.path("data", "alignments", x),
                                    genus_name = genera[i],
                                    unique_name = species[i],
                                    palette = if (genera[i] == "Euapta") "Royal2" else "Rushmore")
            message(" ... DONE.")
        })
    }
}
