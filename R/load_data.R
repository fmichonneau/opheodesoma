load_csv <- function(f) {
    read.csv(file = f, stringsAsFactors = FALSE)
}

### Extractions ----------------------------------------------------------------

get_extractions <- function(db, .genus) {
     dplyr::filter(db, genus == .genus) %>%
        dplyr::select(extraction_number) %>%
        .[[1]] %>%
        strsplit(split = ",") %>%
         unlist
}

euapta_extractions <- function(db, genus = "Euapta") {
   get_extractions(db, genus)
}

opheodesoma_extractions <- function(db, genus = "Opheodesoma") {
    get_extractions(db, genus)
}



### Utils ----------------------------------------------------------------------

labels_from_extractions <- function(lst) {
    dt <- load_csv("data/opheodesoma_euapta_data.csv")
    wch <- vapply(lst, function(x) grep(x, dt[["extraction_number"]]),
                  integer(1))
    paste(dt[["genus"]][wch], dt[["species"]][wch],
          dt[["catalogue_number"]][wch], sep = "_")
}
