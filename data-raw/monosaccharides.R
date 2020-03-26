monosaccharides <-
  readr::read_csv("data-raw/monosaccharides.csv")

usethis::use_data(monosaccharides, overwrite = TRUE)
