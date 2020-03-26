amino_acids <-
  readr::read_csv("data-raw/amino_acids.csv")

usethis::use_data(amino_acids, overwrite = TRUE)
