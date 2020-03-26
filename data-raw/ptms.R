ptms <-
  readr::read_csv("data-raw/ptmlist.csv") %>%
  dplyr::select(abbreviation, name = description, formula)

usethis::use_data(ptms, overwrite = TRUE)
