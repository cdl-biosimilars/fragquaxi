library(tidyverse)

average_masses <-
  read_csv("data-raw/average_mass.csv", col_types = "cciccc") %>%
  mutate(
    average = case_when(
      is.na(conventional_atomic_weight) ~ abridged_standard_atomic_weight,
      TRUE ~ conventional_atomic_weight
    ) %>%
      str_replace_all("[\\(\\)–]", "") %>%
      parse_double()
  ) %>%
  select(element = symbol, name = element_name, Z = atomic_number, average)

all_mono_masses <-
  read_fwf(
    "data-raw/mass16.txt",
    fwf_cols(
      cc = 1, NZ = 3, N = 5, Z = 5, A = 5, space1 = 1, el = 3, o = 4,
      space2 = 1, mass = 13, unc_mass = 11, binding = 11, unc_binding = 9,
      space3 = 1, B = 2, beta = 11, unc_beta = 9, space4 = 5,
      atomic_mass = 12, unc_atomic_mass = 11
    ),
    col_types = cols(A = "i"),
    skip = 39
  ) %>%
  select(element = el, A, atomic_mass)

isotope_abundances <-
  read_csv("data-raw/iso_comp16.csv") %>%
  fill(E) %>%
  select(element = E, A, best_measurement) %>%
  mutate(
    A = parse_integer(A, na = "–"),
    abundance = case_when(
      best_measurement %in% c("1", "–") ~ 1,
      TRUE ~ best_measurement %>%
        str_extract("[0-9 \\.]*") %>%
        str_replace_all(" ", "") %>%
        parse_double()
    )
  ) %>%
  group_by(element) %>%
  filter(abundance == max(abundance))

mono_masses <-
  all_mono_masses %>%
  semi_join(isotope_abundances, by = c("element", "A")) %>%
  mutate(monoisotopic = A + parse_double(atomic_mass) / 1e6)

atomic_masses <-
  average_masses %>%
  left_join(mono_masses, by = "element") %>%
  select(-atomic_mass)

atomic_masses

usethis::use_data(atomic_masses, overwrite = TRUE)
