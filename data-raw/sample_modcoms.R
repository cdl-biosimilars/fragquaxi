sample_modcoms <- tibble::tribble(
  ~modcom_name, ~Hex, ~HexNAc, ~Fuc,
  "G0F/G0",        6,       8,    1,
  "G0F/G0F",       6,       8,    2,
  "G0F/G1F",       7,       8,    2,
  "G1F/G1F",       8,       8,    2,
  "G1F/G2F",       9,       8,    2,
  "G2F/G2F",      10,       8,    2,
)

usethis::use_data(sample_modcoms, overwrite = TRUE)
