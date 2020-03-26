test_that("charge calculates correct m/z values", {
  expect_equal(
    charge(150000, 20:22),
    c("20" = 7501.008000, "21" = 7143.865143, "22" = 6819.189818)
  )
  expect_equal(
    charge(150000, 20:22, charge_agent_mass = 10),
    c("20" = 7510.000000, "21" = 7152.857143, "22" = 6828.181818)
  )
  expect_error(charge(150000, -20:-22))
})

test_that("molecular formulas are correctly calculated for FASTA files", {
  expect_equal(
    {
      mab_sequence <- system.file(
        "extdata", "mab_sequence.fasta",
        package = "fragquaxi"
      )
      load_protein_sequence(mab_sequence)
    },
    molecular_formula("C6464 H9982 N1706 O2014 S44")
  )
  expect_equal(
    {
      mab_sequence <- system.file(
        "extdata", "mab_sequence.fasta",
        package = "fragquaxi"
      )
      load_protein_sequence(mab_sequence, disulfides = 16)
    },
    molecular_formula("C6464 H9950 N1706 O2014 S44")
  )
})
