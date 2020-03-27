test_that("mol_form correctly parses spaces and element counts of 1", {
  expect_equal(
    molecular_formula("C35H48Cl1O10S1"),
    molecular_formula("C35 H48 Cl1 O10 S1")
  )
  expect_equal(
    molecular_formula("C35H48ClO10S"),
    molecular_formula("C35 H48 Cl1 O10 S1")
  )
})

test_that("mol_form is correctly converted to character in Hill order", {
  expect_equal(format(molecular_formula("C35 Cl-3 O10 S1")), "C35 Cl-3 O10 S")
  expect_equal(format(molecular_formula("Cl4 C1")), "C Cl4")
  expect_equal(format(molecular_formula("H12 C6 O6")), "C6 H12 O6")
  expect_equal(format(molecular_formula("S1 O4 H2")), "H2 O4 S")
  expect_equal(format(molecular_formula("")), "empty")
})

test_that("mol_form arithmetic is correct", {
  expect_equal(
    molecular_formula("C6 H4 O6") + molecular_formula("C3 H12 S2"),
    molecular_formula("C9 H16 O6 S2")
  )
  expect_equal(
    molecular_formula("C6 H4 O6") - molecular_formula("C3 H12 S2"),
    molecular_formula("C3 H-8 O6 S-2")
  )
  expect_equal(
    molecular_formula("C6 H4 O6") - molecular_formula("C6 H4 O6"),
    molecular_formula("")
  )
  expect_equal(
    molecular_formula("") - molecular_formula("C6 H4 O6"),
    molecular_formula("C-6 H-4 O-6")
  )
  expect_equal(
    -molecular_formula("C6 H4 O6"),
    molecular_formula("C-6 H-4 O-6")
  )
  expect_equal(
    +molecular_formula("C6 H4 O6"),
    molecular_formula("C6 H4 O6")
  )
  expect_equal(
    2 * molecular_formula("C6 H-4 O"),
    molecular_formula("C12 H-8 O2")
  )
  expect_equal(
    molecular_formula("C6 H-4 O") * 2,
    molecular_formula("C12 H-8 O2")
  )
})

test_that("duplicate elements in mol_form throw an error", {
  expect_error(molecular_formula("C6 C6"))
  expect_error(molecular_formula("C6 H12 C-3"))
})

test_that("mass calculations are correct", {
  expect_equal(get_mass(molecular_formula("C6 H12 O6")), 180.156)
  expect_equal(get_mass(molecular_formula("C6 H12 O6"), "monoisotopic"), 186.06339)
  expect_true(is.na(get_mass(molecular_formula("C6 H12 Rd"))))
})

test_that("new mass sets are defined correctly", {
  expect_equal(
    {
      my_masses <- new_mass_set(c(H = 1, C = 10, O = 100))
      get_mass(molecular_formula("C6 H9 O5"), my_masses)
    },
    569
  )
  expect_true(
    is.na({
      my_masses <- new_mass_set(c(H = 1, C = 10, O = 100))
      get_mass(molecular_formula("C6 H9 S5"), my_masses)
    })
  )
  expect_equal(
    {
      my_masses <- new_mass_set(
        c(H = 1, C = 10, O = 100),
        inherits_from = "average"
      )
      get_mass(molecular_formula("C6 H9 S5"), my_masses)
    },
    229.3
  )
})
