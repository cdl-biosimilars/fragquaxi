test_that("mol_form constructor", {
  expect_equal(vec_size(new_molecular_formula()), 0)
  expect_equal(vec_size(new_molecular_formula(list(integer(0)))), 1)
  expect_equal(vec_size(new_molecular_formula(list(NA_integer_))), 1)
  expect_equal(vec_size(new_molecular_formula(rep(list(NA_integer_), 4))), 4)
  expect_equal(
    format(new_molecular_formula(list(
      NA_integer_,
      integer(0),
      rlang::set_names(integer(0), character(0)),
      c(H = 1L)
    ))), c(NA, "<empty>", "<empty>", "H")
  )
  expect_equal(vec_size(molecular_formula()), 0)
  expect_equal(
    new_molecular_formula(list(c(C = 0L, H = 1L), c(C = 1L, H = 2L), c(C = 2L, H = 3L))),
    molecular_formula(c("H", "C H2", "C2 H3"))
  )
  expect_equal(
    format(molecular_formula(c("", "H2 O", "C6 H12 O6", "H2SO4", NA))),
    c("<empty>", "H2 O", "C6 H12 O6", "H2 O4 S", NA)
  )
  expect_error(molecular_formula(c("", "H2 O", "C6 H12 O6", "H2SO4", "error")))
  expect_equal(molecular_formula("C6 H0"), molecular_formula("C6"))
  expect_equal(format(molecular_formula(c("", ""))), rep("<empty>", 2))
  expect_equal(format(molecular_formula("")), "<empty>")
  expect_true(is_molecular_formula(molecular_formula()))
})

test_that("mol_form parses spaces and element counts of 1", {
  expect_equal(
    molecular_formula("C35H48Cl1O10S1"),
    molecular_formula("C35 H48 Cl1 O10 S1")
  )
  expect_equal(
    molecular_formula("C35H48ClO10S"),
    molecular_formula("C35 H48 Cl1 O10 S1")
  )
})

test_that("duplicate elements in mol_form throw an error", {
  expect_error(molecular_formula("C6 C6"))
  expect_error(molecular_formula("C6 H12 C-3"))
})

test_that("mol_form appears in Hill order", {
  expect_equal(hill_order(c("C", "Cl")), c("C", "Cl"))
  expect_equal(hill_order(c("C", "H", "O")), c("C", "H", "O"))
  expect_equal(hill_order(c("H", "S", "O")), c("H", "O", "S"))
  expect_equal(format(molecular_formula("C35 Cl-3 O10 S1")), "C35 Cl-3 O10 S")
  expect_equal(format(molecular_formula("Cl4 C1")), "C Cl4")
  expect_equal(format(molecular_formula("H12 C6 O6")), "C6 H12 O6")
  expect_equal(format(molecular_formula("S1 O4 H2")), "H2 O4 S")
  expect_equal(format(molecular_formula("")), "<empty>")
})

test_that("mol_form coercing", {
  expect_error(vec_ptype2(integer(), molecular_formula()))
  expect_error(vec_ptype2(molecular_formula(), integer()))
  expect_equal(vec_ptype2(character(), molecular_formula()), character(0))
  expect_equal(vec_ptype2(molecular_formula(), character()), character(0))
  expect_equal(vec_ptype2(NA, molecular_formula()), molecular_formula())
  expect_equal(vec_ptype2(molecular_formula(), NA), molecular_formula())
})

test_that("mol_form casting", {
  s <- c("O2 Ti", "C6 H12 O6")
  m <- molecular_formula(s)
  s0 <- c("<empty>", "<empty>")
  m0 <- molecular_formula(s0)

  expect_equal(vec_cast(m, molecular_formula()), m)
  expect_equal(vec_cast(m, character()), s)
  expect_equal(vec_cast(molecular_formula(s0), character()), s0)
  expect_equal(vec_cast(s, molecular_formula()), m)
  expect_equal(vec_cast(s0, molecular_formula()), m0)
})

test_that("mol_form concatenation", {
  expect_equal(
    vec_c(molecular_formula(c("H2O", "H2 SO4")), molecular_formula(c("I4", NA)),
          .ptype = character()),
    c("H2 O", "H2 O4 S", "I4", NA)
  )
  expect_equal(
    vec_c(molecular_formula(c("H2O", "H2 SO4")), molecular_formula("I4")),
    molecular_formula(c("H2O", "H2 SO4", "I4"))
  )
  expect_equal(
    vec_c(molecular_formula(c("H2O", "C2 SO4")), molecular_formula("TiO2")),
    molecular_formula(c("H2O", "C2 SO4", "TiO2"))
  )
  expect_equal(
    vec_c("I4", molecular_formula(c("H2 SO4", "O3"))),
    c("I4", "H2 O4 S", "O3")
  )
})

test_that("mol_form comparisons", {
  s <- c(
    "C1 H2 O3",
    "C1 H2 O4",
    "C1 H3 O2",
    "C1 H3 O11",
    "C2 H1 O3",
    "C2 H1 O4",
    "C2 H2 O2",
    "C3 H2 O11",
    "H2 O3",
    "H2 O4",
    "H3 O1",
    "H3 O11",
    "O Cl1",
    "O Cl2",
    "O Cl3"
  )
  m <- molecular_formula(s)
  ms <- molecular_formula(sample(s))
  x <- molecular_formula("H2 O4")
  m_na <- molecular_formula(NA)

  expect_equal(sort(m), m)
  expect_equal(sort(ms), m)
  expect_equal(sum(m == x), 1)
  expect_true(max(ms) == molecular_formula("Cl3 O"))
  expect_true(min(ms) == molecular_formula("O3CH2"))
  expect_equal(sum(x < m), 5)
  expect_true(all((m <= x) == (x >= m)))
  expect_equal(x == m_na, NA)
  expect_equal(x < m_na, NA)
  expect_equal(x < NA, NA)
  expect_false(is.na(x))
  expect_true(is.na(m_na))
})

test_that("mol_form summation", {
  expect_equal(
    sum(molecular_formula(c("C2 H5 O", "C-2 H-5 O-1"))),
    molecular_formula("")
  )
  expect_equal(
    sum(molecular_formula(c("", "", "H"))),
    molecular_formula("H")
  )
  expect_equal(
    sum(molecular_formula(c("C2 H5 O", "Cl14", "H2 SO4"))),
    molecular_formula("C2 H7 Cl14 O5 S")
  )
  expect_equal(
    sum(molecular_formula(c("C2 H5 O", "Cl14", NA))),
    molecular_formula(NA)
  )
  expect_equal(
    sum(molecular_formula(c("C2 H5 O", "Cl14", NA)), na.rm = TRUE),
    molecular_formula("C2 H5 O Cl14")
  )
})

test_that("mol_forms add and subtract", {
  a <- molecular_formula("C6 H12")
  b <- molecular_formula("C7 O6")
  c <- c(a, b)
  d <- c(b, NA)
  e0 <- molecular_formula()
  e1 <- molecular_formula("")
  e2 <- molecular_formula(c("", ""))
  en <- molecular_formula(NA)
  expect_equal(a + b, molecular_formula("C13 H12 O6"))
  expect_equal(a - b, molecular_formula("C-1 H12 O-6"))
  expect_equal(a - a, e1)
  expect_equal(a + e0, e0)
  expect_equal(a - e0, e0)
  expect_equal(a + e1, a)
  expect_equal(a - e1, a)
  expect_equal(e1 - a, -a)
  expect_equal(e0 + e0, e0)
  expect_equal(e0 - e0, e0)
  expect_equal(e0 + e1, e0)
  expect_equal(e0 - e1, e0)
  expect_equal(e1 + e1, e1)
  expect_equal(e1 - e1, e1)
  expect_equal(-a, molecular_formula("C-6 H-12"))
  expect_equal(+a, molecular_formula("C6 H12"))
  expect_equal(a + a - a, a)
  expect_equal(a + en, en)
  expect_equal(a + c, c + a)
  expect_equal(c + d, molecular_formula(c("C13 H12 O6", NA)))
  expect_equal(e2 + e2, e2)
  expect_equal(c + e2, c)
  expect_equal(d + e2, d)
})

test_that("mol_forms multiply and divide", {
  a <- molecular_formula("C6 H-4 O")
  two_a <- molecular_formula("C12 H-8 O2")
  b <- molecular_formula(c("Cl4 I3", "C6 H-4 O"))
  two_b <- molecular_formula(c("Cl8 I6", "C12 H-8 O2"))
  d <- molecular_formula(c("C6 H-4 O", NA))
  two_d <- molecular_formula(c("C12 H-8 O2", NA))
  zero_d <- molecular_formula(c("", NA))
  e0 <- molecular_formula()
  e1 <- molecular_formula("")
  e2 <- molecular_formula(c("", ""))
  na1 <- molecular_formula(NA)
  na2 <- molecular_formula(c(NA, NA))

  expect_equal(2L * a, two_a)
  expect_equal(a * 2L, two_a)
  expect_equal(a * 1.6, a)
  expect_equal(a * 0, e1)
  expect_equal(0 * a, e1)
  expect_equal(a * NA_integer_, na1)
  expect_equal(NA_integer_ * a, na1)
  expect_equal(two_a / 2L, a)
  expect_equal(two_a / NA_integer_, na1)
  expect_error(2 / a)
  expect_equal(two_a %/% 2L, a)
  expect_error(2 %/% a)

  expect_equal(2L * b, two_b)
  expect_equal(b * 2L, two_b)
  expect_equal(b * 1.6, b)
  expect_equal(b * 0, e2)
  expect_equal(0 * b, e2)
  expect_equal(b * NA_integer_, na2)
  expect_equal(NA_integer_ * b, na2)
  expect_equal(two_b / 2L, b)
  expect_equal(two_b / NA_integer_, na2)
  expect_error(2 / b)
  expect_equal(two_b %/% 2L, b)
  expect_error(2 %/% b)

  expect_equal(2L * d, two_d)
  expect_equal(d * 2L, two_d)
  expect_equal(d * 1.6, d)
  expect_equal(d * 0, zero_d)
  expect_equal(0 * d, zero_d)
  expect_equal(d * NA_integer_, na2)
  expect_equal(NA_integer_ * d, na2)
  expect_equal(two_d / 2L, d)
  expect_equal(two_d / NA_integer_, na2)
  expect_error(2 / d)
  expect_equal(two_d %/% 2L, d)
  expect_error(2 %/% d)

  expect_equal(e0 * 2L, e0)
  expect_equal(e1 * 2L, e1)
  expect_equal(e2 * 2L, e2)
  expect_equal(e0 * 0L, e0)
  expect_equal(e1 * 0L, e1)
  expect_equal(e2 * 0L, e2)
  expect_equal(e0 * NA_integer_, e0)
  expect_equal(e1 * NA_integer_, na1)
  expect_equal(e2 * NA_integer_, na2)

  expect_equal(na1 * 2L, na1)
  expect_equal(na2 * 2L, na2)
  expect_equal(na1 * 0L, na1)
  expect_equal(na2 * 0L, na2)
  expect_equal(na1 * NA_integer_, na1)
  expect_equal(na2 * NA_integer_, na2)
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
