#' Sort a molecular formula according to the Hill system.
#'
#' * If the formula contains carbon and hydrogen,
#'   list the number of carbon atoms,
#'   followed by the number of hydrogen atoms,
#'   followed by the remaining elements in alphabetical order.
#' * If the formula contains carbon but not hydrogen,
#'   list the number of carbon atoms,
#'   followed by the remaining elements in alphabetical order.
#' * If the formula contains neither carbon nor hydrogen,
#'   list all elements in alphabetical order.
#'
#' @param f Molecular formula to be sorted.
#'
#' @return A molecular formular in Hill order.
#' @keywords internal
hill_order <- function(f) {
  if ("C" %in% names(f)) {
    other <-
      names(f) %>%
      setdiff(c("C", "H")) %>%
      stringr::str_sort()
    if ("H" %in% names(f))
      f[c("C", "H", other)]
    else
      f[c("C", other)]
  } else
    f[stringr::str_sort(names(f))]
}

#' Generate a molecular formula.
#'
#' Convert a character vector to a molecular formula. `molecular_formula()` is
#' vectorized over `s`, which allows its usage in functions such as
#' [dplyr::mutate()].
#'
#' A molecular formula object is a named integer vector with class `"mol_form"`,
#' whose names (values) correspond to element symbols (counts).
#'
#' @param s Character vector containing a molecular formula.
#' @param simplify If `simplify = TRUE`, return a single molecular formula
#'   object if possible (that is, if `s` contains a single value).
#'
#' @return If `s` is a scalar and `simplify = TRUE`, the function returns a
#'   single molecular formula object. Otherwise, it returns a list of molecular
#'   formulas.
#'
#' @export
#'
#' @examples
#' molecular_formula("C35 H48 Cl1 N3 O10 S1")
#'
#' # spaces and counts with value 1 are optional
#' molecular_formula("C35H48ClN3O10S")
#'
#' # counts may be negative
#' molecular_formula("H-3 N-1")
#'
#' # returns a single molecular formula object
#' molecular_formula("H2SO4")
#'
#' # returns a list with a single element
#' molecular_formula("H2SO4", simplify = FALSE)
#'
#' # returns a list with one element for each parsed string
#' molecular_formula(c("H2 O", "C6 H12 O6"))
molecular_formula <- function(s, simplify = TRUE) {
  res <-
    s %>%
    stringr::str_match_all("([:upper:][:lower:]?)([-+]?[:digit:]*)") %>%
    purrr::imap(
      function(x, i) {
        formula <-
          readr::parse_integer(x[,3]) %>%
          rlang::set_names(x[,2]) %>%
          tidyr::replace_na(1L)
        if (any(duplicated(names(formula))))
          stop("Formula '", s[i], "' contains duplicate elements")
        structure(
          hill_order(formula),
          class = "mol_form"
        )
      }
    )
  if (simplify & length(res) == 1)
    res[[1]]
  else
    res
}

#' Format a molecular formula.
#'
#' Format a molecular formula as a list of all elements and their counts in
#' [Hill order][hill_order()].
#'
#' @param x Molecular formula to be converted.
#' @param prettify If `prettify = TRUE`, include control characters for pretty
#'   printing (blue, positive element counts; red: negative element counts;
#'   gray: empty formulas).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character representation of the molecular formula
#' @export
#'
#' @examples
#' molecular_formula("N-1 H-3") %>% format()
#' molecular_formula("N-1 H-3") %>% format(prettify = TRUE)
#'
#' molecular_formula("P+1 O+3 H+1") %>% format()
#' molecular_formula("P+1 O+3 H+1") %>% format(prettify = TRUE)
#'
#' # arguments in ... are forwarded to the default format() method
#' molecular_formula("H2O") %>% format(justify = "right", width = 50)
format.mol_form <- function(x, prettify = FALSE, ...) {
  if (length(x) == 0) {
    if (prettify)
      res <- crayon::silver("empty")
    else
      res <- "empty"
  } else {
    res <- purrr::imap_chr(
      x,
      ~if (.x > 1) {
        if (prettify)
          paste0(.y, crayon::blue(.x))
        else
          paste0(.y, .x)
      } else if (.x < 0) {
        if (prettify)
          paste0(.y, crayon::red(.x))
        else
          paste0(.y, .x)
      } else {
        .y
      }
    ) %>%
      paste0(collapse = " ")
  }
  format(res, ...)
}

#' Print a molecular formula.
#'
#' Printing yields the string `Molecular formula <atoms>`, where `atoms` lists
#' all elements and their counts in [Hill order][hill_order()].
#'
#' @param x molecular formula to be printed.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function returns the printed molecular formula invisibly.
#' @export
#'
#' @examples
#' molecular_formula("N-1 H-3")
#'
#' molecular_formula("P+1 O+3 H+1")
print.mol_form <- function(x, ...) {
  cat("Molecular formula <", format(x, prettify = TRUE), ">", sep = "")
  invisible(x)
}

#' Add or substract molecular formulas.
#'
#' Depending on the value of op, molecular formulas e1 and e2 are either added
#' (op = `+`) or substracted (op = `-`).
#'
#' @param e1 left operand.
#' @param e2 right operand.
#' @param op operator.
#'
#' @return The sum of or difference between two molecular formulas.
#' @keywords internal
add_sub_formula <- function(e1, e2, op) {
  if (missing(e2)) {
    structure(
      op(unclass(e1)),
      class = "mol_form"
    )
  } else {
    both <- intersect(names(e1), names(e2))
    only_e1 <- setdiff(names(e1), names(e2))
    only_e2 <- setdiff(names(e2), names(e1))
    res <- c(
      purrr::map_int(both, ~op(e1[.x], e2[.x])) %>%
        rlang::set_names(both) %>%
        purrr::keep(~!dplyr::near(.x, 0)),
       e1[only_e1],
       op(e2[only_e2])
    )
    if (is.null(names(res)))
      names(res) <- character(0)
    structure(
      hill_order(res),
      class = "mol_form"
    )
  }
}

#' @describeIn molecular_formula Add two molecular formulas `e1` and `e2`.
#' @param e1 left operand
#' @param e2 right operand
#' @export
`+.mol_form` <- function(e1, e2) {
  add_sub_formula(e1, e2, `+`)
}

#' @describeIn molecular_formula Subtract molecular formula `e2` from `e1`.
#' @export
`-.mol_form` <- function(e1, e2) {
  add_sub_formula(e1, e2, `-`)
}

#' @describeIn molecular_formula Multiply a molecular formula `e1` with an
#'   integer `e2` (or conversely). The multiplier is always converted to an
#'   integer.
#' @export
`*.mol_form` <- function(e1, e2) {
  if (class(e1) == "mol_form") {
    res <- unclass(e1) * as.integer(e2)
  } else {
    res <- unclass(e2) * as.integer(e1)
  }
  structure(
    res[!dplyr::near(res, 0)],
    class = "mol_form"
  )
}

#' Calculate the mass of a molecular formula.
#'
#' @param f Molecular formula whose mass should be calculated.
#' @param mass_set Atomic masses that should be used for mass calculation.
#'   Either a string specifying a predefined mass set (`"average"`, average
#'   atomic masses as defined by IUPAC; `"monoisotopic"`, masses of the most
#'   abundant isotope of each element), or a variable containing a custom mass
#'   set defined by [new_mass_set()].
#'
#' @return A scalar double representing the molecular mass (in dalton).
#' @export
#'
#' @seealso [new_mass_set()] for defining custom mass sets, [`atomic_masses`]
#'   for predefined mass sets
#'
#' @examples
#' molecular_formula("C6 H10 O5") %>% get_mass()
#'
#' molecular_formula("C6 H10 O5") %>% get_mass("monoisotopic")
get_mass <- function(f, mass_set = "average") {
  if (is.character(mass_set)) {
    if (!mass_set %in% names(fragquaxi::atomic_masses))
      stop("Mass set '", mass_set, "' unknown.")
    mass_set <-
      fragquaxi::atomic_masses %>%
      dplyr::select(.data$element, !!rlang::sym(mass_set)) %>%
      tibble::deframe()
  }

  purrr::imap_dbl(
    f,
    ~.x * mass_set[.y]
  ) %>%
    sum()
}


#' Create a new mass set.
#'
#' The custom mass set will be available in subsequent calls of `get_mass()`.
#' You may specify masses of selected elements only. In this case, however,
#' calculating the mass of a formula that contains an element whose mass was not
#' specified will return `NA`.
#'
#' @param masses Named numeric vector of the form `c(element_1 = mass_1, ...,
#'   element_n = mass_n)`.
#' @param inherits_from Mass set from which masses should be inherited if they
#'   are not specified in `masses`. Either `"average"`, `"monoisotopic"` or a
#'   variable containing a mass set previously defined by `new_mass_set()`.
#'
#' @return A named numeric vector containing supplied and inherited element
#'   masses.
#' @export
#'
#' @seealso [get_mass()] for calculating molecular masses, [`atomic_masses`] for
#'   predefined mass sets `"average"` and `"monoisotopic"`
#'
#' @examples
#' my_masses <- new_mass_set(c(H = 10, C = 100, O = 1))
#' molecular_formula("C6 H8 O5") %>% get_mass(my_masses)
#'
#' # returns NA since my_masses does not specify the mass of sulfur
#' molecular_formula("H2 S O4") %>% get_mass(my_masses)
#'
#' # inherit masses from another mass set to provide default values
#' my_masses_complete <- new_mass_set(
#'   c(H = 10, C = 100, O = 1),
#'   inherits_from = "average"
#' )
#' molecular_formula("H2 S O4") %>% get_mass(my_masses_complete)
new_mass_set <- function(masses, inherits_from = NULL) {
  res <-
    masses %>%
    tibble::enframe(name = "element", value = "col_new_masses") %>%
    dplyr::right_join(
      fragquaxi::atomic_masses,
      by = "element"
    )

  if (!is.null(inherits_from)) {
    if (is.character(inherits_from)) {
      if (inherits_from %in% names(fragquaxi::atomic_masses)) {
        res <-
          res %>%
          dplyr::mutate(
            col_new_masses = dplyr::coalesce(
              .data$col_new_masses,
              !!rlang::sym(inherits_from)
            )
          )
      } else {
        warning(
          "Mass set '", inherits_from,
          "' is unknown. Masses will not be inherited."
        )
      }
    } else {
      res <-
        res %>%
        dplyr::left_join(
          tibble::enframe(
            inherits_from,
            name = "element",
            value = "col_inherit_masses"
          ),
          by = "element"
        ) %>%
        dplyr::mutate(
          col_new_masses = dplyr::coalesce(
            .data$col_new_masses,
            .data$col_inherit_masses
          )
        )
    }
  }

  res %>%
    dplyr::select(.data$element, .data$col_new_masses) %>%
    tibble::deframe()
}
