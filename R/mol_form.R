#' Internal vctrs methods
#'
#' @import vctrs
#' @keywords internal
#' @name fragquaxi-vctrs
NULL

#' Sort a character vector of chemical elements according to the Hill system.
#'
#' * If the vector contains carbon and hydrogen,
#'   the result is `c("C", "H", others)`.
#' * If the vector contains carbon but not hydrogen,
#'   the result is `c("C", others)`.
#' * If the formula contains neither carbon nor hydrogen,
#'   the result is `c(others)`. In any case, `others` is sorted alphabetically.
#'
#' @param f Character vector to be sorted.
#'
#' @return A character vector sorted in Hill order.
#' @keywords internal
hill_order <- function(el) {
  other <-
    el %>%
    setdiff(c("C", "H")) %>%
    stringr::str_sort()
  c(if ("C" %in% el) "C", if ("H" %in% el) "H", other)
}

new_molecular_formula <- function(df = data.frame()) {
  if (nrow(df) == 0)
    atoms <- list(C = integer(), H = integer())
  else {
    if (!"C" %in% names(df))
      df["C"] <- 0L
    if (!"H" %in% names(df))
      df["H"] <- 0L
    atoms <-
      as.list(df[hill_order(names(df))]) %>%
      purrr::map(tidyr::replace_na, 0L)
  }
  purrr::iwalk(atoms, ~vec_assert(.x, ptype = integer(), arg = .y))
  new_rcrd(atoms, class = "mol_form")
}

methods::setOldClass(c("mol_form", "vctrs_rcrd"))

#' Vector of molecular formulas.
#'
#' Generate a vector of molecular formulas from a character vector.
#'
#' @param s Character vector describing molecular formulas.
#'
#' @return An S3 vector of class `mol_form`.
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
#' # the function is vectorized
#' molecular_formula(c("H2 O", "C6 H12 O6", "C Cl4"))
molecular_formula <- function(s = character()) {
  s %>%
    stringr::str_match_all("([:upper:][:lower:]?)([-+]?[:digit:]*)") %>%
    purrr::imap_dfr(
      function(x, i) {
        if (length(x) == 0) {
          if (s[i] %in% c("", "<empty>"))
            tibble::tibble(.rows = 1)
          else
            stop("Formula '", s[i], "' is invalid.",
                 call. = FALSE)
        } else {
          if (any(duplicated(x[,2])))
            stop("Formula '", s[i], "' contains duplicate elements.",
                 call. = FALSE)
          readr::parse_integer(x[,3]) %>%
            rlang::set_names(x[,2]) %>%
            tidyr::replace_na(1L) %>%
            purrr::keep(~. != 0L) %>%
            {tibble::tibble(!!!., .rows = 1, .name_repair = "minimal")}
        }
      }
    ) %>%
    new_molecular_formula()
}

#' @export
#' @param x An object to test.
#' @rdname molecular_formula
is_molecular_formula <- function(x) {
  inherits(x, "mol_form")
}

#' Format a molecular formula.
#'
#' Format a molecular formula as a list of all elements and their counts in
#' [Hill order][hill_order()].
#'
#' @param x Molecular formula to be formatted
#' @param prettify If `prettify = TRUE`, include control characters for pretty
#'   printing (blue, positive element counts; red: negative element counts;
#'   gray: empty formulas).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character representation of the molecular formula.
#' @export
#'
#' @examples
#' molecular_formula("N-1 H-3") %>% format()
#' molecular_formula("N-1 H-3") %>% format(prettify = TRUE)
#'
#' molecular_formula("P+1 O+3 H+1") %>% format()
#' molecular_formula("P+1 O+3 H+1") %>% format(prettify = TRUE)
format.mol_form <- function(x, prettify = FALSE, ...) {
  purrr::pmap_chr(
    vec_data(x),
    function(...) {
      elements <- purrr::discard(c(...), ~. == 0L)
      if (length(elements) == 0) {
        if (prettify)
          crayon::silver("<empty>")
        else
          "<empty>"
      } else {
        purrr::imap_chr(
          elements,
          ~if (.x > 1L) {
            if (prettify)
              paste0(.y, crayon::blue(.x))
            else
              paste0(.y, .x)
          } else if (.x < 0L) {
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
    }
  )
}

#' @export
vec_ptype_abbr.mol_form <- function(x, ...) {
  "mol"
}

#' @export
vec_ptype_full.mol_form <- function(x, ...) {
  "molecular formula"
}

#' @method vec_ptype2 mol_form
#' @export
#' @export vec_ptype2.mol_form
#' @rdname fragquaxi-vctrs
vec_ptype2.mol_form <- function(x, y, ...) {
  UseMethod("vec_ptype2.mol_form", y)
}

#' @method vec_ptype2.mol_form default
#' @export
vec_ptype2.mol_form.default <- function(x, y, ..., x_arg = "x", y_arg = "y") {
  vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.mol_form mol_form
#' @export
vec_ptype2.mol_form.mol_form <- function(x, y, ...) {
  new_molecular_formula()
}

#' @method vec_ptype2.mol_form character
#' @export
vec_ptype2.mol_form.character <- function(x, y, ...) {
  character()
}

#' @method vec_ptype2.character mol_form
#' @export
vec_ptype2.character.mol_form <- function(x, y, ...) {
  character()
}

#' @method vec_cast mol_form
#' @export
#' @export vec_cast.mol_form
#' @rdname fragquaxi-vctrs
vec_cast.mol_form <- function(x, to, ...) {
  UseMethod("vec_cast.mol_form")
}

#' @method vec_cast.mol_form mol_form
#' @export
vec_cast.mol_form.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}

#' @method vec_cast.mol_form default
#' @export
vec_cast.mol_form.mol_form <- function(x, to, ...) {
  x
}

#' @method vec_cast.mol_form character
#' @export
vec_cast.mol_form.character <- function(x, to, ...) {
  molecular_formula(x)
}

#' @method vec_cast.character mol_form
#' @export
vec_cast.character.mol_form <- function(x, to, ...) {
  format(x)
}

#' @method vec_proxy_equal mol_form
#' @export
vec_proxy_equal.mol_form <- function(x, ...) {
  format(x)
}

#' @method vec_proxy_compare mol_form
#' @export
vec_proxy_compare.mol_form <- function(x, ...) {
  vec_proxy(x) %>%
    dplyr::rename(AC = any_of("C"), AH = any_of("H")) %>%
    purrr::pmap_chr(
      function(...)
        c(...) %>%
        purrr::discard(~. == 0L) %>%
        {stringr::str_glue("{names(.)}{sprintf('%010d', .)}")} %>%
        stringr::str_c(collapse = "")
    )
}

#' @method vec_math mol_form
#' @export
vec_math.mol_form <- function(.fn, .x, ...) {
  switch(.fn,
    sum =
      vec_proxy(.x) %>%
      dplyr::summarise_all(sum) %>%
      dplyr::select_if(~ sum(.) != 0L) %>%
      new_molecular_formula(),
    stop("'", .fn, "()' not supported by molecular formulas.",
         call. = FALSE)
  )
}

#' @method vec_arith mol_form
#' @export
#' @export vec_arith.mol_form
#' @rdname fragquaxi-vctrs
vec_arith.mol_form <- function(op, x, y, ...) {
  UseMethod("vec_arith.mol_form", y)
}

#' @method vec_arith.mol_form default
#' @export
vec_arith.mol_form.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' Add or subtract molecular formulas.
#'
#' Depending on the value of `op`, molecular formulas `x` and `y` are either added
#' (op = `+`) or subtracted (op = `-`).
#'
#' @param x Left operand.
#' @param y Right operand.
#' @param op Operator.
#'
#' @return A molecular formula resulting from the binary operation.
#' @keywords internal
add_sub_mol_form <- function(op, x, y) {
  dplyr::bind_rows(
    dplyr::mutate(vec_proxy(x), id = dplyr::row_number()),
    dplyr::mutate(op(vec_proxy(y)), id = dplyr::row_number())
  ) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarise_all(sum, na.rm = TRUE) %>%
    dplyr::select(-.data$id) %>%
    dplyr::select_if(~ sum(.) != 0L) %>%
    new_molecular_formula()
}

#' Multiply or divide molecular formulas.
#'
#' Depending on the value of `op`, molecular formula `x` is either multiplied by
#' numeric `y` (op = `*`) or integer divided (op = `%/%`). `y` is always cast to
#' integer.
#'
#' @param x Left operand.
#' @param y Right operand.
#' @param op Operator.
#'
#' @return A molecular formula resulting from the binary operation.
#' @keywords internal
mul_div_mol_form <- function(op, x, y) {
  vec_proxy(x) %>%
    dplyr::mutate_all(~op(., allow_lossy_cast(vec_cast(y, integer())))) %>%
    dplyr::select_if(~sum(.) != 0L) %>%
    new_molecular_formula()
}

#' @method vec_arith.mol_form mol_form
#' @export
vec_arith.mol_form.mol_form <- function(op, x, y, ...) {
  switch(op,
    "+" = add_sub_mol_form(`+`, x, y),
    "-" = add_sub_mol_form(`-`, x, y)
  )
}

#' @method vec_arith.mol_form MISSING
#' @export
vec_arith.mol_form.MISSING <- function(op, x, y, ...) {
  switch(op,
    "+" = x,
    "-" = mul_div_mol_form(`*`, x, -1L),
    stop_incompatible_op(op, x, y)
  )
}

#' @method vec_arith.mol_form numeric
#' @export
vec_arith.mol_form.numeric <- function(op, x, y, ...) {
  switch(op,
    "*" = mul_div_mol_form(`*`, x, y),
    "%/%" = mul_div_mol_form(`%/%`, x, y),
    stop_incompatible_op(op, x, y)
  )
}

#' @method vec_arith.numeric mol_form
#' @export
vec_arith.numeric.mol_form <- function(op, x, y, ...) {
  switch(op,
    "*" = mul_div_mol_form(`*`, y, x),
    stop_incompatible_op(op, x, y)
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
