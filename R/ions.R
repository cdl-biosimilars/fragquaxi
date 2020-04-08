#' Calculate mass-to-charge ratios.
#'
#' The function is vectorized over `z` and thus may calculate several \eqn{m/z}
#' values for a single mass.
#'
#' @param m Mass.
#' @param z Vector of charges.
#' @param charge_agent_mass Mass of the charge carrier. If `NULL`, use the
#'   average IUPAC mass of the proton.
#'
#' @return A named vector of mass-to-charge ratios; names correspond to the
#'   values of `z`, converted to character.
#' @export
#'
#' @examples
#' charge(150000, z = 24)
#'
#' charge(150000, z = 20:24)
#'
#' charge(150000, z = 20:24, charge_agent_mass = 3)
charge <- function(m, z, charge_agent_mass = NULL) {
  if (is.null(charge_agent_mass))
    charge_agent_mass <-
      fragquaxi::atomic_masses %>%
      dplyr::filter(.data$element == "H") %>%
      dplyr::pull(.data$average)

  if (any(z < 1))
    stop("z < 1. ",
         "Use a negative charge carrier to obtain negative charge states.")

  rlang::set_names(
    m / z + charge_agent_mass,
    z
  )
}

#' Create a protein specification.
#'
#' This function assembles a protein specification from several FASTA files.
#' It converts the amino acid sequence stored in each FASTA file to a
#' molecular formula. If the protein has multiple chains, each chain must be
#' specified with its own header line, *i.e.*, a line starting with `>`. For
#' each chain, a water molecule is added to the formula.
#'
#' @param ... FASTA files containing protein sequences. Argument names will be
#'   used as protein names.
#' @param .disulfides A vector describing the number of disulfide bridges,
#'   recycled to the number of the given FASTA files. For each disulfide bridge,
#'   two hydrogen atoms are subtracted from the molecular formula.
#'
#' @return A data frame that describes one protein per row and comprises three
#'   columns: \describe{\item{`protein_name`}{Name of each protein, derived from
#'   the argument name in `...` (if present; otherwise, a consecutive
#'   number).}\item{`protein_data`}{Data frame storing the FASTA file name
#'   (`file`) and number of disulfides (`disulfides`) used for calculating its
#'   formula.}\item{`protein_formula`}{Formula calculated from its sequence.}}
#' @export
#'
#' @examples
#' mab_sequence <- system.file("extdata", "mab_sequence.fasta",
#'                             package = "fragquaxi")
#' define_proteins(mab_sequence)
#'
#' define_proteins(mab_sequence, .disulfides = 16)
define_proteins <- function(..., .disulfides = 0L) {
  files <- c(...)
  disulfides <- vec_recycle(.disulfides, vec_size(files), x_arg = ".disulfides")

  amino_acid_compositions <-
    fragquaxi::amino_acids$formula %>%
    molecular_formula() %>%
    rlang::set_names(fragquaxi::amino_acids$abbreviation)

  purrr::map2(
    files,
    disulfides,
    function(path, n_disulfides) {
      sequences <- Biostrings::readAAStringSet(path)
      res <-
        sequences %>%
        Biostrings::letterFrequency(names(amino_acid_compositions)) %>%
        colSums() %>%
        purrr::imap(~amino_acid_compositions[[.y]] * .x) %>%
        vec_unchop() %>%
        sum()
      res +
        molecular_formula("H2O") * length(sequences) -
        molecular_formula("H2") * n_disulfides
    }
  ) %>%
    vec_unchop() %>%
    tibble::enframe(name = "protein_name", value = "protein_formula") %>%
    dplyr::mutate(file = files, disulfides = disulfides) %>%
    tidyr::nest(protein_data = c(.data$file, .data$disulfides)) %>%
    dplyr::select(.data$protein_name, .data$protein_data, .data$protein_formula)
}

#' Create a proteoform specification.
#'
#' ## Proteoform specification
#'
#' The set of proteoforms is described by a data frame passed to the first
#' argument. In this data frame, each row corresponds to a single proteoform,
#' and the columns specify each proteoform as follows:
#'
#' * A mandatory column labeled `name` must contain a short, informative
#' description.
#' * An optional column labeled `protein_formula` specifies a formula for the
#' protein associated with the respective proteoform (useful if the data frame
#' describes proteoforms derived from different proteins). If this column is
#' absent, all rows use the value of the argument `protein_formula`.
#' * Several optional columns labeled by valid PTM names (see below) contain the
#' number of the respective PTM found in the respective proteoform.
#'
#' Valid PTM names include
#'
#' * monosaccharides as provided in the data set [monosaccharides],
#' * PTMs as provided in the data set [ptms],
#' * user-defined PTMs supplied by the user via the argument `other_ptms`.
#'
#' ## Protein formula
#'
#' Each of the \eqn{n} elements in `protein_formula` is interpreted as a
#' distinct protein for which all \eqn{p} proteoforms specified *via* the first
#' argument should be constructed. Hence, the result will contain \eqn{n \times
#' p}{n * p} rows.
#'
#' ## Composite proteoform names
#'
#' In order to ensure unique identifiers in column `name`, its contents will be
#' rewritten if `protein_formula` contains more than one element. The new names
#' will be assembled according to the [glue::glue()] string literal defined by
#' `name_template`. In this string literal,
#'
#' * `{name}` will be replaced by the original value in column `name`, and
#' * `{fname}` will be replaced by the names of vector `protein_formula` (if
#'   present) or by consecutive numbers.
#'
#' @param proteoforms A data frame describing proteoforms (see details).
#' @param protein_formula A string or [molecular formula][molecular_formula()]
#'   describing the molecular formula of the protein (see details). If `NULL`,
#'   the data frame passed to `proteoforms` must contain a column
#'   `protein_formula` with formulas.
#' @param mass_set Atomic masses that should be used for mass calculation.
#' @param other_ptms A named character vector of the form `c(ptm_1 = "correction
#'   formula 1", ..., ptm_n = "correction formula n")` defining additional PTMs.
#' @param name_template A [glue::glue()] string literal that describes the
#'   construction of composite proteoform names (see details).
#' @return A data frame containing a complete proteoform specification. It
#'   comprises all columns of the supplied data frame plus the additional
#'   columns `protein_formula` (if absent in the input), `ptm_formula` (sum
#'   formula of all PTMs), and `mass` (calculated proteoform masses).
#' @export
#'
#' @seealso [`get_mass()`] for predefined mass sets, [monosaccharides] and
#'   [ptms] for predefined monosaccharides and PTMs, respectively
#'
#' @examples
#' proteoforms <- tibble::tribble(
#'   ~name,       ~Hex, ~HexNAc, ~Fuc, ~PHOS, ~foo,
#'   "G0F/G0",       6,       8,    1,    0,     0,
#'   "G0F/G0F",      6,       8,    2,    0,     0,
#'   "G0F/G1F",      7,       8,    2,    0,     0,
#'   "G1F/G1F",      8,       8,    2,    0,     0,
#'   "G1F/G2F",      9,       8,    2,    0,     0,
#'   "G2F/G2F",     10,       8,    2,    0,     0,
#'   "G2F/G2F+P",   10,       8,    2,    1,     0,
#'   "G2F/G2F+foo", 10,       8,    2,    0,     1,
#' )
#'
#' mab_atoms <- "C6464 H9950 N1706 O2014 S44"
#' my_ptms <- c(
#'   foo = "C42",
#'   bar = "H42"
#' )
#'
#' calculate_proteoform_masses(proteoforms, mab_atoms, other_ptms = my_ptms)
#'
#' mab_formula <- molecular_formula(mab_atoms)
#' calculate_proteoform_masses(proteoforms, mab_formula, other_ptms = my_ptms)
calculate_proteoform_masses <- function(proteoforms,
                                        protein_formula = NULL,
                                        mass_set = "average",
                                        other_ptms = NULL,
                                        name_template = "{fname}_{name}") {
  requested_mods <-
    names(proteoforms) %>%
    setdiff(c("name", "protein_formula"))

  unknown_mods <-
    requested_mods %>%
    setdiff(fragquaxi::monosaccharides$abbreviation) %>%
    setdiff(fragquaxi::ptms$abbreviation) %>%
    setdiff(names(other_ptms))
  if (length(unknown_mods) > 0)
    stop(
      "Unknown modifications: ",
      paste(unknown_mods, collapse = " ")
    )

  duplicate_mods <-
    fragquaxi::monosaccharides$abbreviation %>%
    intersect(fragquaxi::ptms$abbreviation) %>%
    intersect(names(other_ptms))
  if (length(duplicate_mods) > 0)
    stop(
      "Modifications with duplicate definition: ",
      paste(duplicate_mods, collapse = " ")
    )

  all_mods <- c(
    fragquaxi::monosaccharides %>%
      dplyr::filter(.data$abbreviation %in% requested_mods) %>%
      dplyr::select(.data$abbreviation, .data$formula) %>%
      tibble::deframe(),
    fragquaxi::ptms %>%
      dplyr::filter(.data$abbreviation %in% requested_mods) %>%
      dplyr::select(.data$abbreviation, .data$formula) %>%
      tibble::deframe(),
    other_ptms[names(other_ptms) %in% requested_mods]
  )
  all_mods <-
    molecular_formula(all_mods) %>%
    rlang::set_names(names(all_mods))

  if ("protein_formula" %in% names(proteoforms)) {
    if (!is.null(protein_formula))
      warning("Protein formulas already specified for proteoforms. ",
              "Ignoring value of argument 'protein_formula'.",
              call. = FALSE)

    if (vec_is(proteoforms$protein_formula, character()))
      proteoforms <-
        proteoforms %>%
        dplyr::mutate(
          protein_formula = molecular_formula(.data$protein_formula)
        )
  } else {
    if (is.null(protein_formula))
      stop("No protein formula specified. ",
           "Either specify a formula via the argument 'protein_formula' or ",
           "include a column 'protein_formula' in the data frame ",
           "passed to 'proteoforms'.",
           call. = FALSE)

    if (vec_is(protein_formula, character()))
      protein_formula <- molecular_formula(protein_formula)

    protein_formula <- tibble::enframe(
      protein_formula,
      name = "fname",
      value = "protein_formula"
    )

    proteoforms <-
      proteoforms %>%
      dplyr::mutate(protein_formula = list(.env$protein_formula)) %>%
      tidyr::unnest(.data$protein_formula)

    if (nrow(protein_formula) > 1)
      proteoforms <-
        proteoforms %>%
        dplyr::mutate(name = stringr::str_glue(name_template))

    proteoforms <-
      proteoforms %>%
      dplyr::select(-.data$fname)
  }

  proteoforms %>%
    dplyr::mutate(
      ptm_formula =
        purrr::pmap(
          .,
          function(name, protein_formula, ...) {
            purrr::imap(
              list(...),
              ~.x * all_mods[.y]
            ) %>%
              vec_unchop() %>%
              sum()
          }
        ) %>%
        vec_unchop(),
      mass =
        get_mass(.data$protein_formula, mass_set) +
        get_mass(.data$ptm_formula, mass_set)
    ) %>%
    dplyr::select(.data$name, names(.env$all_mods), dplyr::everything())
}

#' Calculate mass-to-charge ratios for proteoform ions.
#'
#' For each value \eqn{z} in `charge_states`, the function calculates
#' * the corresponding mass-to-charge ratio \eqn{m / z + mc}
#'    (\eqn{mc}: charge agent mass);
#' * upper and lower tolerances for this ratio, so that these values
#'   form a window of width \eqn{mz * ppm / 1000000} centered at \eqn{mz}.
#'
#' The data frame passed to `proteoforms` must contain at least one column named
#' `mass`; each of its value represents the mass of a single proteoform. The
#' data frame may contain any number of additional columns, which will be
#' returned unchanged.
#'
#' @param proteoforms A data frame describing proteoform masses (see details).
#' @param charge_states Vector of charge states for which m/z values should be
#'   calculated.
#' @param ppm Width of the m/z window in parts per million.
#' @param charge_agent_mass Mass of the charge carrier. If `NULL`, use the
#'   average IUPAC mass of the proton.
#'
#' @return A data frame containing all columns of the supplied data frame plus
#'   four additional columns: `z` (charge state), `mz` (exact mass-to-charge
#'   ratio), `mz_min`, and `mz_max`. The latter two columns represent lower and
#'   upper tolerances for the mass-to-charge ratio and are used as integration
#'   boundaries by [`quantify_ions()`]. If the supplied data frame has \eqn{r}
#'   rows and \eqn{c} charge states were specified, the resulting data frame has
#'   \eqn{r * c} rows.
#' @export
#'
#' @examples
#' masses <- data.frame(name = c("p1", "p2"), mass = c(140000, 145000))
#' ionize(masses, 36L:40L)
ionize <- function(proteoforms,
                   charge_states,
                   ppm = 400,
                   charge_agent_mass = NULL) {
  proteoforms %>%
    dplyr::mutate(mz_data = purrr::map(.data$mass, charge, charge_states)) %>%
    tidyr::unnest_longer(.data$mz_data, indices_to = "z", values_to = "mz") %>%
    dplyr::mutate(
      z = readr::parse_guess(.data$z, guess_integer = TRUE),
      mz_min = .data$mz * (1 - ppm / 2e6),
      mz_max = .data$mz * (1 + ppm / 2e6)
    ) %>%
    dplyr::select(
      -tidyselect::all_of(c("z", "mz", "mz_min", "mz_max")),
      .data$z, .data$mz, .data$mz_min, .data$mz_max
    ) %>%
    dplyr::arrange(.data$name, .data$z)
}
