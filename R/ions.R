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
#'   number).}\item{`protein_data`}{Nested column storing the FASTA file name
#'   (`file`) and number of disulfides (`disulfides`) used for calculating its
#'   formula.}\item{`protein_formula`}{Formula calculated from its sequence.}}
#' @export
#'
#' @examples
#' mab_sequence <- system.file("extdata", "mab_sequence.fasta",
#'                             package = "fragquaxi")
#' define_proteins(mab = mab_sequence)
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

#' Create a specification of PTM compositions.
#'
#' This function assembles a post-translational modification (PTM)
#' specification from PTM compositions.
#'
#' ## PTM compositions
#'
#' PTM compositions ("modcoms") are supplied in a data frame where each row
#' corresponds to a single composition. This data frame must comprise a column
#' labeled `modcom_name` with a short description. The remaining columns must be
#' labeled by valid PTM names (see below) and contain the number of the
#' respective PTM found in the respective composition.
#'
#' Valid PTM names include
#' * monosaccharides as provided in the data set [monosaccharides],
#' * PTMs as provided in the data set [ptms],
#' * user-defined PTMs supplied by the user via the argument `other_ptms`.
#'
#' @param x A data frame with PTM compositions (see details).
#' @param other_ptms A named vector (of strings or molecular formulas) defining
#'   additional PTMs like `c(ptm_1 = correction formula 1, ..., ptm_n =
#'   correction formula n)` .
#'
#' @return A data frame that describes one PTM composition per row and comprises
#'   three columns: \describe{\item{`modcom_name`}{Name of each PTM composition
#'   (same as in the input data frame).}\item{`modcom_data`}{Nested column
#'   storing the remaining columns from the input data
#'   frame.}\item{`modcom_formula`}{Correction formula associated with the
#'   respective PTM composition.}}
#' @export
#'
#' @seealso [monosaccharides] and [ptms] for predefined monosaccharides and
#'   PTMs, respectively
#'
#' @examples
#' ptm_compositions <- tibble::tribble(
#'   ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~PHOS, ~foo,
#'   "G0F/G0",        6,       8,    1,    0,     0,
#'   "G0F/G0F",       6,       8,    2,    0,     0,
#'   "G0F/G1F",       7,       8,    2,    0,     0,
#'   "G1F/G1F",       8,       8,    2,    0,     0,
#'   "G1F/G2F",       9,       8,    2,    0,     0,
#'   "G2F/G2F",      10,       8,    2,    0,     0,
#'   "G2F/G2F+P",    10,       8,    2,    1,     0,
#'   "G2F/G2F+foo",  10,       8,    2,    0,     1,
#' )
#'
#' my_ptms <- c(
#'   foo = "C42",
#'   bar = "H42"
#' )
#'
#' define_ptm_compositions(ptm_compositions, my_ptms)
define_ptm_compositions <- function(x, other_ptms = NULL) {
  if (vec_is(other_ptms, character()))
    other_ptms <- molecular_formula(other_ptms)

  requested_mods <-
    names(x) %>%
    setdiff("modcom_name")

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
      tibble::deframe()
  ) %>%
    molecular_formula()

  all_mods <- c(
    all_mods,
    other_ptms[names(other_ptms) %in% requested_mods]
  )

  dplyr::mutate(
    x,
    modcom_formula =
      purrr::pmap(
        x,
        function(modcom_name, ...) {
          purrr::imap(
            list(...),
            ~.x * all_mods[.y]
          ) %>%
            vec_unchop() %>%
            sum()
        }
      ) %>%
      vec_unchop()
  ) %>%
    tidyr::nest(modcom_data = -c(.data$modcom_name, .data$modcom_formula)) %>%
    dplyr::select(.data$modcom_name, .data$modcom_data, .data$modcom_formula)
}

#' Assemble proteoforms and calculate their masses.
#'
#' This function combines a set of proteins with a set of PTM compositions
#' ("modcoms") and calculates the masses of the resulting proteoforms (i.e., the
#' different molecular forms in which the protein product of a single gene can
#' be found due to changes introduced by posttranslational modifications).
#'
#' @param p A protein specification as returned by [define_proteins()].
#' @param m A modcom specification as returned by [define_ptm_compositions()].
#' @param mass_set Atomic masses that should be used for mass calculation.
#' @return A data frame containing \eqn{p \times m}{p * m} proteoforms, obtained
#'   by combining information on \eqn{p} proteins (i.e., number of rows in `p`)
#'   with information on \eqn{m} PTM compositions (i.e., number of rows in `m`).
#'   This data frame comprises two columns derived from the input data frames
#'   (`protein_name` and `modcom_name`, respectively) and two new columns
#'   `formula` (overall molecular formula) and `mass` (calculated proteoform
#'   masses).
#' @export
#'
#' @seealso [`get_mass()`] for predefined mass sets
#'
#' @examples
#' proteins <- define_proteins(
#'   system.file("extdata", "mab_sequence.fasta", package = "fragquaxi"),
#'   .disulfides = 16
#' )
#'
#' modcoms <- define_ptm_compositions(sample_modcoms)
#'
#' assemble_proteoforms(proteins, modcoms)
assemble_proteoforms <- function(p, m, mass_set = "average") {
  p %>%
    dplyr::mutate(modcoms = list(m)) %>%
    tidyr::unnest(.data$modcoms) %>%
    dplyr::mutate(
      formula = .data$protein_formula + .data$modcom_formula,
      mass = get_mass(.data$formula, mass_set)
    ) %>%
    dplyr::select(.data$protein_name, .data$modcom_name,
                  .data$formula, .data$mass)
}

#' Calculate mass-to-charge ratios for proteoform ions.
#'
#' For each value \eqn{z} in `charge_states`, the function calculates
#' * the corresponding mass-to-charge ratio \eqn{m / z + m_c}{m/z + mc}
#'    (\eqn{m_c}{mc}: charge agent mass);
#' * upper and lower tolerances for this ratio, so that these values
#'   form a window of width \eqn{mz \times ppm / 10^6}{mz * ppm / 1e6} centered at \eqn{mz}.
#'
#' The data frame passed to the first argument must contain at least one column
#' named `mass`; each of its value represents the mass of a single proteoform.
#' The data frame may contain any number of additional columns, which will be
#' returned unchanged.
#'
#' @param p A data frame describing proteoform masses (see details).
#' @param charge_states Vector of charge states for which m/z values should be
#'   calculated.
#' @param ppm Width of the m/z window in parts per million.
#' @param charge_agent_mass Mass of the charge carrier. If `NULL`, use the
#'   average IUPAC mass of the proton.
#'
#' @return A data frame containing all columns of the input data frame plus four
#'   additional columns: `z` (charge state), `mz` (exact mass-to-charge ratio),
#'   `mz_min`, and `mz_max`. The latter two columns represent lower and upper
#'   tolerances for the mass-to-charge ratio and are used as integration
#'   boundaries by [`quantify_ions()`]. If the supplied data frame has \eqn{r}
#'   rows and \eqn{c} charge states were specified, the resulting data frame has
#'   \eqn{r * c} rows.
#' @export
#'
#' @examples
#' masses <- data.frame(name = c("p1", "p2"), mass = c(140000, 145000))
#' ionize(masses, 36L:40L)
ionize <- function(p,
                   charge_states,
                   ppm = 400,
                   charge_agent_mass = NULL) {
  p %>%
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
    )
}
