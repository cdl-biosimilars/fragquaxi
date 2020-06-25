#' Find retention time limits for XIC integration.
#'
#' This internal function is called by [quantify_ions()]. It
#' * converts `rt_limits` to a data frame (since it may be a numeric vector with
#'   two elements or a list).
#' * replaces any infinite value with the minimum or maximum retention time.
#' * finds the scan indices enclosed by the retention time limits.
#'   If a retention time window does not include any scans,
#'   scan indices are `NA`.
#'
#' @param ms_header Header of an `mzR` object, possibly pre-filtered for MS1
#'   scans.
#' @param rt_limits Value of the correspondent argument of [quantify_ions()].
#'
#' @return A data frame with one row per retention time window and the
#'   calculated retention time limits and scan indices in columns `rt_min`,
#'   `rt_max`, and `scans`, respectively.
#' @keywords internal
find_rt_limits <- function(ms_header, rt_limits) {
  if (is.numeric(rt_limits) && length(rt_limits) == 2) {
    rt_limits <- tibble::tibble(
      rt_min = rt_limits[1],
      rt_max = rt_limits[2]
    )
  } else {
    rt_limits <- tibble::as_tibble(rt_limits)
  }

  rt <- ms_header$retentionTime

  rt_limits %>%
    dplyr::mutate(
      rt_min = dplyr::recode(.data$rt_min, `-Inf` = min(rt)),
      rt_max = dplyr::recode(.data$rt_max, `+Inf` = max(rt)),
      id_min = purrr::map_int(
        .data$rt_min,
        function(rt_min)
          purrr::detect_index(rt, ~.x >= rt_min)
      ),
      id_max = purrr::map_int(
        .data$rt_max,
        function(rt_max)
          purrr::detect_index(rt, ~.x <= rt_max, .dir = "backward")
      ),
      scans = purrr::map2(
        .data$id_min,
        .data$id_max,
        function(id_min, id_max) {
          if (id_min > id_max)
            NA_integer_
          else
            ms_header$seqNum[id_min:id_max]
        }
      )
    ) %>%
    dplyr::select(.data$rt_min, .data$rt_max, .data$scans)
}

#' Quantify proteoform ions
#'
#' This function quantifies proteoform ions in a mass spectrometric (MS) data
#' file, which comprises several mass spectra recorded at different time points
#' during a chromatographic run. To this end, `quantify_ions()` employs a
#' two-step integration approach: First, each of the \eqn{i} ions is quantified
#' in each of the \eqn{n} mass spectra by integrating intensity in the m/z
#' tolerance window associated with the respective ion. This step yields an
#' \eqn{i \times n}{i-by-n} matrix of abundances (areas under the curve), with
#' one row per ion and one column per scan. Rows in this matrix correspond to
#' extracted ion currents (XICs). Second, each XIC is integrated within the
#' \eqn{r} given retention time limits. This step yields an \eqn{i \times
#' r}{i-by-r} matrix of abundances, with one row per ion and one column per
#' retention time window.
#'
#' ## Ion specification
#'
#' Ions are provided *via* a data frame that must contain two columns `mz_min`
#' and `mz_max` specifying the lower and upper integration limit for each ion.
#' The data frame may contain additional columns, which are ignored.
#'
#' ## Retention time limits
#'
#' `rt_limits` accepts one of the following:
#' * a numeric vector with two elements specifying a single retention time
#'   window
#' * a data frame with two colums `rt_min` and `rt_max` (other columns are
#'   ignored) specifying one retention time window per row
#' * a data type that can be converted to a data frame, such as a named list
#'
#' The values `-Inf` and `+Inf` are replaced by the minimum and maximum
#' retention time, respectively.
#'
#' ## Integration methods
#'
#' The integration method may be selected separately for the spectrum
#' integration step (by `ifun_spectrum`) and the XIC integration step (by
#' `ifun_xic`).
#' \describe{
#' \item{`"adaptive"`}{Default method for spectrum integration (step 1).
#' Adaptive quadrature via [approxfun()] and [stats::integrate()]. This method
#' is faster, but may yield wrong results if the integrand is zero over nearly
#' all its range (which typically occurs if XICs are integrated over the whole
#' chromatogram). }
#' \item{`"trapezoidal"`}{Default method for XIC integration (step 2).
#' Trapezoidal integration, i.e., the exact area under the polygonal chain
#' formed by (mz or retention time, intensity) tuples. This method is slower,
#' but always gives correct results for XIC integration.}
#' }
#'
#' @param ms_data Mass spectrometric data stored in an `mzR` object as returned
#'   by [mzR::openMSfile()].
#' @param ions A data frame specifying ions (see details).
#' @param rt_limits A specification of retention time limits (see details).
#' @param filter_ms1 If true, only scans from MS level 1 are integrated.
#' @param ifun_spectrum Integration method for spectra (see details).
#' @param ifun_xic Integration method for XICs (see details).
#'
#' @return A list of class `quaxi` with seven elements:
#'   \describe{
#'   \item{`ms_file`}{Path of the MS data file.}
#'   \item{`ifuns`}{Character vector denoting the integration methods used for
#'   spectra and XICs, respectively.}
#'   \item{`ions`}{Data frame passed to the argument `ions` describing \eqn{i}
#'   ions, with an additional column `ion_id` prepended (this column contains a
#'   unique identifier for each ion).}
#'   \item{`rt_limits`}{Data frame with \eqn{r} retention time limits (columns
#'   `rt_min` and `rt_max`) and scan indices enclosed by these limits (columns
#'   `scan_min` and `scan_max`).}
#'   \item{`rt`}{Numeric vector containing \eqn{n} retention times associated
#'   with the \eqn{n} scans in the MS data file.}
#'   \item{`xics`}{An \eqn{i \times n}{i-by-n} matrix containing XIC intensity
#'   values, with one row per ion and one column per scan.}
#'   \item{`abundances`}{An \eqn{i \times r}{i-by-r} matrix containing
#'   abundances, with one row per ion and one column per retention time window.}
#'   }
#'
#' @export
#' @examples
#' ms_data <- mzR::openMSfile(
#'   system.file("extdata", "mzml", "mab1.mzML", package = "fragquaxi")
#' )
#'
#' proteins <- define_proteins(
#'   system.file("extdata", "mab_sequence.fasta", package = "fragquaxi"),
#'   .disulfides = 16
#' )
#'
#' modcoms <- define_ptm_compositions(sample_modcoms)
#'
#' pfm_ions <-
#'   assemble_proteoforms(proteins, modcoms) %>%
#'   ionize(36L:40L)
#'
#' # no rentention time limits: use all spectra
#' quantify_ions(ms_data, pfm_ions)
#'
#' # define a single retention time window
#' quantify_ions(ms_data, pfm_ions, rt_limits = c(240, 330))
#'
#' # define several retention time windows;
#' # instead of a list, you may also pass a data frame to rt_limits
#' quantify_ions(ms_data, pfm_ions,
#'               rt_limits = list(rt_min = c(0, 240), rt_max = c(60, 330)))
quantify_ions <- function(ms_data,
                          ions,
                          rt_limits = c(-Inf, +Inf),
                          filter_ms1 = TRUE,
                          ifun_spectrum = c("adaptive", "trapezoidal"),
                          ifun_xic = c("trapezoidal", "adaptive")) {

  ions <-
    ions %>%
    dplyr::mutate(ion_id = paste0("id_", dplyr::row_number())) %>%
    dplyr::relocate(.data$ion_id)

  ms_header <- mzR::header(ms_data)
  if (filter_ms1) {
    ms_header <-
      ms_header %>%
      dplyr::filter(.data$msLevel == 1L)
  }
  scans <- mzR::spectra(ms_data, ms_header$seqNum)
  rt_limits <- find_rt_limits(ms_header, rt_limits)

  ifun_spectrum <- match.arg(ifun_spectrum)
  if (ifun_spectrum == "adaptive") {
    integrate_spectrum <- function(spectrum) {
      fun <- stats::approxfun(spectrum[,1], spectrum[,2], rule = 2)
      purrr::pmap_dbl(
        ions,
        function(mz_min, mz_max, ...) {
          stats::integrate(fun, mz_min, mz_max, stop.on.error = FALSE)$value
        }
      )
    }
  } else {
    integrate_spectrum <- function(spectrum) {
      purrr::pmap_dbl(
        ions,
        function(mz_min, mz_max, ...) {
          trapz_limits(spectrum[,1], spectrum[,2], mz_min, mz_max)
        }
      )
    }
  }

  ifun_xic <- match.arg(ifun_xic)
  if (ifun_xic == "adaptive") {
    integrate_xic <- function(intensity) {
      fun <- stats::approxfun(ms_header$retentionTime, intensity, rule = 2)
      purrr::pmap_dbl(
        rt_limits,
        function(rt_min, rt_max, ...) {
          stats::integrate(fun, rt_min, rt_max, stop.on.error = FALSE)$value
        }
      )
    }
  } else {
    integrate_xic <- function(intensity) {
      purrr::pmap_dbl(
        rt_limits,
        function(rt_min, rt_max, ...) {
          trapz_limits(ms_header$retentionTime, intensity, rt_min, rt_max)
        }
      )
    }
  }

  xics <-
    scans %>%
    purrr::map(integrate_spectrum) %>%
    purrr::flatten_dbl() %>%
    matrix(ncol = nrow(ms_header)) %>%
    magrittr::set_rownames(paste0("id_", seq_len(nrow(.))))

  abundances <-
    xics %>%
    purrr::array_branch(1) %>%
    purrr::map(integrate_xic) %>%
    purrr::flatten_dbl() %>%
    matrix(nrow = nrow(rt_limits)) %>%
    t() %>%
    magrittr::set_rownames(paste0("id_", seq_len(nrow(.))))

  list(
    ms_file = ms_data@fileName,
    ifuns = c(ifun_spectrum, ifun_xic),
    ions = ions,
    rt_limits = rt_limits,
    rt = ms_header$retentionTime,
    xics = xics,
    abundances = abundances
  ) %>%
    structure(class = "quaxi")
}

#' Printing quantification results.
#'
#' Print a concise summary of quantification results.
#'
#' @param x Object to print.
#' @param ... Other arguments passed on to individual methods.
#'
#' @return The function returns quantification results invisibly.
#' @export
print.quaxi <- function(x, ...) {
  ions <-
    x$ions %>%
    format(n = 5) %>%
    paste0(collapse = "\n")

  rt_limits <-
    x$rt_limits %>%
    format(n = 5) %>%
    paste0(collapse = "\n")

  abundances <-
    dplyr::bind_cols(
      x$rt_limits %>% dplyr::select(.data$rt_min, .data$rt_max),
      x$abundances %>% t() %>% tibble::as_tibble()
    ) %>%
    format(n = 5, width = 80, n_extra = 5) %>%
    paste0(collapse = "\n")

  cli::cli_alert_info(
    c(
      "Abundances of {.val {nrow(x$ions)}} ion{?s} ",
      "quantified in {.val {length(x$rt)}} mass spectr{?um/a} ",
      "using {.val {nrow(x$rt_limits)}} retention time window{?s}."
    )
  )
  cli::cli_h2("Parameters")
  cli::cli_text("MS data file: {.file {x$ms_file}}")
  cli::cli_text("")
  cli::cli_text("Ions:")
  cli::cli_verbatim(ions)
  cli::cli_text("")
  cli::cli_text("Retention time limits:")
  cli::cli_verbatim(rt_limits)
  cli::cli_h2("Results")
  cli::cli_verbatim(abundances)

  invisible(x)
}

#' Coerce quantification results to a data frame.
#'
#' This function converts quantification results obtained from [quantify_ions()]
#' to a nested data frame, which provides data on abundances and extracted ion
#' currents (XIC) for each quantified ion.
#'
#' @param x A `quaxi` object that should be coerced to a data frame.
#' @param abundance_col Name of the column that will contain abundance data. If
#'   `NULL`, do not include abundance data in the output.
#' @param xic_col Name of the column that will contain XIC data. If `NULL`, do
#'   not include XIC data in the output.
#' @param ... Other arguments passed on to individual methods.
#'
#' @return The data frame stored in `x$ions` without column `ion_id` (i.e., a
#'   data frame with one ion per row) and up to two nested columns containing
#'   data on abundances and XICs, respectively.
#'
#'   The former column is labeled according to the value of `abundance_col` and
#'   contains the variables `rt_min`, `rt_max`, `scan_min`, `scan_max`, and
#'   `abundance`.
#'
#'   The latter column is labeled according to the value of `xic_col` and
#'   contains the variables `scan` (scan index), `rt` (corresponding retention
#'   time), and `int` (intensity).
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
#' library(tibble)
#'
#' ms_data <- mzR::openMSfile(
#'   system.file("extdata", "mzml", "mab1.mzML", package = "fragquaxi")
#' )
#'
#' proteins <- define_proteins(
#'   system.file("extdata", "mab_sequence.fasta", package = "fragquaxi"),
#'   .disulfides = 16
#' )
#'
#' modcoms <- define_ptm_compositions(sample_modcoms)
#'
#' pfm_ions <-
#'   assemble_proteoforms(proteins, modcoms) %>%
#'   ionize(36L:40L)
#'
#' abundances <- quantify_ions(ms_data, pfm_ions, c(300, 350))
#'
#' # by default, include data on abundances and XICs in the output
#' as_tibble(abundances)
#'
#' # rename the column containing abundance data
#' # (and omit several columns for readability)
#' as_tibble(abundances, abundance_col = "abundances") %>%
#'   dplyr::select(-formula, -mz_min, -mz_max)
#'
#' # omit XIC data in the output
#' as_tibble(abundances, xic_col = NULL) %>%
#'   dplyr::select(-formula, -mz_min, -mz_max)
as_tibble.quaxi <- function(x,
                            abundance_col = "abundance_data",
                            xic_col = "xic_data",
                            ...) {
  res <- x$ions

  if (!is.null(abundance_col)) {
    abundance_data <-
      dplyr::bind_cols(
        x$rt_limits,
        x$abundances %>%
          t() %>%
          tibble::as_tibble(.name_repair = "unique")
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("id_"),
        names_to = "ion_id",
        values_to = "abundance"
      ) %>%
      tidyr::nest({{abundance_col}} := -.data$ion_id)

    res <- dplyr::left_join(res, abundance_data, by = "ion_id")
  }

  if (!is.null(xic_col)) {
    xic_data <-
      dplyr::bind_cols(
        tibble::tibble(scan = seq_along(x$rt), rt = x$rt),
        x$xics %>%
          t() %>%
          tibble::as_tibble(.name_repair = "unique")
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("id_"),
        names_to = "ion_id",
        values_to = "int"
      ) %>%
      tidyr::nest({{xic_col}} := -.data$ion_id)

    res <- dplyr::left_join(res, xic_data, by = "ion_id")
  }

  dplyr::select(res, -.data$ion_id)
}
