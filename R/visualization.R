#' Plot the total ion current.
#'
#' @param ms_data Mass spectrometry data stored in an mzR object as returned by
#'   [mzR::openMSfile()].
#' @param time_scale Units used for the x-axis.
#'
#' @return A ggplot object describing the created plot.
#' @export
#'
#' @examples
#' ms_data <- mzR::openMSfile(
#'   system.file("extdata", "mzml", "mab1.mzML", package = "fragquaxi")
#' )
#'
#' plot_tic(ms_data)
#'
#' plot_tic(ms_data, time_scale = "s")
plot_tic <- function(ms_data, time_scale = c("min", "s")) {
  vis_data <-
    mzR::tic(ms_data) %>%
    rlang::set_names(c("time", "int"))

  time_scale = match.arg(time_scale)

  if (time_scale == "min") {
    x_name = "time (min)"
  } else {
    vis_data <-
      vis_data %>%
      dplyr::mutate(time = .data$time * 60)
    x_name <- "time (s)"
  }

  p <-
    ggplot2::ggplot(vis_data, ggplot2::aes(.data$time, .data$int)) +
    ggplot2::geom_line() +
    ggplot2::xlab(x_name) +
    ggplot2::ylab("intensity") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )

  p
}


#' Visualize integration boundaries for proteoform ions.
#'
#' This function plots one or several mass spectra overlaid with mass-to-charge
#' ratios and integration boundaries of proteoform ions.
#'
#' @param ms_data Mass spectrometry data stored in an mzR object as returned by
#'   [mzR::openMSfile()].
#' @param ions A data frame describing mass-to-charge ratios of proteoform ions
#'   as returned by [ionize()].
#' @param scans Scan numbers to show (vector of integers). If `NULL`, show all
#'   spectra in the file.
#' @param xlim x-axis limits, given as `c(lower, upper)`. If `NULL`, use the
#'   default scale range.
#' @param plot_elements Character vector specifying which plot elements should
#'   be displayed: \describe{ \item{`"mz"`}{Solid red vertical lines indicating
#'   m/z values.} \item{`"mz_min"`, `"mz_max"`}{Dashed red vertical lines
#'   indicating lower and upper integration boundaries, respectively.}
#'   \item{`"ion_labels"`}{Secondary x-axis, labeled by the names and charges of
#'   the proteoform ions shown.} }
#'
#' @return A ggplot object describing the created plot.
#' @export
#'
#' @examples
#' ms_data <- mzR::openMSfile(
#'   system.file("extdata", "mzml", "mab1.mzML", package = "fragquaxi")
#' )
#' pfm_ions <- tibble::tribble(
#'   ~name,       ~Hex, ~HexNAc, ~Fuc,
#'   "G0F/G0",       6,       8,    1,
#'   "G0F/G0F",      6,       8,    2,
#'   "G0F/G1F",      7,       8,    2,
#'   "G1F/G1F",      8,       8,    2,
#'   "G1F/G2F",      9,       8,    2,
#'   "G2F/G2F",     10,       8,    2,
#' ) %>%
#'   calculate_proteoform_masses("C6464 H9950 N1706 O2014 S44") %>%
#'   ionize(36L)
#'
#' # show all scans
#' plot_ions(ms_data, pfm_ions, xlim = c(4100, 4150))
#'
#' # show selected scans
#' plot_ions(ms_data, pfm_ions, scans = 126:136, xlim = c(4100, 4150))
#'
#' # show a single scan
#' plot_ions(ms_data, pfm_ions, scans = 126, xlim = c(4100, 4150))
#'
#' # only show selected plot elements
#' plot_ions(ms_data, pfm_ions,
#'           scans = 126, xlim = c(4100, 4150),
#'           plot_elements = c("mz", "ion_labels"))
plot_ions <- function(ms_data,
                      ions,
                      scans = NULL,
                      xlim = NULL,
                      plot_elements = c("mz", "mz_min",
                                        "mz_max", "ion_labels")) {
  if (is.null(scans))
    all_spectra <- mzR::spectra(ms_data)
  else if (length(scans) == 1)
    all_spectra <- list(mzR::spectra(ms_data, scans))
  else
    all_spectra <- mzR::spectra(ms_data, scans)

  plot_elements <- match.arg(plot_elements, several.ok = TRUE)

  vis_data <-
    all_spectra %>%
    purrr::map(
      ~.x %>%
        magrittr::set_colnames(c("mz", "int"))
      %>% tibble::as_tibble()
    ) %>%
    rlang::set_names(scans) %>%
    dplyr::bind_rows(.id = "scan")

  if (!is.null(xlim))
    vis_data <-
      vis_data %>%
      dplyr::filter(dplyr::between(.data$mz, xlim[1], xlim[2]))

  if ("mz" %in% plot_elements)
    vline_mz <- ggplot2::geom_vline(
      data = ions,
      mapping = ggplot2::aes(xintercept = .data$mz),
      color = "red",
      alpha = .5,
      na.rm = TRUE
    )
  else
    vline_mz <- NULL

  if ("mz_min" %in% plot_elements)
    vline_mz_min <- ggplot2::geom_vline(
      data = ions,
      mapping = ggplot2::aes(xintercept = .data$mz_min),
      color = "red",
      alpha = .25,
      linetype = "longdash",
      na.rm = TRUE
    )
  else
    vline_mz_min <- NULL

  if ("mz_max" %in% plot_elements)
    vline_mz_max <- ggplot2::geom_vline(
      data = ions,
      mapping = ggplot2::aes(xintercept = .data$mz_max),
      color = "red",
      alpha = .25,
      linetype = "longdash",
      na.rm = TRUE
    )
  else
    vline_mz_max <- NULL

  if ("ion_labels" %in% plot_elements) {
    label_data <-
      ions %>%
      dplyr::mutate(label = stringr::str_glue("{name} (z={z})"))

    ion_label_axis <- ggplot2::sec_axis(
      trans = ~.,
      name = "proteoform",
      breaks = label_data$mz,
      labels = label_data$label
    )
  } else
    ion_label_axis <- ggplot2::waiver()

  p <-
    ggplot2::ggplot(vis_data, ggplot2::aes(.data$mz, .data$int)) +
    ggplot2::geom_line(ggplot2::aes(group = scan), alpha = .5) +
    vline_mz +
    vline_mz_min +
    vline_mz_max +
    ggplot2::scale_x_continuous(
      name = "m/z",
      sec.axis = ion_label_axis,
      limits = xlim
    ) +
    ggplot2::ylab("intensity") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
  p
}
