#' Plot the total ion current.
#'
#' @param ms_data Mass spectrometry data stored in an mzR object as returned by
#'   [mzR::openMSfile()].
#' @param time_scale Units used for the x-axis.
#' @param filter_ms1 If true, only scans from MS level 1 are taken into account.
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
plot_tic <- function(ms_data, time_scale = c("min", "s"), filter_ms1 = TRUE) {
  vis_data <-
    mzR::tic(ms_data) %>%
    rlang::set_names(c("time", "int"))

  if (filter_ms1) {
    ms1_scans <-
      mzR::header(ms_data) %>%
      dplyr::filter(msLevel == 1L) %>%
      dplyr::pull(seqNum)
    vis_data <-
      vis_data %>%
      magrittr::extract(ms1_scans,)
  }

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
#' @param rt_limits A numeric vector with two elements specifying retention time
#'   limits (in seconds) for selecting scans. If values are specified for both
#'   `scans` and `rt_limits`, the function issues a warning and ignores the
#'   value of the latter.
#' @param xlim x-axis limits, given as `c(lower, upper)`. If `NULL`, use the
#'   default scale range.
#' @param plot_elements Character vector specifying which plot elements should
#'   be displayed: \describe{ \item{`"mz"`}{Solid red vertical lines indicating
#'   m/z values.} \item{`"mz_min"`, `"mz_max"`}{Dashed red vertical lines
#'   indicating lower and upper integration boundaries, respectively.}
#'   \item{`"ion_labels"`}{Secondary x-axis, labeled by the names and charges of
#'   the proteoform ions shown.} }
#' @param filter_ms1 If true, only scans from MS level 1 are selected by the
#'   retention time limits.
#'
#' @return A ggplot object describing the created plot.
#' @export
#'
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
                      rt_limits = NULL,
                      xlim = NULL,
                      plot_elements = c("mz", "mz_min",
                                        "mz_max", "ion_labels"),
                      filter_ms1 = TRUE) {
  if (!is.null(rt_limits)) {
    if (!is.null(scans)) {
      warning("Values specified for both 'scans' and 'rt_limits'. ",
              "Ignoring the value of 'rt_limits'.",
              call. = FALSE)
    } else {
      scans <-
        mzR::header(ms_data) %>%
        filter(retentionTime %>% dplyr::between(rt_limits[1], rt_limits[2]))

      if (filter_ms1) {
        scans <-
         scans %>%
          dplyr::filter(msLevel == 1L)
      }

      scans <-
        scans %>%
        dplyr::pull(seqNum)

      if (length(scans) == 0)
        stop("No scans within the given retention time limits.", call. = FALSE)

      message("Showing the following scans ",
              "based on the given retention time limits: ",
              paste(scans, collapse = ", "))
    }
  }

  if (is.null(scans))
    all_spectra <- mzR::spectra(ms_data)
  else if (length(scans) == 1)
    all_spectra <- list(mzR::spectra(ms_data, scans))
  else
    all_spectra <- mzR::spectra(ms_data, scans)

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
      dplyr::mutate(
        label = stringr::str_glue("{protein_name}_{modcom_name} (z={z})")
      )

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

#' Diagnostic plot for selecting optimal peak integration boundaries.
#'
#' This function draws proteoform abundances vs the PPM value used by
#' [ionize()], which in turn determines the integration boundaries for each ion.
#'
#' By default, the function creates four plots: Total proteoform abundance vs
#' PPM (upper left), first derivative of this curve (lower left), individual
#' proteoform abundances vs PPM (upper right) and first derivative of the latter
#' curve (lower right).
#'
#' Optimal integration boundaries capture each peak completely and do not lead
#' to overlaps with neighboring peaks. This corresponds to the PPM value where
#' the curves in the upper plots reach the first plateau, and their derivatives
#' become zero.
#'
#' @param ms_data Mass spectrometry data stored in an mzR object as returned by
#'   [mzR::openMSfile()].
#' @param proteoforms Proteoform masses as returned by [assemble_proteoforms()].
#' @param charge_states Vector of charge states in which the proteoforms should
#'   be quantified.
#' @param ppm_values Vectors of PPM values where abundances should be sampled.
#' @param subplots Subplots that should be drawn; one of `"all"` (draw all
#'   subplots), `"total"` (draw only the total abundance plot), and `"modcom"`
#'   (draw only the subplot with individual proteoform abundances).
#'
#' @return A ggplot object describing the created plot.
#'
#' @export
#'
#' @examples
#' ms_data <- mzR::openMSfile(
#'   system.file("extdata", "mzml", "mab1.mzML", package = "fragquaxi")
#' )
#' proteins <- define_proteins(
#'   system.file("extdata", "mab_sequence.fasta", package = "fragquaxi"),
#'   .disulfides = 16
#' )
#' modcoms <- define_ptm_compositions(sample_modcoms)
#' pfm_masses <- assemble_proteoforms(proteins, modcoms)
#'
#' plot_optimal_ppm(ms_data, pfm_masses, 33L:40L)
plot_optimal_ppm <- function(ms_data,
                             proteoforms,
                             charge_states,
                             ppm_values = seq(0, 2500, 50),
                             subplots = c("all", "total", "modcom")) {
  subplots <- match.arg(subplots)

  pfm_ions <-
    purrr::map(
      ppm_values,
      ~ionize(proteoforms, charge_states = charge_states, ppm = .x)
    ) %>%
    rlang::set_names(ppm_values)

  abundances <-
    purrr::map_dfr(
      pfm_ions,
      ~quantify_ions(
        ms_data,
        .x,
        rt_limits = c(300, 350)
      ) %>%
        tibble::as_tibble() %>%
        dplyr::select(.data$modcom_name, .data$abundance_data),
      .id = "ppm"
    ) %>%
    dplyr::mutate(ppm = as.integer(.data$ppm)) %>%
    tidyr::unnest(.data$abundance_data) %>%
    dplyr::group_by(.data$ppm, .data$modcom_name) %>%
    dplyr::summarise(abundance = sum(.data$abundance))

  delta_abundances <-
    abundances %>%
    dplyr::group_split() %>%
    {purrr::map2(
      .[1:length(.) - 1],
      .[2:length(.)],
      ~dplyr::left_join(.y, .x, by = c("modcom_name")) %>%
        dplyr::mutate(
          delta_abundance = .data$abundance.x - .data$abundance.y
        ) %>%
        dplyr::select(
          ppm = .data$ppm.x,
          .data$modcom_name,
          .data$delta_abundance
        )
    )} %>%
    dplyr::bind_rows()

  plot_total <-
    dplyr::left_join(
      abundances,
      delta_abundances,
      by = c("ppm", "modcom_name")
    ) %>%
    dplyr::group_by(.data$ppm) %>%
    dplyr::summarise(
      abundance = sum(.data$abundance),
      delta_abundance = sum(.data$delta_abundance)
    ) %>%
    dplyr::mutate(modcom_name = "total") %>%
    tidyr::pivot_longer(
      dplyr::contains("abundance"),
      names_to = "type",
      values_to = "area"
    ) %>%
    dplyr::mutate(
      type = dplyr::recode(.data$type, delta_abundance = "first derivative")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(.data$ppm, .data$area)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$type),
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::ggtitle("Total abundance")

  plot_modcom <-
    dplyr::left_join(
      abundances,
      delta_abundances,
      by = c("ppm", "modcom_name")
    ) %>%
    tidyr::pivot_longer(
      dplyr::contains("abundance"),
      names_to = "type",
      values_to = "area"
    ) %>%
    dplyr::mutate(
      type = dplyr::recode(.data$type, delta_abundance = "first derivative")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(.data$ppm, .data$area)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$modcom_name)) +
    ggplot2::scale_color_hue(name = "proteoform") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$type),
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::ggtitle("Proteoform abundances")

  if (subplots == "total")
    p <- plot_total
  else if (subplots == "modcom")
    p <- plot_modcom
  else
    p <- patchwork::wrap_plots(plot_total, plot_modcom)

  p
}

