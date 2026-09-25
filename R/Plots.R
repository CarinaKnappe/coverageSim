coverage_plot <- function(cov, y_breaks, plot.margin) {
  if ("frame" %in% colnames(cov)) cov[, frame := factor(frame)]
  cov_plot <- ggplot(cov, aes(y = count, x = position, fill = frame)) +
    geom_col() + ylab("Coverage") +
    scale_y_continuous(breaks=y_breaks) +
    facet_wrap(~ fraction, ncol = 1, strip.position = "right") +
    theme_classic() + theme(legend.position = "none",
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            plot.margin = plot.margin)
  return(cov_plot)
}

#' Plot transcript coordinates window
#' @param tx_window a GRangesList of window
#' @param sim_reads GRanges of reads
#' @param uORFs GRangesList of uORFs
#' @param cds_region GRangesList of cds window
#' @param cov_real data.table of real coverage
#' @param cov_sim data.table of simulated coverage
#' @param break_string spacing x axis, default: "    "
#' @param y_breaks 50
#' @param prop c(5, 2)
#' @param x_lab ggplot x label, default: "Position"
#' @param plot.margin margin(0, 0, 0, 0, "cm")
#' @return a ggplot object
#' @import ggplot2
#' @importFrom grid unit
#' @export
ORF_frame_plot_window <- function(tx_window, sim_reads, uORFs,
                                  cds_region = startRegion(cds_RPL, upstream = 0, downstream = 100),
                                  cov_real = coveragePerTiling(tx_window, merged_gr, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Real"),
                                  cov_sim = coveragePerTiling(tx_window, sim_reads, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Sim"),
                                  break_string = "    ", y_breaks = 50, prop = c(5, 2), x_lab = "Position",
                                  plot.margin = margin(0, 0, 0, 0, "cm")) {
  if (!requireNamespace("RiboCrypt", quietly = TRUE))
    stop("Package RiboCrypt is not installed, which is required for plotting:
         devtools::install_github('m-swirski/RiboCrypt')")
  # Coverage plot
  cov_combined <- rbindlist(list(cov_real, cov_sim))
  cov_plot <- coverage_plot(cov_combined, y_breaks, plot.margin)
  # ORF panel
  display_genes <- c(uORFs, cds_region)
  names(display_genes) <- c(paste0("U", seq(length(uORFs))), "CDS")
  ORF_Panel_data <- gene_model_panel(display_range = tx_window, display_genes,
                                custom_regions = NULL, viewMode = "any")[[1]]
  ORF_Panel <- geneModelPanelPlot(ORF_Panel_data)
  ORF_Panel <- ORF_Panel +
    ylab("") +
    xlab(x_lab) +
    xlim(0, widthPerGroup(tx_window, F)) +
    scale_y_continuous(labels = function(breaks) {rlang::rep_along(breaks, break_string)}) +
    facet_wrap(~ "", ncol = 1, strip.position = "right") +
    theme_classic() +
    theme(strip.text.y = element_text(),
          strip.background = element_blank(),
          plot.margin = plot.margin,
          axis.ticks.y=element_blank(), axis.line.y =element_blank())
  # Merge the two (2 rows)
  plot <- gridExtra::arrangeGrob(cov_plot, ORF_Panel,
                                 layout_matrix = t(matrix(c(rep(1, prop[1]), rep(2, prop[2])), ncol = sum(prop))),
                                 padding = unit(0.0, "line"))
  return(plot)
}

#' From RiboCrypt
#' @param display_range GRanges object
#' @param annotation GRangesList of annotation
#' @param frame numeric, default 1
#' @param custom_regions custom annotation, i.e. uORFs etc
#' @param viewMode character, "tx" for transcript view, "genomic" for genomic view
#' @return list of length 2, data.table of boxes and numeric vector of line locations
gene_model_panel <- function (display_range, annotation, frame = 1, custom_regions,
                              viewMode) {
  if (!is.null(custom_regions)) {
    same_names <- names(custom_regions) %in% names(annotation)
    names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names],
                                               "_1", sep = "")
    annotation <- c(annotation, custom_regions)
  }
  overlaps <- subsetByOverlaps(annotation, display_range,
                               type = ifelse(viewMode == "tx", "within", "any"))
  if (length(overlaps) > 0) {
    plot_width <- widthPerGroup(display_range)
    onames <- rep(names(overlaps), numExonsPerGroup(overlaps,
                                                    FALSE))
    overlaps <- unlistGrl(overlaps)
    names(overlaps) <- onames
    overlaps$rel_frame <- RiboCrypt:::getRelativeFrames(overlaps)
    rel_frame <- RiboCrypt:::getRelativeFrames(overlaps)
    overlaps <- subsetByOverlaps(overlaps, display_range)
    intersections <- RiboCrypt:::trimOverlaps(overlaps, display_range)
    intersections <- groupGRangesBy(intersections)
    locations <- pmapToTranscriptF(intersections, display_range)
    layers <- RiboCrypt:::geneTrackLayer(locations)
    locations <- unlistGrl(locations)
    rel_frame <- RiboCrypt:::getRelativeFrames(overlaps)
    names(rel_frame) <- names(overlaps)
    if (length(rel_frame) != length(locations))
      rel_frame <- RiboCrypt:::selectFrames(rel_frame, locations)
    locations$rel_frame <- rel_frame
    cols <- RiboCrypt:::colour_bars(locations, overlaps, display_range)
    locations <- ranges(locations)
    blocks <- c(start(locations), end(locations))
    names(blocks) <- rep(names(locations), 2)
    blocks <- sort(blocks)
    lines_locations <- blocks[!(blocks %in% c(1, plot_width))]
    rect_locations <- locations
    locations <- resize(locations, width = 1, fix = "center")
    labels_locations <- start(locations)
    too_close <- labels_locations < 0.02 * plot_width
    too_far <- labels_locations > 0.98 * plot_width
    labels_locations[too_close] <- 1
    labels_locations[too_far] <- plot_width
    hjusts <- rep("center", length(labels_locations))
    hjusts[too_close] <- "left"
    hjusts[too_far] <- "right"
    gene_names <- names(locations)
    custom_region_names <- which(names(lines_locations) %in%
                                   names(custom_regions))
    names(lines_locations) <- rep("black", length(lines_locations))
    names(lines_locations)[custom_region_names] <- "orange4"
    result_dt <- data.table(layers = layers, rect_starts = start(rect_locations),
                            rect_ends = end(rect_locations), cols = cols, labels_locations = labels_locations,
                            hjusts = hjusts, gene_names = gene_names)
  }
  else {
    result_dt <- data.table()
    lines_locations <- NULL
  }
  return(list(result_dt, lines_locations))
}

#' From RiboCrypt
#' @param dt a data.table
#' @param frame integer, default 1
#' @param text_size integer, default 2.
#' @return a ggplot object
geneModelPanelPlot <- function(dt, frame = 1, text_size = 2) {
  base_gg <- ggplot(frame = frame) + ylab("") + xlab("") +
    theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), plot.margin = unit(c(0,
                                                               0, 0, 0), "pt")) + scale_x_continuous(expand = c(0,
                                                                                                                0)) + scale_y_continuous(expand = c(0, 0)) + theme(panel.background = element_rect(fill = "white"))
  if (nrow(dt) > 0) {
    suppressWarnings({
      result_plot <- base_gg + geom_rect(data = dt, mapping = aes(ymin = 0 -
                                                                    layers, ymax = 1 - layers, xmin = rect_starts,
                                                                  xmax = rect_ends), fill = dt$cols, alpha = 0.5) +
        geom_text(data = dt, mapping = aes(y = 0.5 -
                                             layers, x = labels_locations, label = gene_names),
                  color = "black", hjust = dt$hjusts, size = text_size)
    })
  }
  else result_plot <- base_gg
  return(result_plot)
}
