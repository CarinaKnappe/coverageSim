devtools::load_all(".")
library(data.table)
root <- Sys.getenv("COVSIM_RUST_LARGE_OUTPUT", file.path(
  Sys.getenv("COVSIM_DATA_ROOT", unset = file.path(dirname(getwd()), "coverageSim_data")),
  "runs", "2026-09-09_300tx_codon-bias_run2"))
analyses <- c(control = "control_28/rust_correct", strong_single = "strong_28/rust_correct",
              strong_mixed = "strong_mixed/rust_correct", wrong_offset = "strong_mixed/rust_common15")
profiles <- lapply(file.path(root, analyses), function(x) fread(file.path(x, "RUST_A_site_profiles.csv")))
stats <- lapply(file.path(root, analyses), function(x) fread(file.path(x, "RUST_transcript_stats.csv")))
names(profiles) <- names(stats) <- names(analyses)
targets <- c("AAA", "GAA", "CCA", "CGT", "TTC", "GGT")
positions <- -40:19
get_div <- function(p) {
  observed <- as.matrix(p[, as.character(positions), with = FALSE])
  q <- p$expected / sum(p$expected)
  apply(observed, 2, function(x) { x <- x / sum(x); keep <- x > 0 & q > 0; sum(abs(x[keep] * log2(x[keep] / q[keep]))) })
}
summary <- rbindlist(lapply(names(stats), function(n) {
  s <- stats[[n]]; p <- profiles[[n]]; d <- get_div(p)
  data.table(analysis = n, input_reads = 1500000L, transcripts = nrow(s),
    included = sum(s$rust_included), retained_reads = sum(s$reads),
    frame0_fraction = weighted.mean(s$frame_0_ratio, s$reads),
    sense_codons = nrow(p), a_site_divergence = d[positions == 0],
    target_ratio = mean((p[["0"]] / p$expected)[p$codon %in% targets]))
}))
fwrite(summary, file.path(root, "comparison_summary.tsv"), sep = "\t")
colours <- c("grey45", "#D55E00", "#009E73", "#0072B2")
comparison <- function() {
  par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))
  div <- vapply(profiles, get_div, numeric(60))
  matplot(positions, div, type = "l", lty = 1, lwd = 2, col = colours,
    xlab = "Distance from A-site (codons)", ylab = "Absolute divergence", main = "RUST signal versus control")
  legend("topright", names(profiles), col = colours, lty = 1, cex = .65)
  abline(v = 0, lty = 3)
  barplot(summary$frame0_fraction, names.arg = names(profiles), las = 2,
    col = colours, ylim = c(0, 1), ylab = "Fraction of retained reads", main = "Reading frame 0", cex.names = .7)
  boxplot(lapply(stats, function(s) s$reads), col = colours, las = 2,
    ylab = "Reads per trimmed CDS", main = "300 transcripts per analysis", cex.axis = .7)
  ratios <- vapply(profiles, function(p) p[["0"]] / p$expected, numeric(61))
  matplot(seq_len(61), ratios, type = "p", pch = 16, col = colours,
    xlab = "Sense codon index (alphabetical)", ylab = "Observed / expected",
    main = "A-site codon enrichment")
  abline(h = 1, lty = 3)
}
pdf(file.path(root, "RUST_comparison.pdf"), width = 11, height = 9)
comparison()
for (n in names(profiles)) {
  par(mfrow = c(1,1), mar = c(5, 5, 4, 2))
  p <- profiles[[n]]
  values <- log2(as.matrix(p[, as.character(positions), with = FALSE]) / p$expected)
  values[!is.finite(values)] <- -4
  image(positions, seq_len(61), t(pmax(pmin(values, 4), -4)),
    col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), zlim = c(-4,4),
    xlab = "Distance from A-site (codons)", ylab = "Codon", yaxt = "n",
    main = paste(n, "- log2 RUST ratio (blue -4 / white 0 / red +4)"))
  axis(2, at = seq_len(61), labels = p$codon, las = 2, cex.axis = .5)
  abline(v = 0, lty = 2)
  par(mfrow = c(2,1), mar = c(4, 5, 3, 2))
  ratios <- as.matrix(p[, as.character(positions), with = FALSE]) / p$expected
  line_colours <- rep("grey75", nrow(p))
  line_colours[p$codon %in% targets] <- rainbow(length(targets), s = .8, v = .75)
  matplot(positions, t(ratios), type = "l", lty = 1, col = line_colours,
    xlab = "Distance from A-site (codons)", ylab = "Observed / expected",
    main = paste(n, "- all codons, full ratio range"))
  abline(h = 1, lty = 3); abline(v = 0, lty = 2)
  legend("topright", p$codon[p$codon %in% targets],
    col = line_colours[p$codon %in% targets], lty = 1, cex = .7, ncol = 3)
  plot(positions, get_div(p), type = "l", lwd = 2, col = "#0072B2",
    xlab = "Distance from A-site (codons)", ylab = "Absolute divergence",
    main = "Full divergence range")
  abline(v = 0, lty = 2)
}
dev.off()
png(file.path(root, "RUST_comparison.png"), width = 1500, height = 1200, res = 140)
comparison(); dev.off()
print(summary)
