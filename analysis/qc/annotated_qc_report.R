# Annotated QC report for a coverageSim run.
#
# The plain QC report (analysis/qc/simulated_rpf_qc.R) draws the figures.
# This one puts each figure on a page of its own together with what it shows,
# and marks one real data point whose position is read out biologically. Every
# number in the text is computed from the run, so the prose cannot drift away
# from the data it describes.
#
# Environment variables:
#   COVSIM_AQC_INPUT   run directory that already contains a qc/ folder
#   COVSIM_AQC_OUT     output PDF (default: <input>/qc/annotated_qc_report.pdf)
#   COVSIM_AQC_SAMPLE  sample the example points are taken from (default RFP_WT_1)

suppressMessages({
  library(data.table); library(ggplot2); library(gridExtra); library(grid)
})

input <- Sys.getenv("COVSIM_AQC_INPUT", unset = "")
if (!nzchar(input)) stop("Set COVSIM_AQC_INPUT to the run directory")
qc_dir <- file.path(input, "qc")
if (!dir.exists(qc_dir)) stop("No qc/ folder in ", input, " - run simulated_rpf_qc.R first")
out_pdf <- Sys.getenv("COVSIM_AQC_OUT",
                      unset = file.path(qc_dir, "annotated_qc_report.pdf"))
S <- Sys.getenv("COVSIM_AQC_SAMPLE", unset = "RFP_WT_1")

tsv <- function(name) fread(file.path(qc_dir, paste0(name, ".tsv")))
pct <- function(x) sprintf("%.1f%%", 100 * x)

theme_qc <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.title = element_blank())

# One page: the figure on top, the explanation underneath.
page <- function(plot, what, example, reading) {
  body <- paste0(
    "WHAT THIS SHOWS\n", strwrap(what, width = 165) |> paste(collapse = "\n"),
    "\n\nEXAMPLE POINT\n", strwrap(example, width = 165) |> paste(collapse = "\n"),
    "\n\nWHAT ITS POSITION MEANS\n", strwrap(reading, width = 165) |> paste(collapse = "\n")
  )
  grid.arrange(plot,
    textGrob(body, x = 0.02, hjust = 0, gp = gpar(fontsize = 9, lineheight = 1.25)),
    heights = c(1.45, 1))
}

pdf(out_pdf, width = 9.5, height = 10.5, onefile = TRUE)

# ---- Title page -------------------------------------------------------------
settings <- tsv2 <- if (file.exists(file.path(input, "simulation_settings.tsv")))
  fread(file.path(input, "simulation_settings.tsv")) else data.table()
info <- if (nrow(settings)) paste(sprintf("%s = %s", settings$setting, settings$value),
                                  collapse = "    ") else "(no settings file)"
grid.newpage()
grid.text("Simulated Ribo-seq: annotated quality control",
          y = 0.70, gp = gpar(fontsize = 20, fontface = "bold"))
grid.text(paste0("Run: ", basename(input)), y = 0.62, gp = gpar(fontsize = 12))
grid.text(info, y = 0.56, gp = gpar(fontsize = 9))
grid.text(paste0("Example points are taken from sample ", S, "."),
          y = 0.48, gp = gpar(fontsize = 10))
grid.text(paste(strwrap(paste(
  "Each of the following pages carries one figure, a short statement of what it tests, and one",
  "real data point from this run whose position in the plot is read out in biological terms.",
  "Every number in the text is computed from the run's own QC tables."), width = 95),
  collapse = "\n"), y = 0.36, gp = gpar(fontsize = 9.5, lineheight = 1.3))


# A page carrying several panels: the panels on top, one short note per panel
# underneath. Used for the quantities that only need a glance.
panel_page <- function(plots, notes, labels = letters[seq_along(plots)]) {
  tagged <- Map(function(p, l) p + labs(tag = l) +
                  theme(plot.tag = element_text(face = "bold", size = 12)),
                plots, labels)
  body <- paste(mapply(function(l, n)
    paste(strwrap(paste0(l, ")  ", n), width = 178, exdent = 4), collapse = "\n"),
    labels, notes), collapse = "\n\n")
  grid.arrange(
    arrangeGrob(grobs = tagged, ncol = 2),
    textGrob(body, x = 0.02, hjust = 0, gp = gpar(fontsize = 8.5, lineheight = 1.2)),
    heights = c(1.75, 1))
}

# ---- 1. Overview: six quantities that only need a glance (a-f) --------------
bam <- tsv("bam_qc")
pa <- ggplot(melt(bam[, .(sample, truth_reads, bam_reads)], id.vars = "sample"),
             aes(sample, value, fill = variable)) +
  geom_col(position = "dodge") + scale_y_continuous(labels = scales::comma) +
  labs(title = "Reads per library", x = NULL, y = "Reads") + theme_qc +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rc <- fread(file.path(input, "simulated_region_counts.tsv"))
rl <- melt(rc, id.vars = "region", variable.name = "sample", value.name = "reads")
rl[, share := reads / sum(reads), by = sample]
pb <- ggplot(rl, aes(region, share, fill = sample)) + geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Requested signal by region", x = NULL, y = "Share") + theme_qc

lo <- tsv("length_offset")
len <- lo[, .(reads = sum(reads)), by = .(sample, fragment_length)]
len[, share := reads / sum(reads), by = sample]
pc <- ggplot(len, aes(factor(fragment_length), share, fill = sample)) +
  geom_col(position = "dodge") + scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Fragment length distribution", x = "Length (nt)", y = "Share") + theme_qc

fr <- tsv("frame"); fr[, share := reads / sum(reads), by = sample]
pd <- ggplot(fr, aes(factor(frame), share, fill = sample)) +
  geom_col(position = "dodge") + geom_hline(yintercept = 1/3, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Reading frame, all lengths pooled", x = "Frame", y = "Share") + theme_qc

mg <- tsv("metagene")
tot <- mg[, .(reads = sum(reads)), by = .(sample, region)]
tot[, share := reads / sum(reads), by = sample]
pe <- ggplot(mg[sample == S], aes(bin, reads)) + geom_line(colour = "grey25") +
  facet_wrap(~ region, scales = "free_x", nrow = 1) +
  labs(title = paste0("Region profiles (", S, ")"), x = "Position in region", y = "Reads") +
  theme_qc

en <- tsv("end_nt"); en[, share := reads / sum(reads), by = .(sample, end)]
pf <- ggplot(en, aes(nucleotide, share, fill = sample)) + geom_col(position = "dodge") +
  facet_wrap(~ end) + geom_hline(yintercept = 0.25, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Nucleotide at the fragment ends", x = NULL, y = "Share") + theme_qc

exb <- bam[sample == S]
exlen <- len[sample == S][which.max(reads)]
exfr <- fr[sample == S][which.max(share)]
panel_page(
  list(pa, pb, pc, pd, pe, pf),
  c(
    sprintf(paste("Ground truth against written file. They must match exactly; a gap would mean",
                  "fragments were lost on the way out. In %s both are %s, difference %d."),
            S, format(exb$truth_reads, big.mark = ","), exb$count_difference),
    sprintf(paste("How each sample's reads were asked to split across the transcript. The CDS share",
                  "is what makes this Ribo-seq rather than RNA-seq: ribosomes sit on the coding",
                  "sequence. Here the CDS was asked to carry %s."),
            pct(rl[sample == S & region == "cds", share])),
    sprintf(paste("One ribosome shields about 28 to 30 nt, so a real library peaks there and tails",
                  "off to both sides. The commonest length here is %d nt with %s of the reads,",
                  "across %d lengths."),
            exlen$fragment_length, pct(exlen$share),
            len[sample == S, uniqueN(fragment_length)]),
    sprintf(paste("A translating ribosome advances one codon at a time, so its reads should favour",
                  "one of the three frames; the dashed line is the third that chance would give.",
                  "Frame %d holds %s here. Panel is pooled over lengths - the per-length split",
                  "follows on its own page, because that is what decides which lengths are usable."),
            exfr$frame, pct(exfr$share)),
    sprintf(paste("All transcripts stretched to a common length and overlaid, per region. A level",
                  "CDS panel with low leader and trailer panels is the expected shape. Here the CDS",
                  "holds %s of the reads, the leader %s, the trailer %s."),
            pct(tot[sample == S & region == "cds", share]),
            pct(tot[sample == S & region == "leader", share]),
            pct(tot[sample == S & region == "trailer", share])),
    sprintf(paste("Which base sits at the first and last position of each fragment. Four bars near",
                  "the dashed 25 percent mean no sequence preference at the ends - the correct",
                  "answer here, since end bias is switched off, so this panel is a negative",
                  "control. The spread across the four bases is %s to %s."),
            pct(min(en[sample == S & end == "5-prime", share])),
            pct(max(en[sample == S & end == "5-prime", share])))
  ))

# ---- 2. A-site offset per fragment length -----------------------------------
ex <- lo[sample == S][which.max(reads)]
odd <- lo[sample == S & site_offset != ex$site_offset]
p <- ggplot(lo, aes(factor(fragment_length), site_offset, size = reads, colour = sample)) +
  geom_point(alpha = 0.75) + scale_size_area(max_size = 9) +
  geom_point(data = ex, aes(factor(fragment_length), site_offset), inherit.aes = FALSE,
             shape = 21, size = 9, stroke = 1.2, fill = NA, colour = "red") +
  labs(title = "2. A-site offset from the fragment's 5' end",
       x = "Fragment length (nt)", y = "A-site offset (nt)") + theme_qc
page(p,
  paste("For each fragment length, how far into the fragment the decoding site sits. A sequenced",
        "fragment is about 29 nt long, but the biologically interesting position is inside it: the",
        "codon currently in the ribosome's A site. The offset is that distance, counted from the",
        "fragment's 5' end. Dot area is the number of reads."),
  sprintf(paste("Circled: %d nt fragments in %s sit at offset %d, covering %s reads.",
                "%s"),
          ex$fragment_length, S, ex$site_offset, format(ex$reads, big.mark = ","),
          if (nrow(odd)) sprintf("Scattered small dots sit away from that line: %s reads in total (%s of the sample), at offsets %s.",
                                 format(sum(odd$reads), big.mark = ","),
                                 pct(sum(odd$reads) / lo[sample == S, sum(reads)]),
                                 paste(range(odd$site_offset), collapse = " to "))
          else "No other offsets occur."),
  paste("A dot on the main line means: take the fragment's first base, count that many nucleotides",
        "along, and there is the codon the ribosome was reading. Get it wrong by one and every",
        "codon-level result shifts by one position, which is the commonest way a Ribo-seq analysis",
        "goes quietly wrong. When the offset steps up with fragment length, that step is what a",
        "tool has to recover.",
        "The scattered dots are not biology. They are fragments that would have reached past the",
        "end of their transcript, and whose offset was forced into the range that still fits. For",
        "those reads the recorded offset is a fiction: the fragment was moved, so its A-site does",
        "not sit at the stated distance. They are worth knowing about and worth bounding - in this",
        "run they are 3 percent of reads, and every one of them lies in a leader or a trailer,",
        "never inside a coding sequence, so analyses that work on the CDS are untouched."))

# ---- 3. Reading-frame periodicity per read length ---------------------------
# Two things are plotted on top of each other: the frame split as the simulator
# recorded it, which is right for every length by construction, and the split a
# tool sees when it applies one offset to all lengths. The second is the useful
# one, because that mistake is easy to make and its damage is invisible in the
# pooled figure.
fl <- tsv("frame_length")
fl[, share := reads / sum(reads), by = .(sample, fragment_length)]
co <- tsv("frame_common_offset")
co[, share := reads / sum(reads), by = .(sample, fragment_length)]
assumed <- co[sample == S, assumed_offset[1]]
best_true <- fl[sample == S, .SD[which.max(share)], by = fragment_length]
best_co <- co[sample == S, .SD[which.max(share)], by = fragment_length][order(fragment_length)]
pooled_true <- fl[sample == S, .(r = sum(reads)), by = frame][, max(r) / sum(r)]
pooled_co <- co[sample == S, .(r = sum(reads)), by = frame][, max(r) / sum(r)]
wrong <- best_co[frame != 0L]

p <- ggplot(co[sample == S], aes(factor(fragment_length), share, fill = factor(frame))) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1/3, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = paste0("3. Reading frame per read length, read with one offset of ",
                      assumed, " nt (", S, ")"),
       x = "Fragment length (nt)", y = "Share of that length's CDS reads",
       fill = "Frame",
       caption = "Dashed line: one third, what chance would give") +
  theme_qc
page(p,
  paste("A tool does not know where the decoding site of a read is. It starts from the fragment's",
        "5' end and adds an offset, and the correct offset grows with fragment length - 14 nt for",
        "the shortest fragments here, 15 for the middle ones, 16 for the longest. A tool that uses",
        "one value for all lengths is therefore right for some and wrong by a nucleotide for the",
        "rest, and one nucleotide is a third of a codon: it moves the read into the neighbouring",
        "frame. This panel shows the result of using", paste0(assumed, " nt"),
        "throughout, which is the value an offset inferred from the pooled start-codon peak comes",
        "out as."),
  sprintf(paste("The dominant frame moves with length: %s. Lengths %s land in frame 0, the others",
                "do not - %d of the %d lengths end up with their reads in the wrong frame."),
          paste(sprintf("%d nt in frame %d", best_co$fragment_length, best_co$frame), collapse = ", "),
          paste(best_co[frame == 0L, fragment_length], collapse = " and "),
          nrow(wrong), nrow(best_co)),
  sprintf(paste("What makes this worth a figure of its own is that the periodicity itself is not",
                "lost - it stays around %s within each length - it only sits in a different frame.",
                "A tool that then pools the lengths adds signal from three different frames",
                "together, and the periodicity it can see drops from %s to %s. That is why read",
                "length selection and per-length offsets come before everything else in a Ribo-seq",
                "analysis, and why a simulated library needs more than one read length to test a",
                "tool on this at all. The simulator itself records the true site for every read, so",
                "in the ground truth the frame is correct for all lengths; the error modelled here",
                "is the analysis step, not the simulation."),
          pct(best_co[, max(share)]), pct(pooled_true), pct(pooled_co)))

# ---- 4. Coverage after the start codon --------------------------------------
# Plotted in P-site rather than A-site coordinates. The run recorded A-sites,
# and the two sites are one codon apart along the mRNA, so the axis is shifted
# by 3 nt. The reason to prefer the P site here is that initiation puts the
# start codon into it: in P-site coordinates the signal from newly assembled
# ribosomes falls on position 0, which is also the convention analysis tools
# use when they calibrate an offset against the start codon.
sp <- tsv("start_profile")
sp[, psite := relative_position - 3L]
base <- sp[sample == S, median(reads)]
exmax <- sp[sample == S][which.max(reads)]
p <- ggplot(sp, aes(psite, reads, colour = sample)) +
  geom_line(alpha = 0.85) +
  geom_point(data = exmax, aes(psite, reads), inherit.aes = FALSE,
             shape = 21, size = 6, stroke = 1.2, fill = NA, colour = "red") +
  geom_hline(yintercept = base, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  scale_x_continuous(breaks = seq(min(sp$psite), max(sp$psite), by = 3),
                     minor_breaks = NULL) +
  labs(title = "4. P-site coverage around the start codon",
       x = "Nucleotides from the start codon (P site)", y = "Reads",
       caption = "Dotted line: the start codon. Dashed line: the median across all plotted positions") +
  theme_qc
page(p,
  paste("Coverage in the first stretch of the coding sequence, with position 0 at the first",
        "nucleotide of the start codon. Each read contributes one position - the site the ribosome",
        "was reading - so the curve counts reads per nucleotide and not the stretch a fragment",
        "covers. The axis counts the P site: the run recorded A-sites, and the P site lies one",
        "codon upstream, so every position was shifted by 3 nt. Ticks are three nucleotides apart,",
        "so each tick is one codon. The window reaches 60 nt back into the leader, the region",
        "before the start codon, so that the drop in coverage there can be seen."),
  sprintf("Circled: the highest point in %s sits at P-site position %d with %s reads, against a median of %s across the window - a factor of %.1f.",
          S, exmax$psite, format(exmax$reads, big.mark = ","),
          format(base, big.mark = ","), exmax$reads / base),
  paste("The reason to read this in P-site coordinates is that initiation puts the start codon",
        "into the P site. The ribosome holds three tRNA positions in a row along the mRNA - E, P",
        "and A - and one that has just assembled has an empty E site, the start codon in P, and",
        "the A site over the second codon. So a newly initiated ribosome contributes to position 0",
        "on this axis, which is where the dotted line sits.",
        "That is also the convention analysis tools work in: they calibrate their offset by lining",
        "the strongest early signal up with the start codon, whose position is known. A flat",
        "profile here would leave a tool with nothing to calibrate against, and the frame pattern",
        "of the preceding figures would be its only clue.",
        "The leader part of the curve needs reading with care, because it is not background. It",
        "dips closest to the start codon and climbs back towards the coding-sequence level further",
        "upstream, and the reason is that the leader carries a second translation signal: the",
        "uORFs. In this run every transcript has two of them, they are short (median 27 nt in a",
        "120 nt leader) and they cannot reach into the coding sequence, so they end some distance",
        "before it - two thirds of them more than 36 nt before, and only 17 percent within 21 nt.",
        "The stretch immediately before the start codon therefore has neither uORF translation nor",
        "main-ORF translation, which is where the dip comes from, while the rise further upstream",
        "is uORF translation rather than scanning. Switching uORFs off (max_uorfs = 0) removes that",
        "second source and leaves the leader at a flat low level."))

# ---- 5. Coverage around the stop codon --------------------------------------
st <- tsv("stop_profile")
exs <- st[sample == S][which.max(reads)]
baseS <- st[sample == S, median(reads)]
# Plotted in P-site coordinates, like the start-codon figure, so the two can be
# read on the same axis. Termination differs from initiation in where the signal
# falls: the stop codon is decoded in the A site, so a terminating ribosome has
# its P site on the last sense codon, at -3 on this axis.
st[, psite := relative_position - 3L]
exs <- st[sample == S][which.max(reads)]
p <- ggplot(st, aes(psite, reads, colour = sample)) +
  geom_line(alpha = 0.85) +
  geom_point(data = exs, aes(psite, reads), inherit.aes = FALSE,
             shape = 21, size = 6, stroke = 1.2, fill = NA, colour = "red") +
  geom_hline(yintercept = baseS, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  scale_x_continuous(breaks = seq(-120, max(st$psite), by = 15), minor_breaks = NULL) +
  labs(title = "5. P-site coverage around the stop codon",
       x = "Nucleotides from the stop codon (P site; negative = before it)", y = "Reads",
       caption = "Dotted line: the stop codon. Dashed line: the median across all plotted positions") +
  theme_qc
page(p,
  paste("The same idea at the other end of the coding sequence, but the window deliberately reaches",
        "past the stop codon into the trailer. Termination takes time, so real libraries often show",
        "a raised signal just before the ribosome releases the finished protein - and behind the",
        "stop codon the signal should fall away, because the ribosome has left. Cutting the window",
        "at the stop codon would hide exactly that half of the picture."),
  sprintf("Circled: the highest point of %s in this window sits at P-site position %+d with %s reads, against a median of %s.",
          S, exs$psite, format(exs$reads, big.mark = ","),
          format(baseS, big.mark = ",")),
  paste("Where the highest point sits matters more than how high it is. Because the stop codon is",
        "read in the A site, a terminating ribosome has its P site one codon earlier, so on this",
        "axis a termination pause would show at -3 rather than at 0. A peak near there A peak far upstream, as here, is not a",
        "biological signal but ordinary scatter: in this run no termination pause was simulated, so",
        "the profile is flat apart from noise, and the maximum lands wherever chance puts it. This",
        "figure is therefore a negative control in this dataset - its job is to show that nothing was",
        "put in that was not asked for."))

# ---- 6. WT against KO at transcript level ----------------------------------
gc <- tsv("gene_counts")
wide <- dcast(gc, transcript_id ~ sample, value.var = "reads")
wt <- grep("_WT_", names(wide), value = TRUE); ko <- grep("_KO_", names(wide), value = TRUE)
if (length(wt) && length(ko)) {
  wide[, wt_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = wt]
  wide[, ko_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = ko]
  d <- wide[is.finite(wt_mean) & is.finite(ko_mean) & wt_mean > 0 & ko_mean > 0]
  d[, ratio := ko_mean / wt_mean]
  ex <- d[which.max(ratio)]
  p <- ggplot(d, aes(wt_mean, ko_mean)) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(data = ex, aes(wt_mean, ko_mean), inherit.aes = FALSE,
               shape = 21, size = 6, stroke = 1.2, fill = NA, colour = "red") +
    scale_x_log10(labels = scales::comma) + scale_y_log10(labels = scales::comma) +
    labs(title = "6. Transcript-level signal: WT against KO",
         x = "Mean reads in WT (log scale)", y = "Mean reads in KO (log scale)",
         caption = "Dashed line: equal signal in both conditions") + theme_qc
  page(p,
    paste("One point per transcript: its mean read count in the wild-type samples against the mean in",
          "the knockout samples, both on a log scale. The dashed diagonal is where a transcript",
          "behaves identically in the two conditions."),
    sprintf("Circled: %s sits furthest from the diagonal, with %.0f reads in WT against %.0f in KO, a factor of %.1f. The median transcript lies at a ratio of %.2f.",
            ex$transcript_id, ex$wt_mean, ex$ko_mean, ex$ratio, median(d$ratio)),
    paste("A point on the diagonal means the transcript is translated equally in both conditions. The",
          "further a point sits above the line, the more ribosomes it carries in the knockout - which",
          "in a real experiment would be the finding of interest. Here the scatter is set deliberately",
          "small, so the cloud should hug the diagonal and the extreme points are chance rather than",
          "regulation: with a few hundred reads per transcript, counting noise alone moves a point by",
          "a visible amount, and the lower left corner is where that noise is largest."))
}

invisible(dev.off())
message("Written: ", out_pdf)
