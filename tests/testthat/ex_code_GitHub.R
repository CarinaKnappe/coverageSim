#library(coverageSim)
#library(ORFik)

## Simple example
# 6 genes on 6 chromosomes (1 active uORF each)
simGenome6 <- simGenome(n = 6, max_uorfs = 1)
# Simulate Ribo-seq only
gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
                                  libtypes = "RFP", print_statistics = FALSE)
region_count_table <- simCountTablesRegions(gene_count_table,
                                            regionsToSample = c("leader", "cds", "trailer"))
df <- simNGScoverage(simGenome6, region_count_table)


# testing copilot : Comment of what is todo: newline: accept suggestion by hitting the tab-key, or by clicking on the suggestion.
# If you want to see more suggestions, press the tab key multiple times. You can also ask for a new suggestion by pressing Ctrl+Enter.

# Create a scatter plot of mpg vs hp from the mtcars dataset
df <- mtcars

# Simulated coverage data
df <- data.frame(
  position = 1:100,             # 100 positions
  coverage = runif(100, 10, 50) # 100 random coverage values between 10 and 50
)

# Base R scatter plot
plot(df$position, df$coverage,
     main = "Simulated Coverage",
     xlab = "Position",
     ylab = "Coverage")

plot(df$position, df$coverage, main = "Simulated Coverage", xlab = "Position", ylab = "Coverage")

# Create a ggplot of mpg vs hp colored by cyl
ggplot(df, aes(x = position, y = coverage)) +
  geom_point() +
  labs(title = "Simulated Coverage", x = "Position", y = "Coverage") +
  theme_minimal()
