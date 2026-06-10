g_end_dropout_label <- function(drop_fraction) {
  paste0("drop_G_", as.integer(round(drop_fraction * 100)))
}

g_end_dropout_awk_program <- function() {
  paste(
    'BEGIN { srand(seed); total = 0; g_end = 0; removed = 0 }',
    '/^@/ { print; next }',
    '{',
    '  total++',
    '  if ($10 ~ /G$/) {',
    '    g_end++',
    '    if (rand() < drop) { removed++; next }',
    '  }',
    '  print',
    '}',
    'END {',
    '  printf("total\\t%d\\n", total) > stats',
    '  printf("g_end\\t%d\\n", g_end) >> stats',
    '  printf("removed\\t%d\\n", removed) >> stats',
    '  printf("kept\\t%d\\n", total - removed) >> stats',
    '}',
    sep = "\n"
  )
}

rewrite_experiment_base <- function(exp_lines, input_base, output_base) {
  gsub(input_base, output_base, exp_lines, fixed = TRUE)
}
