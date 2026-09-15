library(tidyverse)

# ------------------------------------------------------------------
# Read coverage files
# ------------------------------------------------------------------

files <- list.files(
  "COVERAGE",
  pattern = "\\.coverage\\.txt$",
  full.names = TRUE
)

coverage <- map_dfr(files, function(f) {

  x <- read.table(
    f,
    comment.char = "#",
    header = FALSE,
    stringsAsFactors = FALSE
  )

  tibble(
    file = f,
    coverage = x[[6]]
  )
})


# ------------------------------------------------------------------
# Extract metadata from filenames
# ------------------------------------------------------------------

coverage <- coverage %>%
  mutate(
    filename = basename(file),

    # Tissue
    # "traheal" is the spelling used in the supplied filenames
    tissue = case_when(
      str_detect(filename, regex("bronchial", ignore_case = TRUE)) ~ "Bronchial",
      str_detect(filename, regex("tra?cheal|traheal", ignore_case = TRUE)) ~ "Tracheal",
      TRUE ~ NA_character_
    ),

    # Timepoint
    timepoint = case_when(
      str_detect(filename, "6hpi") ~ "6 hpi",
      str_detect(filename, "24hpi") ~ "24 hpi",
      TRUE ~ NA_character_
    ),

    # Infection status
    infection = case_when(
      str_detect(filename, regex("mock", ignore_case = TRUE)) ~ "Mock",
      TRUE ~ "Infected"
    ),

    # Animal
    animal = str_extract(filename, "Animal-[0-9]+") %>%
      str_remove("Animal-")
  )


# ------------------------------------------------------------------
# Check for failed metadata extraction
# ------------------------------------------------------------------

if (any(is.na(coverage$tissue))) {
  warning(
    "Tissue could not be identified for: ",
    paste(
      coverage$filename[is.na(coverage$tissue)],
      collapse = ", "
    )
  )
}

if (any(is.na(coverage$timepoint))) {
  warning(
    "Timepoint could not be identified for: ",
    paste(
      coverage$filename[is.na(coverage$timepoint)],
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------
# Define the 8 groups
# ------------------------------------------------------------------

coverage <- coverage %>%
  mutate(
    group = factor(
      paste(timepoint, tissue, infection, sep = " - "),
      levels = c(
        "6 hpi - Bronchial - Infected",
        "6 hpi - Bronchial - Mock",
        "6 hpi - Tracheal - Infected",
        "6 hpi - Tracheal - Mock",
        "24 hpi - Bronchial - Infected",
        "24 hpi - Bronchial - Mock",
        "24 hpi - Tracheal - Infected",
        "24 hpi - Tracheal - Mock"
      )
    )
  )


# ------------------------------------------------------------------
# Print sample metadata
# ------------------------------------------------------------------

print(
  coverage %>%
    select(
      animal,
      tissue,
      timepoint,
      infection,
      coverage,
      group
    ) %>%
    arrange(
      timepoint,
      tissue,
      infection,
      animal
    )
)


# ------------------------------------------------------------------
# Count samples in each group
# ------------------------------------------------------------------

print(
  coverage %>%
    count(group, .drop = FALSE)
)


# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------

p <- ggplot(
  coverage,
  aes(
    x = group,
    y = coverage
  )
) +

  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.15,
    height = 0,
    size = 2.5,
    alpha = 0.8
  ) +

  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.05))
  ) +

  labs(
    x = NULL,
    y = "Genome coverage (%)"
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

png(
  filename = "coverage.png",
  width = 1000,
  height =1000,
  res = 300 )

print(p)

dev.off()

cat("\nPlot written to: coverage.png\n")# ==================================================================
# PAIRED ANALYSIS: 6 vs 24 hpi
# Ignore mock samples
# Animals are paired between timepoints
# ==================================================================


# ------------------------------------------------------------------
# Keep infected samples only
# ------------------------------------------------------------------

paired <- coverage %>%
  filter(
    infection == "Infected",
    timepoint %in% c("6 hpi", "24 hpi")
  )


# ------------------------------------------------------------------
# Check that animals have both timepoints
# ------------------------------------------------------------------

pair_check <- paired %>%
  count(tissue, animal, timepoint) %>%
  complete(
    tissue,
    animal,
    timepoint = c("6 hpi", "24 hpi"),
    fill = list(n = 0)
  )

print(pair_check)


# ------------------------------------------------------------------
# Keep only complete 6/24 hpi pairs
# ------------------------------------------------------------------

complete_pairs <- pair_check %>%
  group_by(tissue, animal) %>%
  filter(all(n > 0)) %>%
  distinct(tissue, animal)


paired <- paired %>%
  semi_join(
    complete_pairs,
    by = c("tissue", "animal")
  )


# ------------------------------------------------------------------
# Set timepoint order
# ------------------------------------------------------------------

paired <- paired %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("6 hpi", "24 hpi")
    ),
    tissue = factor(
      tissue,
      levels = c("Bronchial", "Tracheal")
    )
  )


# ------------------------------------------------------------------
# Paired 6 -> 24 hpi plot
# ------------------------------------------------------------------

p_paired <- ggplot(
  paired,
  aes(
    x = timepoint,
    y = coverage,
    group = animal
  )
) +

  # Individual animal trajectories
  geom_line(
    alpha = 0.5,
    linewidth = 0.7
  ) +

  # Individual samples
  geom_point(
    size = 3
  ) +

  # Distribution at each timepoint
  geom_boxplot(
    aes(group = timepoint),
    width = 0.25,
    alpha = 0.25,
    outlier.shape = NA
  ) +

  facet_wrap(
    ~ tissue,
    nrow = 1
  ) +

  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.05))
  ) +

  labs(
    x = NULL,
    y = "Genome coverage (%)"
  ) +

  theme_classic(
    base_size = 13
  )

png(
  filename = "coverage_6_vs_24hpi_paired.png",
  width = 1200,
  height = 600,
  res = 300
)

print(p_paired)

dev.off()


# ------------------------------------------------------------------
# Calculate change in coverage for each animal
# ------------------------------------------------------------------

coverage_change <- paired %>%
  select(
    animal,
    tissue,
    timepoint,
    coverage
  ) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = coverage
  ) %>%
  mutate(
    change = `24 hpi` - `6 hpi`,
    fold_change = `24 hpi` / `6 hpi`
  )


# ------------------------------------------------------------------
# Print paired changes
# ------------------------------------------------------------------

print(
  coverage_change %>%
    arrange(tissue, animal)
)


# ------------------------------------------------------------------
# Plot change in coverage
# ------------------------------------------------------------------

p_change <- ggplot(
  coverage_change,
  aes(
    x = tissue,
    y = change
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +

  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.12,
    size = 3,
    alpha = 0.8
  ) +

  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  ) +

  labs(
    x = NULL,
    y = "Change in genome coverage (%)\n24 hpi - 6 hpi"
  ) +

  theme_classic(
    base_size = 13
  )

png(
  filename = "coverage_change_24_minus_6hpi.png",
  width = 1200,
  height = 800,
  res = 300
)

print(p_change)

dev.off()