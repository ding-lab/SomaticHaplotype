################################################################################
# Overview figure
# Data available and QC
################################################################################

library(tidyverse)

data_dir = file.path("data_for_plots/01_overview")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

main = "figures/01_overview/main/"
supp = "figures/01_overview/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# Data available
{
  samples_plot_df <- read_tsv(file.path(data_dir, "samples_plot_df.tsv"),
                              show_col_types = FALSE)

  ggplot(samples_plot_df, aes(x = timepoint, y = patient)) +
    geom_point(aes(color = my_color_100), shape = 16, size = 3, show.legend = FALSE) +
    geom_point(data = samples_plot_df %>% filter(sorted), shape = 16, color = "#ffffff") +
    scale_color_identity() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.spacing = unit(c(0,0,0,0), "null"),
          plot.margin = unit(c(0,0,0,0), "null"),
          strip.text = element_text(size = 8))

  ggsave(str_c(main, "samples.pdf"),
         width = 2.75,
         height = 2.5,
         useDingbats = FALSE)

  rm(samples_plot_df)
}

# Data QC (all metrics)
{
  qc_metrics_plot_df <- read_tsv(file = file.path(data_dir, "qc_metrics_plot_df.tsv"),
                                 show_col_types = FALSE)

  ggplot(qc_metrics_plot_df %>% filter(!str_detect(patient, "NA1")),
         aes(x = normal_sample, y = result*data_multiplier)) +
    geom_violin(show.legend = FALSE, draw_quantiles = .5) +
    geom_jitter(aes(color = my_color_100, shape = my_shape),
                height = 0, width = 0.25, show.legend = FALSE) +
    geom_point(data = qc_metrics_plot_df %>% filter(str_detect(patient, "NA1")),
               aes(color = my_color_100, shape = my_shape),
               show.legend = FALSE) +
    geom_label(data = qc_metrics_plot_df %>%
                 select(category, data_multiplier_label, normal_sample) %>%
                 filter(normal_sample == "Tumor") %>% unique(),
               aes(label = data_multiplier_label),
               y = -Inf, color = "#000000") +
    facet_wrap(~category, ncol = 4, scales = "free_y") +
    scale_color_identity() +
    scale_shape_identity() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 8)
    )

  ggsave(str_c(supp, "qc_metrics.pdf"),
         height = 7.5,
         width = 7.5,
         useDingbats = FALSE)

  rm(qc_metrics_plot_df)
}

# Data QC (top metrics)
{

  highlight_columns <- sort(c("n50_linked_reads_per_molecule",
                              "n50_phase_block",
                              "molecule_length_mean"))

  top_qc_metrics_plot_df <- read_tsv(file = file.path(data_dir, "top_qc_metrics_plot_df.tsv"),
                                     na = c("", "NA", "NaN"),
                                     show_col_types = FALSE) %>%
    mutate(category = factor(category,
                             levels = highlight_columns,
                             labels = c("Molecule Length\\n(mean, Kb)",
                                        "Linked-reads per\\nmolecule (N50)",
                                        "Phase Set Length\\n(N50, Mb)")))

  ggplot(top_qc_metrics_plot_df %>%
           filter(!str_detect(patient, "NA1")),
         aes(x = normal_sample, y = result*data_multiplier)) +
    geom_violin(show.legend = FALSE, draw_quantiles = .5) +
    geom_jitter(aes(color = my_color_100, shape = my_shape),
                height = 0, width = 0.25, show.legend = FALSE) +
    geom_point(data = top_qc_metrics_plot_df %>%
                 filter(str_detect(patient, "NA1")),
               aes(color = my_color_100, shape = my_shape),
               show.legend = FALSE) +
    expand_limits(y = 0) +
    facet_wrap(~category, nrow = 1, scales = "free_y") +
    scale_color_identity() +
    scale_shape_identity() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(main, "top_qc_metrics.pdf"),
         useDingbats = FALSE,
         height = 2.5,
         width = 4.25)

  rm(top_qc_metrics_plot_df, highlight_columns)

}

rm(main, supp, data_dir)
