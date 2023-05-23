################################################################################
# Extend phase sets
################################################################################

library(tidyverse)

data_dir = file.path("data_for_plots/06_extend")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

main = "figures/06_extend/main/"
supp = "figures/06_extend/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

manuscript_numbers[["06_extend"]] <- list()

# data-driven example
{

  data_example_plot_df <- read_tsv(file = file.path(data_dir, "data_example_plot_df.tsv"),
                                   show_col_types = FALSE) %>%
    mutate(recommendation = factor(recommendation, levels = c("Switch", "No Switch", "No Recommendation"), ordered = TRUE))

  x <- data_example_plot_df %>%
    filter(recommendation != "No Recommendation") %>%
    select(ps1, ps1_midpoint, ps2, recommendation)

  arrows <- data_example_plot_df %>%
    filter(recommendation != "No Recommendation") %>%
    select(ps1, ps1_midpoint, ps2, recommendation) %>%
    left_join(x, by = "ps2") %>%
    filter(ps1.x != ps1.y) %>%
    rowwise() %>%
    mutate(start_pos = min(ps1_midpoint.x, ps1_midpoint.y),
           end_pos = max(ps1_midpoint.x, ps1_midpoint.y)) %>%
    mutate(mod_switch = ((recommendation.x == "Switch") + (recommendation.y == "Switch")) %% 2) %>%
    select(start_pos, end_pos, mod_switch) %>%
    mutate(mod_switch = factor(mod_switch, levels = c(1, 0 , -1), labels = c("Switch", "No Switch", "No Recommendation"), ordered = TRUE)) %>%
    unique()

  rm(x)

  y_top_position = 2
  up_and_down = .4

  ggplot(data_example_plot_df) +
    geom_rect(aes(xmin = ps1_FirstVariantPosition,
                  xmax = ps1_LastVariantPosition,
                  ymin = y_top_position - up_and_down,
                  ymax = y_top_position + up_and_down,
                  fill = ps1)) +
    geom_rect(data = data_example_plot_df %>% select(ps2, ps2_FirstVariantPosition, ps2_LastVariantPosition) %>% unique(),
              aes(xmin = ps2_FirstVariantPosition,
                  xmax = ps2_LastVariantPosition,
                  ymin = 1 - up_and_down,
                  ymax = 1 + up_and_down,
                  fill = ps2)) +
    geom_segment(aes(x = ps1_midpoint, xend = ps2_midpoint,
                     y = y_top_position, yend = 1,
                     lty = recommendation),
                 lwd = 0.25) +
    geom_curve(data = arrows,
               aes(x = start_pos, xend = end_pos,
                   y = y_top_position, yend = y_top_position,
                   lty = mod_switch),
               curvature = -.5,
               show.legend = FALSE) +
    geom_point(aes(x = ps1_midpoint, y = y_top_position)) +
    geom_point(data = data_example_plot_df %>% select(ps2_midpoint) %>% unique(),
               aes(x = ps2_midpoint), y = 1) +
    scale_linetype_manual(values = c(1, 5, 3)) +
    scale_y_continuous(breaks = c(1, y_top_position),
                       labels = c("27522 (P)", "27522 (Rel)"),
                       limits = c(1 - up_and_down, y_top_position + up_and_down)) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_fill_brewer(palette = "Pastel1") +
    guides(fill = FALSE) + #, lty = FALSE) +
    labs(x = "chr1 Position (Mb)") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8), #, angle = 270, hjust = .5),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,.3,0.3,0), "lines"),
          legend.spacing = unit(c(0,0,0,0), "lines"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-10,-10))

  ggsave(str_c(main, "data_example.pdf"),
         width = 3.5, height = 2, useDingbats = FALSE)

  rm(arrows, up_and_down, y_top_position, data_example_plot_df)
}

# n segments per cluster
# n shared variants
{
  length_overlap_plot_df <- read_tsv(file = file.path(data_dir, "length_overlap_plot_df.tsv"),
                                     show_col_types = FALSE) %>%
    mutate(recommendation = factor(recommendation,
                                   levels = c("Switch", "No Switch", "No Recommendation"),
                                   ordered = TRUE)) %>%
    arrange(desc(recommendation))

  p <- length_overlap_plot_df %>%
    ggplot(aes(x = log10(length_overlap),
               fill = recommendation)) +
    geom_histogram(boundary = 0, binwidth = 0.25) +
    geom_vline(xintercept = log10(10^5), lty = 2) +
    labs(x = "Phase Set Overlap Length (bp, log10)", y = NULL) +
    scale_fill_viridis_d(option = "C") +
    scale_y_continuous(expand = c(0.02,0.02)) +
    scale_x_continuous(expand = c(0.02,0.02)) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"))

    ggsave(str_c(main, "length_overlap.with_legend.pdf"), p,
           useDingbats = FALSE, width = 2.25, height = 2.25)
    ggsave(str_c(main, "length_overlap.no_legend.pdf"), p + guides(fill = FALSE),
           useDingbats = FALSE, width = 2.25, height = 2.25)

    rm(p)

    length_overlap_plot_df %>%
      ggplot(aes(x = log10(n_variants_overlap + 1),
                 fill = recommendation)) +
      geom_histogram(boundary = 0, binwidth = 0.25, show.legend = FALSE) +
      geom_vline(xintercept = log10(10^1), lty = 2) +
      labs(x = "Overlapping Phased Variants (log10)", y = NULL) +
      scale_fill_viridis_d(option = "C") +
      scale_y_continuous(expand = c(0.02,0.02)) +
      scale_x_continuous(expand = c(0.02,0.02)) +
      theme_bw() +
      theme(axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_line(size = 0.5),
            strip.background = element_blank(),
            strip.text = element_text(size = 8),
            plot.margin = unit(c(0,0,0,0), "lines"))

    ggsave(str_c(main, "variants_overlap.pdf"),
           useDingbats = FALSE, width = 2.25, height = 2.25)

    rm(length_overlap_plot_df)
}

# example of extended phase sets
{

  extension_landscape_plot_df <- read_tsv(file = file.path(data_dir, "extension_landscape_plot_df.tsv"),
                                          show_col_types = FALSE) %>%
    mutate(recommendation = factor(recommendation,
                                   levels = c("Switch", "No Switch", "No Recommendation"),
                                   ordered = TRUE)) %>%
    filter((sample == "27522_1" & extended_by == "27522_2" & ps1_Chromosome == "chr1")) %>%
    arrange(length_after) %>%
    mutate(ps2 = fct_reorder(factor(ps2), length_after))

  ggplot(data = extension_landscape_plot_df, aes(y = fct_reorder(ps2, length_after))) +
    geom_segment(aes(x = left_aligned_first_variant_position/1e6,
                     xend = left_aligned_last_variant_position/1e6,
                     yend = ps2,
                     color = recommendation),
                 lwd = 2,
                 show.legend = FALSE) +
    geom_segment(aes(x = left_aligned_first_variant_position/1e6,
                     xend = left_aligned_first_variant_position/1e6 + 0.1,
                     yend = ps2,
                     color = recommendation),
                 lwd = 4,
                 show.legend = FALSE) +
    geom_segment(aes(x = left_aligned_last_variant_position/1e6 - 0.1,
                     xend = left_aligned_last_variant_position/1e6,
                     yend = ps2,
                     color = recommendation),
                 lwd = 4,
                 show.legend = FALSE) +
    labs(x = "Genomic Length and Relative Position (Mb)",
         y = "Groups of Extendable Phase Sets") +
    scale_color_viridis_d(option = "C", drop = F) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(main, "extension_landscape.pdf"),
         useDingbats = FALSE, width = 4.75, height = 4.75)

  plot_df_dist_after_plot_df <- read_tsv(file = file.path(data_dir, "plot_df_dist_after_plot_df.tsv"),
                                         show_col_types = FALSE)

  plot_df_dist_before_plot_df <- read_tsv(file = file.path(data_dir, "plot_df_dist_before_plot_df.tsv"),
                                          show_col_types = FALSE)

  ggplot() +
    geom_violin(data = plot_df_dist_before_plot_df, aes(x = "Before",
                                                y = log10(ps1_LastVariantPosition - ps1_FirstVariantPosition)), draw_quantiles = 0.5) +
    geom_violin(data = plot_df_dist_after_plot_df, aes(x = "Extended",
                                               y = log10(length_after)), draw_quantiles = 0.5) +
    expand_limits(y = c(0)) +
    labs(y = "Phase Set Length (bp, log10)", x = NULL) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          #plot.margin = unit(c(0,0,0,0), "lines"),
          strip.text = element_text(size = 8))

  ggsave(str_c(main, "extension_distribution.pdf"),
         useDingbats = FALSE, width = 1.75, height = 3.25)

  rm(plot_df_dist_before_plot_df, plot_df_dist_after_plot_df)
}

rm(main, supp, data_dir, extension_landscape_plot_df)
