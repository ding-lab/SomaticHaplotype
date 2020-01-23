################################################################################
# Extend phase sets
################################################################################

main = "figures/06_extend/main/"
supp = "figures/06_extend/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

extend_ps_no_skin <- extend_phase_sets_tbl %>%
  filter(!is.na(cluster_base_phase_set)) %>%
  filter(!str_detect(sample, "skin"), !str_detect(extended_by, "skin"))
extend_stats_no_skin <- extend_stats_tbl %>%
  mutate(pct_switch_hilo = n_variants_flip_to_match/n_variants_overlap) %>%
  filter(!str_detect(sample, "skin"), !str_detect(extended_by, "skin"))

# data-driven example
{
  plot_df <- extend_stats_no_skin %>%
    #filter(ps2 %in% c("chr1:109390814", "chr1:112459045", "chr1:115991285"),
    #filter(ps2 %in% c("chr1:112459045", "chr1:115991285"),
    filter(ps2 %in% c("chr1:112459045"),
           sample == "27522_3", extended_by == "27522_1") %>%
    filter(ps1_LastVariantPosition - ps1_FirstVariantPosition > 1e3,
           ps2_LastVariantPosition - ps2_FirstVariantPosition > 1e3) %>%
    select(ps1, ps1_FirstVariantPosition, ps1_LastVariantPosition,
           ps2, ps2_FirstVariantPosition, ps2_LastVariantPosition,
           min_overlap_position, max_overlap_position,
           n_variants_overlap, n_variants_flip_to_match,
           recommendation, graph_weight) %>%
    mutate(ps1_FirstVariantPosition = ps1_FirstVariantPosition/1e6,
           ps2_FirstVariantPosition = ps2_FirstVariantPosition/1e6,
           ps1_LastVariantPosition = ps1_LastVariantPosition/1e6,
           ps2_LastVariantPosition = ps2_LastVariantPosition/1e6,
           min_overlap_position = min_overlap_position/1e6,
           max_overlap_position = max_overlap_position/1e6) %>%
    mutate(ps1_midpoint = ps1_FirstVariantPosition + (ps1_LastVariantPosition - ps1_FirstVariantPosition)/2,
           ps2_midpoint = ps2_FirstVariantPosition + (ps2_LastVariantPosition - ps2_FirstVariantPosition)/2) %>%
    mutate(recommendation = factor(recommendation, levels = c("Switch", "No Switch", "No Recommendation"), ordered = TRUE))

  x <- plot_df %>%
    filter(recommendation != "No Recommendation") %>%
    select(ps1, ps1_midpoint, ps2, recommendation)

  arrows <- plot_df %>%
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
  ggplot(plot_df) +
    geom_rect(aes(xmin = ps1_FirstVariantPosition,
                  xmax = ps1_LastVariantPosition,
                  ymin = y_top_position - up_and_down,
                  ymax = y_top_position + up_and_down,
                  fill = ps1)) +
    geom_rect(data = plot_df %>% select(ps2, ps2_FirstVariantPosition, ps2_LastVariantPosition) %>% unique(),
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
    geom_point(data = plot_df %>% select(ps2_midpoint) %>% unique(),
               aes(x = ps2_midpoint), y = 1) +
    scale_linetype_manual(values = c(1, 5, 3)) +
    scale_y_continuous(breaks = c(1, y_top_position),
                       labels = c("27522_1", "27522_3"),
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
          #plot.margin = unit(c(0,0,0,0), "lines"),
          legend.spacing = unit(c(0,0,0,0), "lines"),
          legend.position = "bottom",
          #legend.title = element_text(size = 8),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-10,-10)) +
    ggsave(str_c(main, "data_example.pdf"),
           width = 3.5, height = 2, useDingbats = FALSE)

  rm(arrows, plot_df, up_and_down, y_top_position)
}

# n segments per cluster
# n shared variants
{
  plot_df <- extend_stats_no_skin %>%
    mutate(recommendation = case_when(pct_switch_hilo < 0.95 & recommendation == "Switch" ~ "No Recommendation",
                                      pct_switch_hilo > 0.05 & recommendation == "No Switch" ~ "No Recommendation",
                                      TRUE ~ recommendation)) %>%
    mutate(recommendation = factor(recommendation,
                                   levels = c("Switch", "No Switch", "No Recommendation"),
                                   ordered = TRUE)) %>%
    arrange(desc(recommendation))

  plot_df %>%
    ggplot(aes(x = log10(length_overlap),
               fill = recommendation)) +
    geom_histogram(boundary = 0, binwidth = 0.25, show.legend = FALSE) +
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
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "length_overlap.pdf"),
           useDingbats = FALSE, width = 2.25, height = 2.25)

  plot_df %>%
    ggplot(aes(x = log10(n_variants_overlap + 1),
               fill = recommendation)) +
    geom_histogram(boundary = 0, binwidth = 0.25, show.legend = FALSE) +
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
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "variants_overlap.pdf"),
           useDingbats = FALSE, width = 2.25, height = 2.25)

  rm(plot_df)
}

# example of extended phase sets
{

  rm_no_rec <- extend_stats_no_skin %>%
    mutate(recommendation = case_when(pct_switch_hilo < 0.95 & recommendation == "Switch" ~ "No Recommendation",
                                      pct_switch_hilo > 0.05 & recommendation == "No Switch" ~ "No Recommendation",
                                      TRUE ~ recommendation)) %>%
    mutate(recommendation = factor(recommendation,
                                   levels = c("Switch", "No Switch", "No Recommendation"),
                                   ordered = TRUE)) %>%
    filter(recommendation != "No Recommendation")

  plot_df <- rm_no_rec %>%
    group_by(sample, extended_by, ps2) %>%
    summarize(min_position = min(ps1_FirstVariantPosition),
              max_position = max(ps1_LastVariantPosition),
              length_after = max_position - min_position,
              count = n()) %>%
    filter(count > 1) %>%
    ungroup() %>%
    left_join(rm_no_rec,
              by = c("sample", "extended_by", "ps2")) %>%
    ungroup() %>%
    mutate(ps_midpoint = ps1_FirstVariantPosition + (ps1_LastVariantPosition - ps1_FirstVariantPosition)/2,
           left_aligned_ps_midpoint = ps_midpoint - min_position,
           left_aligned_first_variant_position = ps1_FirstVariantPosition - min_position,
           left_aligned_last_variant_position = ps1_LastVariantPosition - min_position) %>%
    filter((sample == "27522_1" & extended_by == "27522_2" & ps1_Chromosome == "chr1") | (sample == "58408_2" & extended_by == "58408_1" & ps1_Chromosome == "chr1")) %>%
    arrange(length_after) %>%
    mutate(ps2 = fct_reorder(factor(ps2), length_after))

  ggplot(data = plot_df, aes(y = fct_reorder(ps2, length_after))) +
    geom_segment(aes(x = left_aligned_first_variant_position/1e6,
                     xend = left_aligned_last_variant_position/1e6,
                     yend = ps2,
                     color = recommendation),
                 lwd = 2,
                 show.legend = FALSE) +
    geom_point(aes(x = left_aligned_ps_midpoint/1e6, color = recommendation),
               size = 3,
               show.legend = FALSE) +
    labs(x = "Relative Genomic Length (Mb)",
         y = "Group of Extended Phase Sets") +
    facet_wrap(~sample, scales = "free_y") +
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
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "extension_landscape.pdf"),
           useDingbats = FALSE, width = 4.75, height = 4.75)

  rm(plot_df, rm_no_rec)
}

rm(extend_stats_no_skin, extend_ps_no_skin)
