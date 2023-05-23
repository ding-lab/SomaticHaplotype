################################################################################
# Use ancestry information
################################################################################

data_dir = file.path("data_for_plots/07_ancestry")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

main = "figures/07_ancestry/main/"
supp = "figures/07_ancestry/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

ibd_overlap_plot_df <- read_tsv(file = file.path(data_dir, "ibd_overlap_plot_df.tsv"),
                                show_col_types = FALSE) %>%
  mutate(match_category = factor(match_category,
                                 levels = c("H1", "H2"),
                                 ordered = TRUE))

x_limits <- c(min(ibd_overlap_plot_df$position, na.rm = TRUE),
              max(ibd_overlap_plot_df$position, na.rm = TRUE))

ggplot(data = ibd_overlap_plot_df,
       aes(x = position/1e6, y = my_row)) +
  geom_segment(aes(x = Phase_set_start_position/1e6, xend = Phase_set_end_position/1e6,
                   y = my_row, yend = my_row)) +
  geom_segment(aes(x = IBD_start_position/1e6, xend = IBD_end_position/1e6,
                   y = my_row, yend = my_row)) +
  geom_jitter(aes(color = match_category), width = 0, height = 0.1, shape = 16, alpha = 0.25) +
  geom_jitter(aes(y = "IBD Segment"), width = 0, height = 0.1, shape = 16, alpha = 0.25) +
  geom_label(aes(x = label_position/1e6, y = "Phase Sets", label = Phase_set_key),
             nudge_y = .5, size = 8/ggplot2:::.pt) +
  coord_cartesian(xlim = x_limits/1e6) +
  scale_colour_manual(values = c("#ae8dc1", "#7fbf7b")) +
  labs(x = "chr18 Position (Mb)", y = NULL, color = "Matching Allele Haplotype") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_line(size = 0.5),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10))

ggsave(str_c(main, "ibd_overlap.pdf"),
       width = 7.25, height = 1.5, useDingbats = FALSE)

rm(main, supp, x_limits, ibd_overlap_plot_df)
