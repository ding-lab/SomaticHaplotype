################################################################################
# Use ancestry information
################################################################################

main = "figures/07_ancestry/main/"
supp = "figures/07_ancestry/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

manuscript_numbers[["07_ancestry"]] <- list()

manuscript_numbers[["07_ancestry"]][["shared_event_between_NA12878_NA10851"]] <- segment_overlap_tbl %>%
  filter(Sample2_ID == "NA10851")

positions_list <- segment_overlap_tbl %>%
  filter(Chromosome == "chr18", IBD_segment_index == 8) %>%
  pull(Overlapping_variant_positions) %>%
  str_split(pattern =  ",")
H1_allele_list <- segment_overlap_tbl %>%
  filter(Chromosome == "chr18", IBD_segment_index == 8) %>%
  pull(Phase_set_H1_alleles) %>%
  str_split(pattern =  ",")
H2_allele_list <- segment_overlap_tbl %>%
  filter(Chromosome == "chr18", IBD_segment_index == 8) %>%
  pull(Phase_set_H2_alleles) %>%
  str_split(pattern =  ",")
IBD_allele_list <- segment_overlap_tbl %>%
  filter(Chromosome == "chr18", IBD_segment_index == 8) %>%
  pull(IBD_alleles) %>%
  str_split(pattern =  ",")

n_alleles_ps1 <- length(positions_list[[1]])
n_alleles_ps2 <- length(positions_list[[4]])

x_limits <- c(min(as.numeric(positions_list[[1]])),
              max(as.numeric(positions_list[[4]])))

ibd <- segment_overlap_tbl %>%
  filter(Chromosome == "chr18", IBD_segment_index == 8) %>%
  select(IBD_start_position, IBD_end_position) %>%
  mutate(my_row = "IBD Segment") %>% unique()

ps <- segment_overlap_tbl %>%
  filter(Chromosome == "chr18", IBD_segment_index == 8) %>%
  select(Phase_set_key, Phase_set_start_position, Phase_set_end_position) %>%
  mutate(my_row = "Phase Sets") %>%
  rowwise() %>%
  mutate(label_position = case_when(Phase_set_start_position < x_limits[1] ~ median(c(x_limits[1], Phase_set_end_position)),
                                    Phase_set_end_position > x_limits[2] ~ median(c(x_limits[2], Phase_set_start_position)),
                                    TRUE ~ median(c(Phase_set_start_position, Phase_set_end_position))))

tibble(my_row = c(rep("Phase Sets", n_alleles_ps1), rep("Phase Sets", n_alleles_ps2)),
       position = c(as.numeric(positions_list[[1]]), as.numeric(positions_list[[4]])),
       matches_h1 = c(IBD_allele_list[[1]] == H1_allele_list[[1]],
                      IBD_allele_list[[4]] == H1_allele_list[[4]]),
       matches_h2 = c(IBD_allele_list[[1]] == H2_allele_list[[1]],
                      IBD_allele_list[[4]] == H2_allele_list[[4]]),
       match_category = as.numeric(matches_h1) + 2*as.numeric(matches_h2)) %>%
  mutate(match_category = factor(match_category,
                                 levels = c(0, -1, -2, 1, 2),
                                 labels = c("Neither", "X", "Y", "H1", "H2"),
                                 ordered = TRUE)) %>%
  bind_rows(ps) %>%
  bind_rows(ibd) %>%
  ggplot(aes(x = position/1e6, y = my_row)) +
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
        legend.box.margin = margin(-10,-10,-10,-10)) +
  ggsave(str_c(main, "ibd_overlap.pdf"),
         width = 7.25, height = 1.5, useDingbats = FALSE)

rm(main, supp, positions_list, H1_allele_list, H2_allele_list, IBD_allele_list,
  n_alleles_ps1, n_alleles_ps2, x_limits, ibd, ps)
