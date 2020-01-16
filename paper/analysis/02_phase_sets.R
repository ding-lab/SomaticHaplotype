################################################################################
# Phase sets figure
# Distribution of phase set lengths
# Relationship of phase set lengths to genomic location
# Heterozygotes/Coverage in unphased regions
################################################################################

main = "figures/02_phase_sets/main/"
supp = "figures/02_phase_sets/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# Phase set N50 by chromosome

{
  phase_set_summary_tbl %>%
    filter(timepoint != "Normal") %>%
    ggplot(aes(x = chromosome, y = N50_broad/1e6)) +
    geom_violin(draw_quantiles = 0.5) +
    geom_jitter(aes(color = my_color_100),
                height = 0, width = 0.25, shape = 16, show.legend = FALSE) +
    labs(x = NULL, y = "Phase Set Length (N50, Mb)") +
    scale_color_identity() +
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
    ggsave(str_c(main, "phase_set_by_chromosome.pdf"),
           useDingbats = FALSE, width = 7.25, height = 1.5)
}

# Phase set N50 by sample

{
  phase_set_summary_tbl %>%
    filter(timepoint != "Normal") %>%
    ggplot(aes(x = sample, y = N50_broad/1e6)) +
    geom_violin(draw_quantiles = 0.5) +
    geom_jitter(aes(color = my_color_100),
                height = 0, width = 0.25, shape = 16, show.legend = FALSE) +
    labs(x = NULL, y = "Phase Set Length (N50, Mb)") +
    scale_color_identity() +
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
    ggsave(str_c(main, "phase_set_by_sample.pdf"),
           useDingbats = FALSE, width = 7.25, height = 1.5)
}

# Phase set length by deletion
{
  phase_sets_tbl %>%
    filter(sample == "27522_1", length_variants > 1e3) %>%
    #mutate(deleted = case_when(chromosome == %in% c("chr13", "chr22") ~ "chr13 or chr22\n(Deletion)",
    #                           TRUE ~ "Other\nchromosomes")) %>%
    mutate(deleted = case_when(chromosome == "chr13" ~ "chr13",
                               chromosome == "chr22" ~ "chr22",
                               TRUE ~ "Others")) %>%
    ggplot(aes(x = deleted, y = length_variants/1e6)) +
    geom_violin(draw_quantiles = 0.5) +
    geom_jitter(aes(color = my_color_100), alpha = 0.25, shape = 16,
                height = 0, width = 0.25, show.legend = FALSE) +
    scale_color_identity() +
    labs(x = NULL, y = "Phase Set Length (Mb)") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "phase_set_deletion.pdf"),
           useDingbats = FALSE, width = 2, height = 1.5)
}

# Phase sets genome coverage
{
  phase_sets_tbl %>% filter(!is.na(length_variants),
                            length_variants > 0,
                            timepoint != "Normal") %>%
    mutate(length_group = plyr::round_any(length_variants,
                                          accuracy = 1e6,
                                          f = ceiling)/1e6) %>%
    mutate(length_label = str_c(length_group - 1, "-", length_group, " Mb")) %>%
    group_by(length_group) %>%
    summarize(total_length = sum(as.numeric(length_variants)),
              count = n()) %>%
    ungroup() %>%
    mutate(length_label = str_c(length_group - 1, "-", length_group)) %>%
    mutate(length_factor = factor(length_label, levels = length_label, ordered = TRUE)) %>%
    ggplot(aes(x = length_factor, y = total_length/1e9)) +
    geom_col(fill = "#bdbdbd") +
    geom_text(aes(label = str_c("", count)), angle = -90,
              hjust = 1, vjust = 0.5,  nudge_y = 0.05, size = 2) +
    scale_fill_identity() +
    labs(x = "Phase Set Length (Mb)", y = "Phase Sets Genome Coverage (Gb)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")
    ) +
    ggsave(str_c(main, "phase_set_genome_coverage.pdf"),
           useDingbats = FALSE, width = 7.25, height = 1.75)
}

# Phase sets by sample
{
  phase_sets_tbl_short <- phase_sets_tbl %>%
    filter(timepoint != "Normal", length_variants > 0, length_variants <= 1e3) %>%
    mutate(phase_set_color = "#636363")

  phase_sets_tbl_long <- phase_sets_tbl %>%
    filter(timepoint != "Normal", length_variants > 1e3) %>%
    mutate(event_number = row_number() %% 2) %>%
    mutate(phase_set_color = case_when(event_number == 0 ~ my_color_100,
                                       TRUE ~ my_color_50))

  plot_df <- phase_sets_tbl_long %>% bind_rows(phase_sets_tbl_short) %>%
    mutate(chr_num = as.numeric(chromosome)) %>%
    mutate(sample_num = as.numeric(factor(sample)))

  sample_names <- plot_df %>% filter(chromosome == "chr13") %>%
    group_by(patient, sample) %>%
    summarize(mean_length = mean(length_variants)) %>%
    group_by(patient) %>%
    summarize(max_sample = sample[which.max(mean_length)]) %>%
    pull(max_sample)
  sample_names <- sort(c(sample_names, "27522_1"))
  n_samples <- length(sample_names)

  ggplot(plot_df %>% filter(chromosome == "chr13",
                            sample %in% sample_names) %>%
           mutate(sample_num = as.numeric(factor(sample))),
         aes(xmin = start/1e6, xmax = end/1e6,
             ymin = sample_num - 0.25, ymax = sample_num + 0.25)) +
    geom_rect(aes(fill = "#ffffff"), show.legend = FALSE) +
    geom_rect(aes(fill = phase_set_color), show.legend = FALSE) +
    scale_y_continuous(breaks = seq(1, n_samples),
                       expand = c(0.02, 0),
                       limits = c(1 - 0.25, n_samples + 0.25),
                       labels = sample_names) +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_identity() +
    labs(x = "chr13 Position (Mb)") +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "chr13_phase_sets.pdf"),
                 height = 2.65, width = 5, useDingbats = FALSE)

  rm(phase_sets_tbl_short, phase_sets_tbl_long, plot_df, sample_names, n_samples)
}

# phase set length vs. n heterozygous SNPs
{
  phase_sets_tbl %>%
    filter(n_variants_total > 0,
           timepoint != "Normal") %>%
    ggplot(aes(x = length_variants/1e6, y = n_variants_total/1e3)) +
    geom_point(aes(color = my_color_100), shape = 16, show.legend = FALSE) +
    geom_smooth(method = "lm") +
    expand_limits(x = 0) +
    scale_color_identity() +
    labs(x = "Phase Set Length (Mb)", y = "Number of Phased Heterozygotes (thousands)") +
    facet_wrap( ~ sample, ncol = 4) + #, scales = "free") +
    theme_bw() +
    theme(panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          #axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "phase_set_length_vs_variants.pdf"),
           width = 7.5, height = 7.5, useDingbats = FALSE)

  phase_sets_tbl %>%
    filter(n_variants_total > 0,
           sample == "27522_1") %>%
    ggplot(aes(x = length_variants/1e6, y = n_variants_total/1e3)) +
    geom_point(aes(color = my_color_100), shape = 16, show.legend = FALSE) +
    geom_smooth(method = "lm") +
    expand_limits(x = 0) +
    scale_color_identity() +
    labs(x = "Phase Set Length (Mb)", y = "Number of Phased Heterozygotes (thousands)") +
    facet_wrap( ~ chromosome, ncol = 4) + #, scales = "free") +
    theme_bw() +
    theme(panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          #axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "phase_set_length_vs_variants.27522_1.pdf"),
           width = 7.5, height = 7.5, useDingbats = FALSE)
}

# TODO draw phase sets for all samples
{
  # For each sample, draw phase blocks on chromosomes
  library(tidyverse)

  combined <- read.delim("data_for_plotting/combined.ps.tsv")
  combined <- combined %>% mutate(chr = factor(chr, levels = str_c("chr", seq(1:22))))

  chromosomes <- read.delim("data_for_plotting/genome.fa.fai", header = FALSE,
                            col.names = c("chr", "length", "x", "y", "z"))
  chromosomes <- chromosomes %>%
    filter(chr %in% str_c("chr", seq(1:22))) %>%
    select(chr, length) %>%
    mutate(chr_num = str_remove(chr, "chr") %>% as.numeric())

  mutations <- read.delim("data_for_plotting/somatic_combined.phasing_variants.tsv")

  samples <- combined %>% pull(sample_id) %>% unique()
  lengths <- c(0, 1e6)
  for (sample in samples) {
    for (min_length in lengths) {
      combined %>%
        filter(sample_id == sample) %>%
        filter(length_variants > min_length & !is.na(length_variants)) %>%
        mutate(alternating_01 = factor(row_number() %% 2)) %>%
        mutate(chr_num = str_remove(chr, "chr") %>% as.numeric()) %>%
        ggplot(aes(y = chr_num)) +
        geom_segment(data = chromosomes, aes(y = chr_num, yend = chr_num,
                                             x = 0, xend = length/1e6)) +
        geom_rect(aes(xmin = first_variant_pos/1e6,
                      xmax = last_variant_pos/1e6,
                      ymin = chr_num - .3,
                      ymax = chr_num + .3,
                      fill = alternating_01),
                  show.legend = FALSE) +
        labs(x = "Position (Mb)",
             y = "Chromosome",
             title = str_c("Sample ", sample)) +
        theme_bw(base_size = 20) +
        theme(panel.grid.minor = element_blank()) +
        scale_y_continuous(breaks = seq(1:22),
                           labels = str_c("chr", seq(1:22))) +
        scale_x_continuous(breaks = seq(0, 250, by = 50)) +
        scale_fill_brewer(type = "qual", palette = "Paired") +
        if (min_length == 1e6) {
          ggsave(str_c("phase_sets/exclude_1Mb/", sample, ".pdf"), width = 15, height = 10, useDingbats = FALSE)
        } else if (min_length == 0) {
          ggsave(str_c("phase_sets/all_phase_sets/", sample, ".pdf"), width = 15, height = 10, useDingbats = FALSE)
        }
    }
  }

  # Flipped orientation
  for (sample in samples) {
    for (min_length in lengths) {
      combined %>%
        filter(sample_id == sample) %>%
        filter(length_variants > min_length & !is.na(length_variants)) %>%
        mutate(alternating_01 = factor(row_number() %% 2)) %>%
        mutate(chr_num = str_remove(chr, "chr") %>% as.numeric()) %>%
        ggplot(aes(y = chr_num)) +
        geom_segment(data = chromosomes, aes(y = chr_num, yend = chr_num,
                                             x = 0, xend = length/1e6)) +
        geom_rect(aes(xmin = first_variant_pos/1e6,
                      xmax = last_variant_pos/1e6,
                      ymin = chr_num - .3,
                      ymax = chr_num + .3,
                      fill = alternating_01),
                  show.legend = FALSE) +
        labs(x = "Position (Mb)",
             y = "Chromosome") +
        theme_bw() + #base_size = 20) +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              #axis.text.x = element_blank(),
              #axis.ticks.x = element_blank(),
              #axis.title.x = element_blank(),
              #axis.line.x = element_blank(),
              axis.title.y = element_blank(),
              #axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              plot.margin = margin(l = 20)) +
        scale_y_reverse(position = "right", expand = c(0.025, 0),
                        breaks = seq(1:22), labels = str_c("chr", seq(1:22))) +
        scale_x_reverse(position = "top", expand = c(0, 0), limits = c(250, 0)) +
        scale_fill_brewer(type = "qual", palette = "Paired") +
        if (min_length == 1e6) {
          ggsave(str_c("phase_sets/exclude_1Mb/", sample, ".flipped.pdf"), width = 15, height = 10, useDingbats = FALSE)
          if (sample == "59114_4") {
            ggsave(str_c("main_figures/", sample, ".flipped.pdf"), width = 15, height = 10, useDingbats = FALSE)
          }
        } else if (min_length == 0) {
          ggsave(str_c("phase_sets/all_phase_sets/", sample, ".flipped.pdf"), width = 15, height = 10, useDingbats = FALSE)
        }
    }
  }
}

rm(main, supp)
