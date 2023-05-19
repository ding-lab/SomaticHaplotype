################################################################################
# Phase sets figure
# Distribution of phase set lengths
# Relationship of phase set lengths to genomic location
# Heterozygotes/Coverage in unphased regions
################################################################################

data_dir = file.path("data_for_plots/02_overview")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

main = "figures/02_phase_sets/main/"
supp = "figures/02_phase_sets/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)
dir.create(str_c(supp, "/phase_sets"), recursive = TRUE, showWarnings = FALSE)

# Phase set N50 by chromosome

{
  phase_set_by_chromosome_plot_df <- read_tsv(file = file.path(data_dir, "phase_set_by_chromosome_plot_df.tsv"),
                                              show_col_types = FALSE)

  ggplot(data = phase_set_by_chromosome_plot_df,
         aes(x = chromosome, y = N50_broad/1e6)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "phase_set_by_chromosome.pdf"),
         useDingbats = FALSE, width = 7.25, height = 1.5)

  rm(phase_set_by_chromosome_plot_df)
}

# Phase set N50 by sample

{
  phase_set_by_sample_plot_df <- read_tsv(file = file.path(data_dir, "phase_set_by_sample_plot_df.tsv"),
                                          show_col_types = FALSE)

  ggplot(data = phase_set_by_sample_plot_df,
         aes(x = fct_reorder(display_name, sample_n), y = N50_broad/1e6)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "phase_set_by_sample.pdf"),
         useDingbats = FALSE, width = 7.25, height = 1.5)

  rm(phase_set_by_sample_plot_df)
}

# coverage of phase sets
# no relationship between phase set length and coverage if coverage above ~30
if (FALSE) {
  get_coverage <- function(my_coverage_tbl, my_sample, my_chromosome, my_start, my_end){
    my_coverage_tbl %>%
      filter(sample == my_sample,
             chromosome == my_chromosome,
             ((my_start <= start & start <= my_end & my_end <= stop) |
                (start <= my_start & my_end <= stop) |
                (start <= my_start & my_start <= stop & my_end >= stop) |
                (my_start <= start & stop <= my_end))) %>%
      mutate(length_overlap = case_when(my_start <= start & start <= my_end & my_end <= stop ~ as.double(my_end - start),
                                        start <= my_start & my_end <= stop ~ as.double(my_end - my_start),
                                        start <= my_start & my_start <= stop & my_end >= stop ~ as.double(stop - my_start),
                                        TRUE ~ as.double(stop - start))) %>%
      summarize(sum_of_coverage = sum(coverage*length_overlap),
                total_overlap = sum(length_overlap),
                average_coverage = sum_of_coverage/total_overlap) %>%
      pull(average_coverage) %>%
      return()
  }

  coverage_of_phase_sets <- phase_sets_tbl %>%
    filter(chr == "chr1",
           length_variants >= 250000) %>%
    select(ps_id, sample, chr, start, end, sample,
           length_variants, my_color_100) %>%
    rowwise() %>%
    mutate(cov = get_coverage(coverage_tbl, sample, chr, start, end))

  ggplot(coverage_of_phase_sets,
         aes(y = length_variants, x = cov, color = my_color_100)) +
    geom_point(show.legend = FALSE) +
    facet_wrap(~sample) +
    scale_color_identity()

  rm(coverage_of_phase_sets, get_coverage)

}

# Phase set length by deletion
{
  phase_set_deletion_plot_df <- read_tsv(file = file.path(data_dir, "phase_set_deletion_plot_df.tsv"),
                                         show_col_types = FALSE)

  ggplot(data = phase_set_deletion_plot_df,
         aes(x = deleted, y = length_variants/1e6)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "phase_set_deletion.pdf"),
           useDingbats = FALSE, width = 2, height = 1.5)

  rm(phase_set_deletion_plot_df)
}

# Phase sets genome coverage
{
  phase_set_genome_coverage_plot_df <- read_tsv(file = file.path(data_dir, "phase_set_genome_coverage_plot_df.tsv"),
                                                show_col_types = FALSE) %>%
    mutate(length_factor = factor(length_factor, levels = length_label, ordered = TRUE))

  ggplot(phase_set_genome_coverage_plot_df,
         aes(x = length_factor, y = total_length/1e9)) +
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
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")
    )

  ggsave(str_c(supp, "phase_set_genome_coverage.pdf"),
         useDingbats = FALSE, width = 7.25, height = 1.75)

  rm(phase_set_genome_coverage_plot_df)
}

# Phase sets by sample
{
  chr13_chr22_phase_sets_plot_df <- read_tsv(file = file.path(data_dir, "chr13_chr22_phase_sets_plot_df.tsv"),
                                             show_col_types = FALSE)

  sample_names <- c("27522 (P)", "27522 (Rem)")
  n_samples <- length(sample_names)

  ggplot(chr13_chr22_phase_sets_plot_df,
         aes(xmin = start/1e6, xmax = end/1e6,
             ymin = sample_num - 0.4, ymax = sample_num + 0.4)) +
    geom_rect(aes(fill = "#ffffff"), show.legend = FALSE) +
    geom_rect(aes(fill = phase_set_color), show.legend = FALSE) +
    scale_y_continuous(breaks = c(3, 4),
                       expand = c(0.05, 0),
                       limits = c(3 - 0.4, 4 + 0.4),
                       labels = sample_names) +
    scale_x_continuous(expand = c(0,0)) +
    expand_limits(x = 0) +
    scale_fill_identity() +
    facet_wrap(~ chr, ncol = 1, scales = "free_x") +
    labs(x = "Position (Mb)") +
    theme_bw() +
    theme(#panel.grid.major.y = element_blank(),
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
      plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "chr13_chr22_phase_sets.pdf"),
         height = 2.65, width = 5, useDingbats = FALSE)

  rm(chr13_chr22_phase_sets_plot_df, sample_names, n_samples)
}

# phase set length vs. n heterozygous SNPs
{
  phase_set_length_vs_variants_plot_df <- read_tsv(file = file.path(data_dir, "phase_set_length_vs_variants_plot_df.tsv"),
                                                   show_col_types = FALSE)

  ggplot(data = phase_set_length_vs_variants_plot_df,
         aes(x = length_variants/1e6, y = n_variants_total/1e3)) +
    geom_point(aes(color = my_color_100), shape = 16, show.legend = FALSE) +
    geom_smooth(method = "lm") +
    expand_limits(x = 0) +
    scale_color_identity() +
    labs(x = "Phase Set Length (Mb)", y = "Number of Phased Heterozygotes (thousands)") +
    facet_wrap( ~ display_name, ncol = 4) +
    theme_bw() +
    theme(panel.background = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "phase_set_length_vs_variants.pdf"),
           width = 7.5, height = 7.5, useDingbats = FALSE)

  rm(phase_set_length_vs_variants_plot_df)

  phase_set_length_vs_variants.27522_1_plot_df <- read_tsv(file = file.path(data_dir, "phase_set_length_vs_variants.27522_1_plot_df.tsv"),
                                                           show_col_types = FALSE)

  ggplot(data = phase_set_length_vs_variants.27522_1_plot_df,
         aes(x = length_variants/1e6, y = n_variants_total/1e3)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "phase_set_length_vs_variants.27522_1.pdf"),
           width = 7.5, height = 7.5, useDingbats = FALSE)

  rm(phase_set_length_vs_variants.27522_1_plot_df)

}


# draw phase sets for all samples
{

  plot_chr <- read_tsv(file = file.path(data_dir, "plot_chr.tsv"),
                       show_col_types = FALSE)

  phase_sets_dir <- file.path(data_dir, "phase_sets")
  dir.create(path = phase_sets_dir, recursive = TRUE, showWarnings = FALSE)

  for (this_sample in patient_sample_names_tbl$sample) {

    if (patient_sample_names_tbl %>% filter(sample == this_sample) %>% pull(timepoint) != "Normal") {

      print_name <- patient_sample_names_tbl %>% filter(sample == this_sample) %>% pull(display_name)

      plot_df <- read_tsv(file = file.path(phase_sets_dir, str_c(this_sample, ".tsv")),
                          show_col_types = FALSE)

      ggplot(plot_df) +
        geom_segment(data = plot_chr,
                     aes(x = 0, xend = size/1e6, y = chr_num, yend = chr_num)) +
        geom_rect(aes(xmin = start/1e6, xmax = end/1e6,
                      ymin = chr_num - 0.25, ymax = chr_num + 0.25,
                      fill = "#ffffff"), show.legend = FALSE) +
        geom_rect(aes(xmin = start/1e6, xmax = end/1e6,
                      ymin = chr_num - 0.25, ymax = chr_num + 0.25,
                      fill = phase_set_color), show.legend = FALSE) +
        scale_y_continuous(expand = c(0.02, 0),
                           breaks = seq(1,22),
                           labels = str_c("chr", seq(1, 22)),
                           limits = c(1 - 0.25, 22 + 0.25)) +
        scale_x_continuous(expand = c(0,0)) +
        expand_limits(x = 0) +
        scale_fill_identity() +
        labs(x = "Genomic Position (Mb)", y = "Chromosome", title = print_name) +
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
              plot.margin = unit(c(0,0,0,0), "lines"))

      ggsave(str_c(supp, "phase_sets/", this_sample, ".phase_sets.pdf"),
             height = 4.875, width = 7.25, useDingbats = FALSE)

      rm(plot_df, print_name)
    }
  }

  rm(this_sample, plot_chr, phase_sets_dir)
}

rm(main, supp, data_dir)
