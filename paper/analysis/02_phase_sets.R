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
dir.create(str_c(supp, "/phase_sets"), recursive = TRUE, showWarnings = FALSE)

manuscript_numbers[["02_phase_sets"]] <- list()

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

  manuscript_numbers[["02_phase_sets"]][["median_phase_set_length_chr_min"]] <- phase_set_summary_tbl %>% filter(timepoint != "Normal") %>% group_by(chromosome) %>% summarize(median_ps_length = median(N50_broad/1e6)) %>% arrange(median_ps_length) %>% head(1)
  manuscript_numbers[["02_phase_sets"]][["median_phase_set_length_chr_max"]] <- phase_set_summary_tbl %>% filter(timepoint != "Normal") %>% group_by(chromosome) %>% summarize(median_ps_length = median(N50_broad/1e6)) %>% arrange(median_ps_length) %>% tail(1)
  manuscript_numbers[["02_phase_sets"]][["phase_set_length_sd_min"]] <- phase_set_summary_tbl %>% filter(timepoint != "Normal") %>% group_by(chromosome) %>% summarize(median_ps_length = median(N50_broad/1e6), sd_ps_length = sd(N50_broad/1e6)) %>% arrange(sd_ps_length) %>% head(1)
  manuscript_numbers[["02_phase_sets"]][["phase_set_length_sd_max"]] <- phase_set_summary_tbl %>% filter(timepoint != "Normal") %>% group_by(chromosome) %>% summarize(median_ps_length = median(N50_broad/1e6), sd_ps_length = sd(N50_broad/1e6)) %>% arrange(sd_ps_length) %>% tail(1)
  manuscript_numbers[["02_phase_sets"]][["chr21_phase_sets_gt20Mb"]] <- phase_set_summary_tbl %>% filter(timepoint != "Normal") %>% filter(chromosome == "chr21", N50_broad > 20*1e6) %>% nrow()
  manuscript_numbers[["02_phase_sets"]][["chr21_phase_sets_gt20Mb_59114"]] <- phase_set_summary_tbl %>% filter(timepoint != "Normal") %>% filter(chromosome == "chr21", N50_broad > 20*1e6, patient == "59114") %>% nrow()
  manuscript_numbers[["02_phase_sets"]][["25183_qc_measures"]] <- lr_summary_tbl %>% filter(timepoint != "Normal") %>% select(sample, molecule_length_mean, mapped_reads) %>% arrange(molecule_length_mean) %>% tail(1)
}

# Phase set N50 by sample

{
  phase_set_summary_tbl %>%
    filter(timepoint != "Normal") %>%
    ggplot(aes(x = fct_reorder(display_name, sample_n), y = N50_broad/1e6)) +
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

  rm(coverage_of_phase_sets)

}



# Phase set length by deletion
{
  phase_sets_tbl %>%
    filter(sample == "27522_1", length_variants > 1e3) %>%
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

  manuscript_numbers[["02_phase_sets"]][["27522_1_chr13_chr22_phase_set_lengths"]] <- phase_set_summary_tbl %>% filter(sample == "27522_1", chromosome %in% c("chr13", "chr22")) %>% select(sample, N50_broad, chromosome)
  manuscript_numbers[["02_phase_sets"]][["27522_1_overall_phase_set_lengths"]] <- phase_set_summary_tbl %>% group_by(sample) %>% summarize(median_n50_ps_length = median(N50_broad)) %>% filter(sample == "27522_1")
}

# Phase sets genome coverage
{
  plot_df <- phase_sets_tbl %>% filter(!is.na(length_variants),
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
    mutate(length_factor = factor(length_label, levels = length_label, ordered = TRUE))

  ggplot(plot_df, aes(x = length_factor, y = total_length/1e9)) +
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
    ) +
    ggsave(str_c(main, "phase_set_genome_coverage.pdf"),
           useDingbats = FALSE, width = 7.25, height = 1.75)

  manuscript_numbers[["02_phase_sets"]][["total_genome_coverage_Gb"]] <- plot_df %>% pull(total_length) %>% sum()/1e9
  manuscript_numbers[["02_phase_sets"]][["average_genome_coverage_Gb"]] <- plot_df %>% pull(total_length) %>% sum()/(23*1e9)
  manuscript_numbers[["02_phase_sets"]][["length_group_most_ps"]] <- plot_df %>% group_by(length_factor) %>% summarize(count = count, total = sum(plot_df$count), count_pct = 100*count/total) %>% arrange(count_pct) %>% tail(1)
  manuscript_numbers[["02_phase_sets"]][["0-1_Mb_percentage_of_coverage"]] <- plot_df %>% group_by(length_factor) %>% summarize(total_length = total_length/1e9, total_length_all = sum(plot_df$total_length)/1e9, total_length_pct = 100*total_length/total_length_all) %>% filter(length_factor == "0-1")
  manuscript_numbers[["02_phase_sets"]][["1-2_Mb_percentage_of_coverage"]] <- plot_df %>% group_by(length_factor) %>% summarize(count = count, total_length = total_length/1e9, total_length_all = sum(plot_df$total_length)/1e9, total_length_pct = 100*total_length/total_length_all) %>% filter(length_factor == "1-2")
  manuscript_numbers[["02_phase_sets"]][["n_phase_sets_gt_30Mb"]] <- plot_df %>% filter(length_group >= 31) %>% pull(count) %>% sum()
  manuscript_numbers[["02_phase_sets"]][["longest_ps"]] <- plot_df %>% arrange(length_group) %>% tail(1)

  rm(plot_df)
}

# Phase sets by sample
{
  phase_sets_tbl_short <- phase_sets_tbl %>%
    filter(sample %in% c("27522_1", "27522_2"),
           chromosome %in% c("chr13", "chr22")) %>%
    filter(timepoint != "Normal", length_variants > 0, length_variants <= 1e3) %>%
    mutate(phase_set_color = "#636363")

  phase_sets_tbl_long <- phase_sets_tbl %>%
    filter(sample %in% c("27522_1", "27522_2"),
           chromosome %in% c("chr13", "chr22")) %>%
    filter(timepoint != "Normal", length_variants > 1e3) %>%
    mutate(event_number = row_number() %% 2) %>%
    mutate(phase_set_color = case_when(event_number == 0 ~ my_color_100,
                                       TRUE ~ my_color_50))

  plot_df <- phase_sets_tbl_long %>% bind_rows(phase_sets_tbl_short) %>%
    mutate(chr_num = as.numeric(chromosome)) %>%
    mutate(sample_num = as.numeric(display_name))

  sample_names <- c("27522 (P)", "27522 (Rem)")
  n_samples <- length(sample_names)

  ggplot(plot_df,
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
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "chr13_chr22_phase_sets.pdf"),
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
    facet_wrap( ~ display_name, ncol = 4) +
    theme_bw() +
    theme(panel.background = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "phase_set_length_vs_variants.pdf"),
           width = 7.5, height = 7.5, useDingbats = FALSE)

  manuscript_numbers[["02_phase_sets"]][["phase_set_length_vs_n_phased_variants"]] <- summary(lm(n_variants_total ~ length_variants, data = phase_sets_tbl %>% mutate(phase_set_length_Mb = length_variants/1e6, n_variants_total_per1000 = n_variants_total/1000) %>% filter(timepoint != "Normal")))

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


# draw phase sets for all samples
{

  for (this_sample in patient_sample_names_tbl$sample) {

    if (patient_sample_names_tbl %>% filter(sample == this_sample) %>% pull(timepoint) != "Normal") {

      print_name <- patient_sample_names_tbl %>% filter(sample == this_sample) %>% pull(display_name)

      phase_sets_tbl_short <- phase_sets_tbl %>%
        filter(sample == this_sample) %>%
        filter(length_variants > 0, length_variants <= 1e3) %>%
        mutate(phase_set_color = "#636363")

      phase_sets_tbl_long <- phase_sets_tbl %>%
        filter(sample == this_sample) %>%
        filter(length_variants > 1e3) %>%
        mutate(event_number = row_number() %% 2) %>%
        mutate(phase_set_color = case_when(event_number == 0 ~ my_color_100,
                                           TRUE ~ my_color_50))

      plot_df <- phase_sets_tbl_long %>% bind_rows(phase_sets_tbl_short) %>%
        mutate(chr_num = as.numeric(chromosome))

      plot_chr <- chromosome_length_tbl %>% mutate(chr_num = as.numeric(contig))

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
              plot.margin = unit(c(0,0,0,0), "lines")) +
        ggsave(str_c(supp, "phase_sets/", this_sample, ".phase_sets.pdf"),
               height = 4.875, width = 7.25, useDingbats = FALSE)

      rm(phase_sets_tbl_short, phase_sets_tbl_long, plot_df, print_name, plot_chr)
    }
  }

  rm(this_sample)
}


manuscript_numbers[["02_phase_sets"]][["HLA_stats"]] <- phase_sets_tbl %>%
  filter(chr == "chr6", timepoint != "Normal") %>%
  filter((start <= 28510120 & 28510120 <= end & end <= 33480577) |
           (28510120 <= start & end <= 33480577) |
           (28510120 <= start & start <= 33480577 & end >= 33480577) |
           (start <= 28510120 & 33480577 <= end)) %>% filter(length_variants >= 1e3) %>% rowwise() %>% mutate(coverage = min(end, 33480577) - max(start, 28510120)) %>% ungroup() %>% group_by(sample) %>% summarize(total_coverage = sum(coverage)/(33480577 - 28510120), total = n()) %>% ungroup() %>% summarize(total_coverage_median = median(total_coverage), total_coverage_min = min(total_coverage), total_coverage_max = max(total_coverage), n_phase_sets_median = median(total), n_phase_sets_min = min(total), n_phase_sets_max = max(total))

rm(main, supp)
