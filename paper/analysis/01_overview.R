################################################################################
# Overview figure
# Data available and QC
################################################################################

main = "figures/01_overview/main/"
supp = "figures/01_overview/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# Data available
{
  plot_df <- patient_sample_names_tbl %>%
    mutate(timepoint = factor(timepoint,
                              levels = c("SMM",
                                         "Primary",
                                         "Pre-transplant",
                                         "Post-transplant",
                                         "Remission",
                                         "Relapse",
                                         "Normal"),
                              labels = c("SMM (S)",
                                         "Primary (P)",
                                         "Pre-transplant (PrT)",
                                         "Post-transplant (PoT)",
                                         "Remission (Rem)",
                                         "Relapse (Rel)",
                                         "Normal (N)"),
                              ordered = TRUE))

  ggplot(plot_df, aes(x = timepoint, y = patient)) +
    geom_point(aes(color = my_color_100), shape = 16, size = 3, show.legend = FALSE) +
    geom_point(data = plot_df %>% filter(sorted), shape = 16, color = "#ffffff") +
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
          strip.text = element_text(size = 8)) +
    ggsave(str_c(main, "samples.pdf"),
           width = 2.75,
           height = 2.5,
           useDingbats = FALSE)

  rm(plot_df)
}

# Data QC (all metrics)
{
  data_columns <- sort(c("gems_detected", "mean_dna_per_gem", "bc_on_whitelist",
                         "bc_mean_qscore", "n50_linked_reads_per_molecule",
                         "corrected_loaded_mass_ng", "snps_phased",
                         "genes_phased_lt_100kb", "longest_phase_block",
                         "n50_phase_block", "molecule_length_mean",
                         "molecule_length_stddev", "number_reads",
                         "median_insert_size", "mean_depth", "zero_coverage",
                         "mapped_reads", "pcr_duplication", "r1_q20_bases_fract",
                         "r1_q30_bases_fract", "r2_q20_bases_fract",
                         "r2_q30_bases_fract", "si_q20_bases_fract",
                         "si_q30_bases_fract", "bc_q20_bases_fract",
                         "bc_q30_bases_fract", "large_sv_calls",
                         "short_deletion_calls"))

  data_multiplier <- c(1,1,1,1,1,1e-6,1,1e-3,1e-6,1,1,1e-5,1,1e-3,1e-3,1,1e-6,1e-9,1,1,1,1,1,1e-3,1,1,1,1e3)
  data_multiplier_label <- c(NA,NA,NA,NA,NA,"1e6",NA,"1e3","Mb",NA,NA,"1e5",NA,"Kb","Kb",NA,"Mb","1e9",NA,NA,NA,NA,NA,"1e3",NA,NA,NA,"1e-3")

  multiplier_tbl <- tibble(category = data_columns, data_multiplier, data_multiplier_label)

  plot_df <- lr_summary_tbl %>%
    bind_rows(lr_summary_1000G_tbl) %>%
    gather(data_columns, key = "category", value = "result") %>%
    select(-c(vcf_column, cnv_maf_status, longranger_version)) %>%
    mutate(normal_sample = case_when(timepoint == "Normal" ~ "Normal",
                                     TRUE ~ "Tumor")) %>%
    left_join(multiplier_tbl, by = "category")

  ggplot(plot_df %>% filter(!str_detect(patient, "NA1")),
         aes(x = normal_sample, y = result*data_multiplier)) +
    geom_violin(show.legend = FALSE, draw_quantiles = .5) +
    geom_jitter(aes(color = my_color_100, shape = my_shape),
                height = 0, width = 0.25, show.legend = FALSE) +
    geom_point(data = plot_df %>% filter(str_detect(patient, "NA1")),
               aes(color = my_color_100, shape = my_shape),
               show.legend = FALSE) +
    geom_label(data = plot_df %>%
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
          ) +
    ggsave(str_c(supp, "qc_metrics.pdf"),
           height = 7.5,
           width = 7.5,
           useDingbats = FALSE)

  rm(plot_df, data_columns, data_multiplier, data_multiplier_label, multiplier_tbl)
}

# Data QC (top metrics)
{
  data_columns <- sort(c("gems_detected", "mean_dna_per_gem", "bc_on_whitelist",
                         "bc_mean_qscore", "n50_linked_reads_per_molecule",
                         "corrected_loaded_mass_ng", "snps_phased",
                         "genes_phased_lt_100kb", "longest_phase_block",
                         "n50_phase_block", "molecule_length_mean",
                         "molecule_length_stddev", "number_reads",
                         "median_insert_size", "mean_depth", "zero_coverage",
                         "mapped_reads", "pcr_duplication", "r1_q20_bases_fract",
                         "r1_q30_bases_fract", "r2_q20_bases_fract",
                         "r2_q30_bases_fract", "si_q20_bases_fract",
                         "si_q30_bases_fract", "bc_q20_bases_fract",
                         "bc_q30_bases_fract", "large_sv_calls",
                         "short_deletion_calls"))

  data_multiplier <- c(1,1,1,1,1,1e-6,1,1e-3,1e-6,1,1,1e-5,1,1e-3,1e-3,1,1e-6,1e-9,1,1,1,1,1,1e-3,1,1,1,1e3)
  data_multiplier_label <- c(NA,NA,NA,NA,NA,"1e6",NA,"1e3","Mb",NA,NA,"1e5",NA,"Kb","Kb",NA,"Mb","1e9",NA,NA,NA,NA,NA,"1e3",NA,NA,NA,"1e-3")

  highlight_columns <- sort(c("n50_linked_reads_per_molecule",
                              "n50_phase_block",
                              "molecule_length_mean"))

  multiplier_tbl <- tibble(category = data_columns,
                           data_multiplier, data_multiplier_label) %>%
    filter(category %in% highlight_columns)

  plot_df <- lr_summary_tbl %>%
    bind_rows(lr_summary_1000G_tbl) %>%
    select(patient, sample, timepoint, n50_linked_reads_per_molecule, n50_phase_block, molecule_length_mean, my_color_100, my_shape) %>%
    gather(highlight_columns, key = "category", value = "result") %>%
    mutate(normal_sample = case_when(timepoint == "Normal" ~ "Normal",
                                     TRUE ~ "Tumor")) %>%
    left_join(multiplier_tbl, by = "category") %>%
    mutate(category = factor(category,
                             levels = highlight_columns,
                             labels = c("Molecule Length\n(mean, Kb)",
                                        "Linked-reads per\nmolecule (N50)",
                                        "Phase Set Length\n(N50, Mb)")))

  ggplot(plot_df %>% filter(!str_detect(patient, "NA1")),
         aes(x = normal_sample, y = result*data_multiplier)) +
    geom_violin(show.legend = FALSE, draw_quantiles = .5) +
    geom_jitter(aes(color = my_color_100, shape = my_shape),
                height = 0, width = 0.25, show.legend = FALSE) +
    geom_point(data = plot_df %>% filter(str_detect(patient, "NA1")),
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
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "top_qc_metrics.pdf"),
           useDingbats = FALSE,
           height = 2.5,
           width = 4.25)

  rm(plot_df, data_columns, data_multiplier, data_multiplier_label, highlight_columns, multiplier_tbl)
}

rm(main, supp)
