################################################################################
# Closer look at allele relationships
################################################################################

main = "figures/04_alleles/main/"
supp = "figures/04_alleles/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

manuscript_numbers[["04_alleles"]] <- list()

# mutation coverage at CNV neutral sites
{
  cnv_neutral_mutation_sites <- variant_pairs_mapq20_tbl %>%
    filter(cnv_position1 > -0.25,
           cnv_position1 < 0.2) %>%
    mutate(Variant = Variant1) %>%
    select(sample, Variant) %>%
    unique()
  cnv_neutral_mutation_sites <- bind_rows(cnv_neutral_mutation_sites,
                                          variant_pairs_mapq20_tbl %>%
                                            filter(cnv_position2 > -0.25,
                                                   cnv_position2 < 0.2) %>%
                                            mutate(Variant = Variant2) %>%
                                            select(sample, Variant) %>%
                                            unique()) %>%
    unique()

  # coverage at mutation sites

  cnv_neutral_mutation_sites %>%
    left_join(phasing_variants_mapq20_tbl,
              by = c("sample", "Variant")) %>%
    group_by(sample, my_color_100, Variant) %>%
    summarize(phased_barcode_coverage = sum(barcode_REF_H1, barcode_REF_H2,
                                            barcode_ALT_H1, barcode_ALT_H2),
              phased_alt_barcode_coverage = sum(barcode_ALT_H1, barcode_ALT_H2)) %>%
    filter(phased_barcode_coverage > 0) %>%
    ungroup() %>%
    ggplot(aes(x = log2(phased_barcode_coverage),
               y = phased_alt_barcode_coverage)) +
    geom_vline(xintercept = log2(10), lty = 2) +
    geom_vline(xintercept = log2(100), lty = 2) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_point(aes(color = my_color_100), shape = 16, alpha = 0.25) +
    scale_color_identity() +
    facet_wrap(~sample, ncol = 1) +
    labs(x = "Phased Barcodes Covering Mutation Site (log2)",
         y = "Phased Barcodes Covering Mutation Site and Supporting Mutant Allele (non-transformed)") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.grid = element_line(size = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 8)) +
    ggsave(str_c(supp, "phased_coverage.pdf"),
           width = 7.25,
           height = 6.75,
           useDingbats = FALSE)

  # distance between mutations

  cnv_neutral_good_coverage_sites <- cnv_neutral_mutation_sites %>%
    left_join(phasing_variants_mapq20_tbl,
              by = c("sample", "Variant")) %>%
    group_by(sample, my_color_100, Variant) %>%
    summarize(phased_barcode_coverage = sum(barcode_REF_H1, barcode_REF_H2,
                                            barcode_ALT_H1, barcode_ALT_H2),
              phased_alt_barcode_coverage = sum(barcode_ALT_H1, barcode_ALT_H2)) %>%
    ungroup() %>%
    filter(phased_barcode_coverage >= 10) %>%
    filter(phased_barcode_coverage <= 100) %>%
    filter(phased_alt_barcode_coverage > 0)

  variant_pairs_mapq20_tbl_cnv_neutral_good_coverage <- NULL

  for (sample_id in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    this_sample <- cnv_neutral_good_coverage_sites %>% filter(sample == sample_id)
    if (is.null(variant_pairs_mapq20_tbl_cnv_neutral_good_coverage)) {
      variant_pairs_mapq20_tbl_cnv_neutral_good_coverage <- variant_pairs_mapq20_tbl %>%
        filter(sample == sample_id,
               Variant1 %in% this_sample$Variant,
               Variant2 %in% this_sample$Variant)
    } else {
      variant_pairs_mapq20_tbl_cnv_neutral_good_coverage <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
        bind_rows(variant_pairs_mapq20_tbl %>%
                    filter(sample == sample_id,
                           Variant1 %in% this_sample$Variant,
                           Variant2 %in% this_sample$Variant))
    }
  }

  median_molecule_length <- lr_summary_tbl %>%
    filter(timepoint != "Normal") %>%
    pull(molecule_length_mean) %>%
    median() %>%
    plyr::round_any(2000)

  max_overlapping_bx <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants <= median_molecule_length) %>%
    filter(distance_between_variants >= 100) %>%
    pull(n_overlapping_barcodes) %>% max() %>% plyr::round_any(10, f = ceiling)

  p <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants <= median_molecule_length) %>%
    filter(distance_between_variants >= 100) %>%
    mutate(overlap_category = case_when(n_overlapping_barcodes == 0 ~ "Zero",
                                        n_overlapping_barcodes <= 10 ~ "<= 10",
                                        n_overlapping_barcodes <= 20 ~ "<= 20",
                                        n_overlapping_barcodes <= max_overlapping_bx ~ str_c("<= ", as.character(max_overlapping_bx)))) %>%
    mutate(overlap_category = fct_rev(factor(overlap_category,
                                             levels = c("Zero", "<= 10", "<= 20", str_c("<= ", as.character(max_overlapping_bx))),
                                             ordered = TRUE))) %>%
    ggplot(aes(x = distance_between_variants/1000)) +
    geom_histogram(aes(fill = overlap_category),
                   boundary = 0, binwidth = 1000/500) +
    scale_fill_viridis(discrete = TRUE, option = "C") +
    labs(x = "Distance Between Somatic Mutations (Kb)",
         y = "Pairs of Somatic Mutations",
         fill = "Overlapping\nBarcodes") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"))

    ggsave(str_c(main, "overlapping_barcodes.with_legend.pdf"), p,
         width = 2, height = 2, useDingbats = FALSE)
    ggsave(str_c(main, "overlapping_barcodes.without_legend.pdf"), p + guides(fill = FALSE),
           width = 2, height = 2, useDingbats = FALSE)

  variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants < 100) %>%
    mutate(overlap_category = case_when(n_overlapping_barcodes == 0 ~ "Zero",
                                        n_overlapping_barcodes <= 10 ~ "<= 10",
                                        n_overlapping_barcodes <= 20 ~ "<= 20",
                                        n_overlapping_barcodes <= 50 ~ "<= 50",
                                        n_overlapping_barcodes <= 100 ~ "<= 100",
                                        n_overlapping_barcodes > 100 ~ "> 100")) %>%
    mutate(overlap_category = fct_rev(factor(overlap_category,
                                             levels = c("Zero", "<= 10", "<= 20",
                                                        "<= 50", "<= 100", "> 100"),
                                             ordered = TRUE))) %>%
    ggplot(aes(x = distance_between_variants)) +
    geom_histogram(aes(fill = overlap_category),
                   boundary = 0.5, binwidth = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Distance Between Somatic Mutations",
         y = "Pairs of Somatic Mutations",
         fill = "Overlapping\nBarcodes") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "overlapping_barcodes.lt_100.pdf"),
           width = 7.25, height = 2, useDingbats = FALSE)

  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_gt62kb"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants > median_molecule_length) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_gt62kb_no_shared_barcodes"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants > median_molecule_length) %>% filter(n_overlapping_barcodes == 0) %>% nrow()
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_gt62kb_no_shared_barcodes"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_gt62kb_no_shared_barcodes"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_gt62kb"]]

  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length) %>% filter(n_overlapping_barcodes >= 1) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap_gte1_mutation"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length) %>% filter(n_overlapping_barcodes >= 1) %>% filter(n_barcodes_with_mutation >= 1) %>% nrow()

  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants >= 100) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_0_overlap"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants >= 100) %>% filter(n_overlapping_barcodes == 0) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_lte10_overlap"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants >= 100) %>% filter(n_overlapping_barcodes <= 10, n_overlapping_barcodes > 0) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_lte20_overlap"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants >= 100) %>% filter(n_overlapping_barcodes <= 20, n_overlapping_barcodes > 10) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_gt20_overlap"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants >= 100) %>% filter(n_overlapping_barcodes > 20) %>% nrow()

  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_0_overlap"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_0_overlap"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp"]]
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_lte10_overlap"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_lte10_overlap"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp"]]
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_lte20_overlap"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_lte20_overlap"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp"]]
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_gt20_overlap"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp_gt20_overlap"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte100bp"]]

  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_lt100bp"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants < 100) %>% nrow()
  manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_lt100bp_0_overlap"]] <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% filter(distance_between_variants <= median_molecule_length, distance_between_variants < 100) %>% filter(n_overlapping_barcodes == 0) %>% nrow()
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_lt100bp_0_overlap"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_lt100bp_0_overlap"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_lt100bp"]]

  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb"]] <- 100*manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage"]]
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap"]] <- 100*(manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb"]])
  manuscript_numbers[["04_alleles"]][["pct_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap_gte1_mutation"]] <- 100*(manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap_gte1_mutation"]]/manuscript_numbers[["04_alleles"]][["n_variant_pairs_cnv_neutral_good_coverage_lte62kb_gte1_overlap"]])

  rm(cnv_neutral_mutation_sites, cnv_neutral_good_coverage_sites, max_overlapping_bx, sample_id, this_sample)
}

# automated allele combination stats
{
  # allele combinations (handmade with known figures)

  n_total <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>% nrow()

  n_within_length <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants <= median_molecule_length) %>% nrow()
  p_within_length <- 100*n_within_length/n_total

  n_within_length_share_barcode <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants <= median_molecule_length) %>%
    filter(n_overlapping_barcodes > 0) %>% nrow()
  p_within_length_share_barcode <- 100*n_within_length_share_barcode/n_within_length

  n_within_length_share_barcode_share_variant <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants <= median_molecule_length) %>%
    filter(n_overlapping_barcodes > 0) %>%
    filter(n_barcodes_with_mutation > 0) %>% nrow()
  p_within_length_share_barcode_share_variant <- 100*n_within_length_share_barcode_share_variant/n_within_length_share_barcode

  plot_df <- tribble(~level1, ~level2, ~value,
                     "Within\nLength",  "Yes",                                   p_within_length,
                     "Within\nLength",   "No",                             100 - p_within_length,
                     "Share\nBarcode",  "Yes",                     p_within_length_share_barcode,
                     "Share\nBarcode",   "No",               100 - p_within_length_share_barcode,
                     "Share\nVariant",  "Yes",       p_within_length_share_barcode_share_variant,
                     "Share\nVariant",   "No", 100 - p_within_length_share_barcode_share_variant) %>%
    mutate(level1 = factor(level1,
                           levels = c("Within\nLength", "Share\nBarcode", "Share\nVariant"),
                           labels = c(str_c("Within\n", median_molecule_length/1000, " Kb"), "Share\nBarcode", "Share\nVariant"),
                           ordered = TRUE))

  bar_width = 0.5
  ggplot(data = plot_df, aes(x = level1, y = value, fill = level2)) +
    geom_bar(stat = "identity", width = bar_width, show.legend = FALSE) +
    geom_text(data = plot_df %>% filter(level2 == "Yes"),
              aes(label = str_c(round(value, 1), "%")),
              vjust = 1, nudge_y = -0.5,
              color = "white",
              size = 6/ggplot2:::.pt,
              fontface = "bold") +
    geom_text(data = plot_df %>% filter(level2 == "No"),
              aes(label = str_c(round(value, 1), "%"), y = 100),
              vjust = 1, nudge_y = -0.5,
              color = "white",
              size = 6/ggplot2:::.pt,
              fontface = "bold") +
    geom_segment(x = 1 + bar_width/2, y = p_within_length, xend = 2 - bar_width/2, yend = 100,
                 lty = 2, color = "#bdbdbd", lwd = 0.5) +
    geom_segment(x = 1 + bar_width/2, y = 0, xend = 2 - bar_width/2, yend = 0,
                 lty = 2, color = "#bdbdbd", lwd = 0.5) +
    geom_segment(x = 2 + bar_width/2, y = p_within_length_share_barcode, xend = 3 - bar_width/2, yend = 100,
                 lty = 2, color = "#bdbdbd", lwd = 0.5) +
    geom_segment(x = 2 + bar_width/2, y = 0, xend = 3 - bar_width/2, yend = 0,
                 lty = 2, color = "#bdbdbd", lwd = 0.5) +
    geom_segment(x = 3 + bar_width/2, y = p_within_length_share_barcode_share_variant, xend = 4 - bar_width/2, yend = 100,
                 lty = 2, color = "#bdbdbd", lwd = 0.5) +
    geom_segment(x = 3 + bar_width/2, y = 0, xend = 4 - bar_width/2, yend = 0, lty = 2,
                 color = "#bdbdbd", lwd = 0.5) +
    scale_fill_manual(values = c("#80cdc1", "#01665e")) +
    scale_x_discrete(expand = c(0,0)) +
    expand_limits(x = 4 - bar_width/2) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    labs(x = NULL, y = "Somatic Mutation Pairs (%)", fill = NULL) +
    ggsave(str_c(main, "n_variant_pairs.pdf"),
           width = 2, height = 2, useDingbats = FALSE)

  ### Other parts

  text_df <- variant_pairs_mapq20_tbl_cnv_neutral_good_coverage %>%
    filter(distance_between_variants <= median_molecule_length) %>%
    filter(n_overlapping_barcodes > 0) %>%
    filter(n_barcodes_with_mutation > 0) %>%
    mutate(category = str_c(as.numeric(n_bx_overlap_00 > 0),
                            as.numeric(n_bx_overlap_01 > 0),
                            as.numeric(n_bx_overlap_10 > 0),
                            as.numeric(n_bx_overlap_11 > 0))) %>%
    group_by(category) %>%
    summarize(count = n()) %>%
    mutate(category = case_when(str_starts(category, "0") ~ "Other",
                                category == "1111" ~ "Other",
                                TRUE ~ category)) %>%
    group_by(category) %>%
    summarize(total = sum(count)) %>%
    mutate(category = factor(category,
                             levels = c("1001", "1010", "1100", "1110",
                                        "1011", "1101", "Other"),
                             ordered = TRUE))

  plot_df <- tribble(~alleles, ~category, ~filled,
                     4, "1001", 1,
                     3, "1001", 0,
                     2, "1001", 0,
                     1, "1001", 1,
                     4, "1010", 1,
                     3, "1010", 0,
                     2, "1010", 1,
                     1, "1010", 0,
                     4, "1100", 1,
                     3, "1100", 1,
                     2, "1100", 0,
                     1, "1100", 0,
                     4, "1110", 1,
                     3, "1110", 1,
                     2, "1110", 1,
                     1, "1110", 0,
                     4, "1011", 1,
                     3, "1011", 0,
                     2, "1011", 1,
                     1, "1011", 1,
                     4, "1101", 1,
                     3, "1101", 1,
                     2, "1101", 0,
                     1, "1101", 1,
                     4, "Other", 2,
                     3, "Other", 2,
                     2, "Other", 2,
                     1, "Other", 2) %>%
    mutate(category = factor(category,
                             levels = c("1001", "1010", "1100", "1110",
                                        "1011", "1101", "Other"),
                             ordered = TRUE))

  ggplot() +
    geom_tile(data = plot_df,
              aes(x = category, y = alleles, fill = as.factor(filled)),
              #color = NA,
              width = 0.8,
              height = 0.95) +
    geom_text(data = text_df, aes(x = category, y = 5.5, label = total),
              size = 2) +
    geom_text(data = text_df, aes(x = category, y = 5, label = str_c(round(100*total/n_within_length_share_barcode_share_variant, 1), "%")),
              size = 2) +
    scale_fill_manual(values = c("#80cdc1", "#01665e", "#bdbdbd")) +
    scale_y_continuous(breaks = c(seq(1,5), 5.5),
                       labels = c("ALT/ALT", "ALT/REF", "REF/ALT", "REF/REF", "", "Total"),
                       position = "left") +
    scale_x_discrete(expand = c(0,0), position = "bottom") +
    labs(x = "Combinations of Linked Somatic Mutations", y = NULL) +
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "allele_combinations.pdf"), height = 2, width = 3, useDingbats = FALSE)

  manuscript_numbers[["04_alleles"]][["allele_pair_combinations"]] <- text_df %>% mutate(pct = 100*total/n_within_length_share_barcode_share_variant)

  rm(n_total, n_within_length, p_within_length,
     n_within_length_share_barcode, p_within_length_share_barcode,
     n_within_length_share_barcode_share_variant, p_within_length_share_barcode_share_variant,
     plot_df, text_df,
     variant_pairs_mapq20_tbl_cnv_neutral_good_coverage,
     median_molecule_length,
     bar_width)
}

# simplified variant pairs
{

  plot_barcode_variants <- function(barcodes_variants,
                                    sample,
                                    variants_of_interest,
                                    var1,
                                    var2,
                                    gene1,
                                    gene2,
                                    protein1,
                                    protein2,
                                    position1,
                                    position2,
                                    output_dir,
                                    mutation_pattern,
                                    output_filename){

    bcv <- barcodes_variants[[sample]] %>%
      filter(Somatic_Variant %in% variants_of_interest)

    bcv_filtered <- bcv %>%
      filter(Allele != "No Coverage" &
               ((Filter %in% c("PASS", "10X_PHASING_INCONSISTENT") &
                   Genotype != "1|1") |
                  Variant %in% variants_of_interest))

    if (bcv_filtered %>% pull(Somatic_Variant) %>% unique() %>% length() == 1) {
      return("Only one somatic variant with passing coverage.")
    }

    variants_more_than_1 <- bcv_filtered %>%
      select(Barcode, Variant) %>%
      unique() %>%
      group_by(Variant) %>%
      summarize(count = n()) %>%
      filter(count > 1) %>%
      mutate(Variant_rank = seq(1:n()))

    barcodes_more_than_1 <- bcv_filtered %>%
      select(Barcode, Variant) %>%
      unique() %>%
      filter(Variant %in% variants_more_than_1$Variant) %>%
      group_by(Barcode) %>%
      summarize(count = n()) %>%
      filter(count > 1) %>%
      pull(Barcode)

    ordered_barcodes <- bcv_filtered %>%
      filter(Barcode %in% barcodes_more_than_1) %>%
      group_by(Barcode) %>%
      summarize(proportion_H1 = mean(Haplotype == "H1"),
                proportion_H2 = mean(Haplotype == "H2")) %>%
      arrange(desc(proportion_H2), proportion_H1) %>%
      select(Barcode) %>%
      mutate(Barcode_rank = seq(1:n())) %>%
      mutate(Barcode_ordered = factor(Barcode_rank, labels = Barcode),
             ordered = TRUE)

    plot_data <- bcv_filtered %>%
      filter(Variant %in% variants_more_than_1$Variant,
             Barcode %in% barcodes_more_than_1) %>%
      left_join(ordered_barcodes, by = "Barcode") %>%
      left_join(variants_more_than_1, by = "Variant")

    horizontal_lines <- plot_data %>%
      group_by(Barcode) %>%
      summarize(y = unique(Barcode_rank),
                yend = unique(Barcode_rank),
                x = min(Variant_rank),
                xend = max(Variant_rank)) %>%
      ungroup()

    n_variants <- plot_data %>% pull(Variant) %>% unique() %>% length()
    n_barcodes <- plot_data %>% pull(Barcode) %>% unique() %>% length()

    plot_data %>%
      left_join(horizontal_lines, by = "Barcode") %>%
      mutate(Haplotype = factor(Haplotype,
                                levels = c("H1", "H2", "Not_Phased/Not_Heterozygote/Allele_Not_Matching"),
                                labels = c("Haplotype 1", "Haplotype 2", "Not Phased"),
                                ordered = TRUE)) %>%
      ggplot(aes(x = Variant_rank, y = Barcode_rank,
                 label = Allele, fill = Haplotype)) +
      geom_vline(xintercept = plot_data %>% filter(Variant == var1) %>%
                   pull(Variant_rank) %>% unique(),
                 color = "red") +
      geom_vline(xintercept = plot_data %>% filter(Variant == var2) %>%
                   pull(Variant_rank) %>% unique(),
                 color = "red") +
      geom_segment(aes(y = y, yend = yend, x = x, xend = xend)) +
      geom_tile() +
      geom_text(size = 2, color = "#ffffff") +
      labs(x = "Variant", y = "Barcode", fill = NULL,
           title = sample,
           subtitle = str_c(gene1, " ", protein1, " ", "(", var1, ")", "\n",
                            gene2, " ", protein2, " ", "(", var2, ")")) +
      scale_fill_brewer(palette = "Dark2", drop = FALSE) +
      scale_x_continuous(breaks = seq(1:n_variants),
                         labels = plot_data %>% pull(Variant) %>% sort() %>% unique(),
                         expand = c(0,0)) +
      scale_y_continuous(breaks = seq(1:n_barcodes),
                         labels = plot_data %>% arrange(Barcode_rank) %>% pull(Barcode) %>% unique(),
                         expand = c(0,0)) +
      coord_equal() +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            axis.text.y = element_text(size = 6),
            axis.title = element_text(size = 8),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_line(size = 0.5),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 8),
            plot.margin = unit(c(0,0,0,0), "lines")
      ) +
      ggsave(str_c(output_dir, output_filename), width = 7.5, height = 10, useDingbats = FALSE)
  }

  pv <- phasing_variants_mapq20_tbl
  vp <- variant_pairs_mapq20_tbl %>%
    filter(#cnv_over_range > -.25,
           #cnv_over_range < .2) ,
           cnv_over_range == cnv_position1,
           cnv_over_range == cnv_position2)

  perfect_alt_on_H1 <- pv %>%
    filter(barcode_ALT_H1 >= 2 | barcode_ALT_H2 >= 2) %>%
    filter(pct_ALT_on_H1 >= 0.91) %>%
    select(Variant, sample) %>% unique() %>% mutate(combo = str_c(Variant, sample))

  perfect_alt_on_H2 <- pv %>%
    filter(barcode_ALT_H1 >= 2 | barcode_ALT_H2 >= 2) %>%
    filter(pct_ALT_on_H1 <= 0.09) %>%
    select(Variant, sample) %>% unique() %>% mutate(combo = str_c(Variant, sample))

  phased_pairs <- vp %>%
    mutate(mutation_pattern = case_when( (n_bx_overlap_01 > 0 & n_bx_overlap_10 > 0 & n_bx_overlap_11 == 0) ~ 'independent',
                                         (n_bx_overlap_01 > 0 & n_bx_overlap_10 == 0 & n_bx_overlap_11 > 0) ~ 'var2>var1',
                                         (n_bx_overlap_01 == 0 & n_bx_overlap_10 > 0 & n_bx_overlap_11 > 0) ~ 'var1>var2',
                                         (n_bx_overlap_01 == 0 & n_bx_overlap_10 == 0 & n_bx_overlap_11 > 0) ~ 'co-occur')) %>%
    filter(!is.na(mutation_pattern)) %>%
    filter(distance_between_variants >= 10,
           n_overlapping_barcodes >= 10,
           n_barcodes_with_mutation >= 5) %>%
    filter(distance_between_variants >= 100 | n_barcodes_with_mutation >= 10) %>%
    mutate(combo1 = str_c(Variant1, sample),
           combo2 = str_c(Variant2, sample)) %>%
    filter( (combo1 %in% perfect_alt_on_H1$combo & combo2 %in% perfect_alt_on_H1$combo) |
              (combo1 %in% perfect_alt_on_H2$combo & combo2 %in% perfect_alt_on_H2$combo) |
              (combo1 %in% perfect_alt_on_H1$combo & combo2 %in% perfect_alt_on_H2$combo) |
              (combo1 %in% perfect_alt_on_H2$combo & combo2 %in% perfect_alt_on_H1$combo))

  for (i in 1:nrow(phased_pairs)) {

    this_mutation_pattern <- phased_pairs$mutation_pattern[i]

    dir.create(str_c(supp, "phased_pairs/", this_mutation_pattern), recursive = TRUE, showWarnings = FALSE)

    this_var1 <- phased_pairs$Variant1[i]
    this_var2 <- phased_pairs$Variant2[i]

    this_variants_of_interest <- c(this_var1, this_var2)

    this_sample <- phased_pairs$sample[i]

    this_gene1 <- "Gene1"
    this_protein1 <- "Protein1"

    this_position1 <- phased_pairs$position1[i]

    this_gene2 <- "Gene2"
    this_protein2 <- "Protein2"

    this_position2 <- phased_pairs$position2[i]

    this_output_filename <- str_c(this_sample, i, "barcode_variants.pdf", sep = ".")
    print(this_output_filename)

    plot_barcode_variants(barcodes_variants = barcodes_variants_mapq20_tbl,
                          sample = this_sample,
                          variants_of_interest = this_variants_of_interest,
                          var1 = this_var1,
                          var2 = this_var2,
                          gene1 = this_gene1,
                          gene2 = this_gene2,
                          protein1 = this_protein1,
                          protein2 = this_protein2,
                          position1 = this_position1,
                          position2 = this_position2,
                          output_dir =  str_c(supp, "phased_pairs/", this_mutation_pattern, "/"),
                          mutation_pattern = this_mutation_pattern,
                          output_filename = this_output_filename)

  }

  {
    pv <- phasing_variants_driver_mapq20_tbl
    vp <- variant_pairs_driver_mapq20_tbl

    perfect_alt_on_H1 <- pv %>%
      filter(barcode_ALT_H1 >= 2 | barcode_ALT_H2 >= 2) %>%
      filter(pct_ALT_on_H1 >= .91) %>%
      select(Variant, sample) %>%
      unique() %>%
      mutate(combo = str_c(Variant, sample))

    perfect_alt_on_H2 <- pv %>%
      filter(barcode_ALT_H1 >= 2 | barcode_ALT_H2 >= 2) %>%
      filter(pct_ALT_on_H1 <= 0.09) %>%
      select(Variant, sample) %>%
      unique() %>%
      mutate(combo = str_c(Variant, sample))

    phased_pairs <- vp %>%
      mutate(mutation_pattern = case_when( (n_bx_overlap_01 > 0 & n_bx_overlap_10 > 0 & n_bx_overlap_11 == 0) ~ 'independent',
                                           (n_bx_overlap_01 > 0 & n_bx_overlap_10 == 0 & n_bx_overlap_11 > 0) ~ 'var2>var1',
                                           (n_bx_overlap_01 == 0 & n_bx_overlap_10 > 0 & n_bx_overlap_11 > 0) ~ 'var1>var2')) %>%
      filter(!is.na(mutation_pattern)) %>%
      mutate(combo1 = str_c(Variant1, sample),
             combo2 = str_c(Variant2, sample)) %>%
      filter( (combo1 %in% perfect_alt_on_H1$combo & combo2 %in% perfect_alt_on_H1$combo) |
                (combo1 %in% perfect_alt_on_H2$combo & combo2 %in% perfect_alt_on_H2$combo) |
                (combo1 %in% perfect_alt_on_H1$combo & combo2 %in% perfect_alt_on_H2$combo) |
                (combo1 %in% perfect_alt_on_H2$combo & combo2 %in% perfect_alt_on_H1$combo))

    for (i in 1:nrow(phased_pairs)) {

      this_mutation_pattern <- phased_pairs$mutation_pattern[i]

      dir.create(str_c(supp, "phased_pairs_drivers/", this_mutation_pattern), recursive = TRUE, showWarnings = FALSE)

      this_var1 <- phased_pairs$Variant1[i]
      this_var2 <- phased_pairs$Variant2[i]

      this_variants_of_interest <- c(this_var1, this_var2)

      this_sample <- phased_pairs$sample[i]

      this_gene1 <- driver_mutations_tbl %>%
        mutate(var = str_c(chr, pos, ref, alt, sep = ":")) %>%
        filter(var == this_var1) %>% pull(gene)
      this_protein1 <- driver_mutations_tbl %>%
        mutate(var = str_c(chr, pos, ref, alt, sep = ":")) %>%
        filter(var == this_var1) %>% pull(protein)
      this_position1 <- phased_pairs$position1[i]

      this_gene2 <- driver_mutations_tbl %>%
        mutate(var = str_c(chr, pos, ref, alt, sep = ":")) %>%
        filter(var == this_var2) %>% pull(gene)
      this_protein2 <- driver_mutations_tbl %>%
        mutate(var = str_c(chr, pos, ref, alt, sep = ":")) %>%
        filter(var == this_var2) %>% pull(protein)
      this_position2 <- phased_pairs$position2[i]

      this_output_filename <- str_c(this_sample, i, "barcode_variants.pdf", sep = ".")
      print(this_output_filename)

      plot_barcode_variants(barcodes_variants = barcodes_variants_driver_mapq20_tbl,
                            sample = this_sample,
                            variants_of_interest = this_variants_of_interest,
                            var1 = this_var1,
                            var2 = this_var2,
                            gene1 = this_gene1,
                            gene2 = this_gene2,
                            protein1 = this_protein1,
                            protein2 = this_protein2,
                            position1 = this_position1,
                            position2 = this_position2,
                            output_dir =  str_c(supp, "phased_pairs_drivers/", this_mutation_pattern, "/"),
                            mutation_pattern = this_mutation_pattern,
                            output_filename = this_output_filename)

    }
  }

}

manuscript_numbers[["04_alleles"]][["27522_NRAS_allele_combinations"]] <- variant_pairs_driver_mapq20_tbl %>% filter(sample == "27522_1", Variant1 == "chr1:114713909:G:T", Variant2 == "chr1:114716124:C:G") %>% select(starts_with("n_bx_overlap"))
manuscript_numbers[["04_alleles"]][["27522_NRAS_VAFs"]] <- driver_mutations_vaf_tbl %>% filter(patient == "27522", gene == "NRAS") %>% select(sample, gene, vaf, protein)

# plot interesting allele pairs
{
  plot_allele_pairs <- function(sample_id, gene, variant1, variant2, protein1, protein2, barcodes_variants){
    barcodes_variants[[sample_id]] %>%
      filter(Somatic_Variant %in% c(variant1, variant2),
             Variant %in% c(variant1, variant2)) %>%
      group_by(Barcode) %>%
      summarize(V1_0 = any(Variant == variant1 & Allele == 0),
                V1_1 = any(Variant == variant1 & Allele == 1),
                V2_0 = any(Variant == variant2 & Allele == 0),
                V2_1 = any(Variant == variant2 & Allele == 1)) %>%
      mutate(V1_allele = case_when(V1_0 == 1 & V1_1 == 0 ~ "REF",
                                   V1_0 == 0 & V1_1 == 1 ~ "ALT",
                                   TRUE ~ "NC"),
             V2_allele = case_when(V2_0 == 1 & V2_1 == 0 ~ "REF",
                                   V2_0 == 0 & V2_1 == 1 ~ "ALT",
                                   TRUE ~ "NC")) %>%
      group_by(V1_allele, V2_allele) %>%
      summarize(count = n()) %>%
      ungroup() %>%
      arrange(desc(count)) %>%
      mutate(allele_pair = str_c(V1_allele, V2_allele, sep = "/")) %>%
      mutate(allele_pair = str_c(allele_pair, " (", count, ")")) %>%
      mutate(V1_allele = factor(V1_allele, levels = c("REF", "ALT", "NC"), ordered = TRUE),
             V2_allele = factor(V2_allele, levels = c("REF", "ALT", "NC"), ordered = TRUE)) %>%
      gather(V1_allele, V2_allele, key = "variant", value = "allele") %>%
      mutate(allele = factor(allele, levels = c("REF", "ALT", "NC"), ordered = TRUE)) %>%
      ggplot(aes(x = variant, y = fct_reorder(allele_pair, count))) +
      geom_point(aes(shape = allele), color = "#000000", fill = "#FFFFFF",
                 size = 3, show.legend = FALSE) +
      scale_shape_manual(values = c(21, 25, 4)) +
      scale_x_discrete(labels = c(protein1, protein2)) +
      scale_y_discrete(expand = c(0.04, 0.04)) +
      labs(x = gene, y = NULL) +
      theme_bw() +
      theme(panel.grid.major.y = element_line(size = 2),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title.x = element_text(size = 8),
            strip.background = element_blank(),
            strip.text = element_text(size = 8),
            plot.margin = unit(c(0,0,0,0), "lines"),
            legend.background = element_blank(),
            legend.text = element_text(size = 6)) +
      ggsave(str_c(main, "allele_pairs.", sample_id, ".", gene, ".pdf"),
             width = 1.5, height = 1.75, useDingbats = FALSE)
  }

  plot_allele_pairs(sample_id = "27522_1",
                    variant1 = "chr1:114713909:G:T",
                    variant2 = "chr1:114716124:C:G",
                    gene = "NRAS",
                    protein1 = "Q61K",
                    protein2 = "G13R",
                    barcodes_variants = barcodes_variants_driver_mapq20_tbl)

  plot_allele_pairs(sample_id = "57075_3",
                    variant1 = "chr21:34792263:A:C",
                    variant2 = "chr21:34792313:T:G",
                    gene = "RUNX1",
                    protein1 = "S439A",
                    protein2 = "E422A",
                    barcodes_variants = barcodes_variants_driver_mapq20_tbl)



}

# 27522_1 NRAS G13 Q61 figure
if (FALSE) {
  bcv <- read.delim("~/Desktop/mmy_paper/27522_1.NRAS/27522_1.NRAS.barcodes_variants.tsv")

  # good_barcodes <- tribble(~Barcode, ~barcode_rank, ~best_haplotype,
  #                          "TGTAGTGAGCTGCCCA-1", 1, "Not Phased",
  #                          "GAATGAATCTGAGAGG-1", 2, "Haplotype 1",
  #                          "CTGTGCTAGGAAGCCT-1", 3, "Not Phased",
  #                          "CCCGGAACATGGAGAC-1", 4, "Haplotype 2",
  #                          "TCAGGATCAATAGCAA-1", 5, "Haplotype 2",
  #                          "TATCGAGGTATCAGTC-1", 6, "Haplotype 2",
  #                          "GGGCATCGTTAGAGCG-1", 7, "Haplotype 2",
  #                          "CCTACGTTCTGGGCAC-1", 8, "Haplotype 1",
  #                          "GGAACCCTCCGCTGGA-1", 9, "Haplotype 2",
  #                          "GTAGCATCACCAACAT-1", 10, "Haplotype 1")

  good_barcodes <- tribble(~Barcode, ~barcode_rank, ~best_haplotype,

                           "CCCGGAACATGGAGAC-1", 1, "Haplotype 2",
                           "TCAGGATCAATAGCAA-1", 2, "Haplotype 2",
                           "TATCGAGGTATCAGTC-1", 3, "Haplotype 2",
                           "GGGCATCGTTAGAGCG-1", 4, "Haplotype 2",
                           "GGAACCCTCCGCTGGA-1", 5, "Haplotype 2",
                           "TGTAGTGAGCTGCCCA-1", 6, "Not Phased",
                           "CTGTGCTAGGAAGCCT-1", 7, "Not Phased",
                           "GAATGAATCTGAGAGG-1", 8, "Haplotype 1",
                           "CCTACGTTCTGGGCAC-1", 9, "Haplotype 1",
                           "GTAGCATCACCAACAT-1", 10, "Haplotype 1")

  good_variants <- tribble(~Variant, ~variant_rank,
                           "chr1:114625476:T:C", 1,
                           "chr1:114637384:A:G", 2,
                           "chr1:114671205:G:A", 3,
                           "chr1:114676462:T:A", 4,
                           "chr1:114684872:G:A", 5,
                           "chr1:114686258:T:G", 6,
                           "chr1:114690425:C:T", 7,
                           "chr1:114692389:C:A", 8,
                           "chr1:114713909:G:T", 9,
                           "chr1:114716124:C:G", 10)

  n_variants <- good_variants %>% nrow()
  n_barcodes <- good_barcodes %>% nrow()

  horizontal_lines <- bcv %>%
    filter(Barcode %in% good_barcodes$Barcode,
           Variant %in% good_variants$Variant) %>%
    left_join(good_barcodes, by = "Barcode") %>%
    left_join(good_variants, by = "Variant") %>%
    group_by(Barcode, best_haplotype) %>%
    summarize(x_min = min(variant_rank), x_max = max(variant_rank),
              y_min = min(barcode_rank), y_max = max(barcode_rank)) %>%
    select(Barcode, x_min, x_max, y_min, y_max, best_haplotype) %>%
    ungroup()

  plot_data <- bcv %>%
    filter(Barcode %in% good_barcodes$Barcode,
           Variant %in% good_variants$Variant) %>%
    droplevels() %>%
    complete(Barcode, Variant) %>%
    left_join(good_barcodes, by = "Barcode") %>%
    left_join(horizontal_lines, by = "Barcode") %>%
    left_join(good_variants, by = "Variant") %>%
    mutate(allele_color = case_when(is.na(Allele) ~ "Missing Information",
                                    Allele == 0 ~ "Reference Allele",
                                    Allele == 1 ~ "Alternate Allele")) %>%
    mutate(allele_color = factor(allele_color,
                                 levels = c("Reference Allele",
                                            "Alternate Allele",
                                            "Missing Information"),
                                 ordered = TRUE)) %>%
    mutate(allele_size = case_when(variant_rank < x_min | variant_rank > x_max ~ 1,
                                   is.na(Allele) ~ 1.5,
                                   TRUE ~ 3)) %>%
    select(Barcode, Variant, barcode_rank, variant_rank, allele_color, allele_size)

  plot_data <- add_row(plot_data, Barcode = "Consensus H1", Variant = good_variants$Variant, barcode_rank = -1, variant_rank = good_variants$variant_rank, allele_color = c("Alternate Allele", "Alternate Allele", "Reference Allele", "Alternate Allele", "Alternate Allele", "Alternate Allele", "Alternate Allele", "Alternate Allele", "Reference Allele", "Reference Allele"), allele_size = 3)

  plot_data <- add_row(plot_data, Barcode = "Consensus H2 - No Mutation", Variant = good_variants$Variant, barcode_rank = -2, variant_rank = good_variants$variant_rank, allele_color = c("Reference Allele", "Reference Allele", "Alternate Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele"), allele_size = 3)

  plot_data <- add_row(plot_data, Barcode = "Consensus H2 - Q61 Mutation", Variant = good_variants$Variant, barcode_rank = -2.5, variant_rank = good_variants$variant_rank, allele_color = c("Reference Allele", "Reference Allele", "Alternate Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Alternate Allele", "Reference Allele"), allele_size = 3)

  plot_data <- add_row(plot_data, Barcode = "Consensus H2 - G13 Mutation", Variant = good_variants$Variant, barcode_rank = -3, variant_rank = good_variants$variant_rank, allele_color = c("Reference Allele", "Reference Allele", "Alternate Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Reference Allele", "Alternate Allele"), allele_size = 3)

  horizontal_lines <- add_row(horizontal_lines, Barcode = "Consensus H1", x_min = 1, x_max = 10, y_min = -1, y_max = -1, best_haplotype = "Haplotype 1")
  horizontal_lines <- add_row(horizontal_lines, Barcode = "Consensus H2 - No Mutation", x_min = 1, x_max = 10, y_min = -2, y_max = -2, best_haplotype = "Haplotype 2")
  horizontal_lines <- add_row(horizontal_lines, Barcode = "Consensus H2 - Q61 Mutation", x_min = 1, x_max = 10, y_min = -2.5, y_max = -2.5, best_haplotype = "Haplotype 2")
  horizontal_lines <- add_row(horizontal_lines, Barcode = "Consensus H2 - G13 Mutation", x_min = 1, x_max = 10, y_min = -3, y_max = -3, best_haplotype = "Haplotype 2")

  p <- ggplot() +
    theme_bw() +
    geom_rect(aes(xmin = 0.5, xmax = 12, ymin = 7.5, ymax = 10.5),
              alpha = 0.25, color = NA, fill = "#4292c6") +
    geom_rect(aes(xmin = 0.5, xmax = 12, ymin = 0.5, ymax = 5.5),
              alpha = 0.25, color = NA, fill = "#ef3b2c") +
    geom_vline(xintercept = seq(1:n_variants), linetype = 1, color = "#f0f0f0", size = 0.5) +
    geom_segment(data = horizontal_lines, aes(y = y_min, yend = y_max, x = x_min, xend = x_max), color = "white", size = 2.5) +
    geom_point(data = plot_data %>% filter(allele_size > 1), aes(x = variant_rank, y = barcode_rank, size = allele_size, fill = allele_color), color = "white", shape = 21, stroke = 1.25) +
    geom_segment(data = horizontal_lines, aes(y = y_min, yend = y_max, x = x_min, xend = x_max, color = best_haplotype), size = 1) +
    geom_point(data = plot_data, aes(x = variant_rank, y = barcode_rank, size = allele_size, fill = allele_color), shape = 21, stroke = 0) +
    coord_equal(ratio = 1) +
    scale_x_continuous(breaks = seq(1, n_variants + 1),
                       labels = c(good_variants$Variant, "N barcodes observed consistent\nwith this NRAS mutation pattern"),
                       expand = c(0, 3)) +
    scale_y_continuous(breaks = seq(-4, n_barcodes + 1),
                       labels = c("Consensus H2 - NRAS G13 Mutation", "Consensus H2 - NRAS Q61 Mutation", "Consensus H2 - No NRAS Mutation", "Consensus H1 - No NRAS Mutation", "", good_barcodes$Barcode, "")) +
    theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5)) +

    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(hjust = 0)) +

    scale_color_manual(values = c("#4292c6", "#ef3b2c", "#bdbdbd")) +
    scale_fill_manual(values = c("#762a83", "#c2a5cf", "#bdbdbd")) +
    labs(color = "Haplotype", fill = "Allele") +
    guides(size = FALSE) +
    annotate("text", x = 11, y = 10, label = 9) +
    annotate("text", x = 11, y = 9, label = 10) +
    annotate("text", x = 11, y = 8, label = 1) +
    annotate("text", x = 11, y = 7, label = 1) +
    annotate("text", x = 11, y = 6, label = 1) +
    annotate("text", x = 11, y = 5, label = 15) +
    annotate("text", x = 11, y = 4, label = 4) +
    annotate("text", x = 11, y = 3, label = 6) +
    annotate("text", x = 11, y = 2, label = 11) +
    annotate("text", x = 11, y = 1, label = 1) +
    annotate("text", x = 10.3, y = -1, label = "No NRAS mutation", hjust = 0,
             color = "#bdbdbd", fontface = "italic") +
    annotate("text", x = 10.3, y = -2, label = "No NRAS mutation", hjust = 0,
             color = "#bdbdbd", fontface = "italic") +
    annotate("text", x = 10.3, y = -2.5, label = "NRAS Q61 mutation", hjust = 0,
             color = "#bdbdbd", fontface = "italic") +
    annotate("text", x = 10.3, y = -3, label = "NRAS G13 mutation", hjust = 0,
             color = "#bdbdbd", fontface = "italic") +

    geom_hline(yintercept = 0, lty = 2) +
    geom_label(aes(x = 0.5, y = 0, label = "Consensus Haplotypes"), hjust = 0,
               fontface = "bold") +
    geom_label(aes(x = 0.5, y = 11, label = "Observed Patterns of Linked-Reads"),
               hjust = 0, fontface = "bold") +

    geom_text(aes(x = 12, y = 3, label = "Haplotype 2"),
              color = "#ef3b2c", angle = -90, vjust = 1, nudge_x = -0.1) +
    geom_text(aes(x = 12, y = 9, label = "Haplotype 1"),
              color = "#4292c6", angle = -90, vjust = 1, nudge_x = -0.1) +

    geom_label(aes(x = 9.5, y = -0.5, label = "Consensus H1"),
               color = "white", fill = "#4292c6", label.padding = unit(0.2, "lines")) +
    geom_label(aes(x = 9.5, y = -1.5, label = "Consensus H2"),
               color = "white", fill = "#ef3b2c", label.padding = unit(0.2, "lines"))

  q <- p + theme(axis.text = element_blank()) + guides(fill = FALSE, color = FALSE)

  #ggsave("main_figures/nras_haplotype.pdf", p, width = 10, height = 10, useDingbats = FALSE)
  #ggsave("main_figures/nras_haplotype.no_legend_or_labels.pdf", q, width = 7, height = 10, useDingbats = FALSE)

  ggsave("~/Desktop/nras_haplotype.3.pdf", p, width = 10, height = 10, useDingbats = FALSE)
  ggsave("~/Desktop/nras_haplotype.no_legend_or_labels.3.pdf", q, width = 7, height = 10, useDingbats = FALSE)

  # r <- good_variants %>%
  #   separate(Variant, into = c("chr", "pos", "ref", "alt"), sep = ":") %>%
  #   ggplot(aes(x = as.numeric(pos), y = 0)) +
  #   geom_segment(x = 114704469, xend = 114716894, y = 0, yend = 0, color = "#c6dbef", lwd = 10) +
  #   geom_hline(yintercept = 0, color = "#c6dbef") +
  #   annotate("text", x = (114704469 + 114716894)/2, y = 0.125, label = "NRAS", vjust = 0) +
  #   geom_point() +
  #   theme_bw() +
  #   theme(axis.line = element_line(colour = "black"),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border = element_blank(),
  #         panel.background = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.line.y = element_blank(),
  #         axis.line.x = element_blank()) +
  #   scale_x_continuous(position = "top",
  #                      breaks = seq(114600000, 114725000, 25000),
  #                      labels = c("114,600,000", "114,625,000", "114,650,000", "114,675,000", "114,700,000", "114,725,000"),
  #                      limits = c(114600000, 114725000)) +
  #   scale_y_continuous(limits = c(-0.1, 0.2)) +
  #   labs(x = "Chromosome 1 Position")
  # ggsave("main_figures/nras_positions.pdf", r, width = 7, height = 1, useDingbats = FALSE)

}

rm(plot_barcode_variants, main, supp)
