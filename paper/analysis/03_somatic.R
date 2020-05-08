################################################################################
# Somatic phasing
# Concordance with longranger
# Mutations on phase sets
# Mutation pairs per phase set
################################################################################

main = "figures/03_somatic_phasing/main/"
supp = "figures/03_somatic_phasing/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

manuscript_numbers[["03_somatic"]] <- list()

manuscript_numbers[["03_somatic"]][["reads_overlapping_mutation_site_assigned_haplotype_proportion"]] <- phasing_variants_mapq20_tbl %>% summarize(a = sum(barcode_REF_H1), b = sum(barcode_REF_H2), c = sum(barcode_ALT_H1), d = sum(barcode_ALT_H2), abcd = a + b + c + d, e = sum(just_barcode_not_phased), abcde = abcd + e, overall = abcd/abcde) %>% pull(overall)

manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_maf"]] <- maf_tbl %>% filter(Variant_Type == "SNP") %>% nrow()
manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_maf_per_sample"]] <- manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_maf"]]/6
manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping"]] <- phasing_variants_mapq20_tbl %>% nrow()
manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping_per_sample"]] <- phasing_variants_mapq20_tbl %>% nrow()/6
manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping_enough_coverage"]] <- phasing_variants_mapq20_tbl %>% filter(enough_coverage) %>% nrow()
manuscript_numbers[["03_somatic"]][["pct_somatic_mutations_from_10Xmapping_enough_coverage"]] <- 100*manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping_enough_coverage"]]/manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping"]]
manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping_enough_coverage_by_sample"]] <- phasing_variants_mapq20_tbl %>% filter(enough_coverage) %>% group_by(sample) %>% summarize(n())

# Precision/recall balance threshold to maximize longranger concordance
{
  proportions <- seq(0.8, 1, by = 0.01)
  tp <- rep(NA, length(proportions))
  tn <- rep(NA, length(proportions))
  fp <- rep(NA, length(proportions))
  fn <- rep(NA, length(proportions))

  for (i in 1:length(proportions)) {
    tp[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(phased_by_linked_alleles != "NC" & phased_by_barcodes != "NC") %>%
      filter( (Genotype == "1|0" & (pct_ALT_on_H1 >= proportions[i] | phased_by_barcodes == "H1")) |
                (Genotype == "0|1" & (pct_ALT_on_H2 >= proportions[i] | phased_by_barcodes == "H2"))) %>%
      nrow()

    fp[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(phased_by_linked_alleles != "NC" & phased_by_barcodes != "NC") %>%
      filter( (Genotype == "1|0" & (pct_ALT_on_H2 >= proportions[i] | phased_by_barcodes == "H2")) |
                (Genotype == "0|1" & (pct_ALT_on_H1 >= proportions[i] | phased_by_barcodes == "H1"))) %>%
      nrow()

    fn[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(phased_by_linked_alleles != "NC" & phased_by_barcodes != "NC") %>%
      filter( pct_ALT_on_H1 < proportions[i] & pct_ALT_on_H2 < proportions[i] & phased_by_barcodes == "NP") %>%
      nrow()
  }

  prec = tp/(tp + fp)
  rec = tp/(tp + fn)

  ggplot(tibble(proportions, prec, rec), aes(x = prec, y = rec)) +
    geom_point() +
    geom_label_repel(aes(label = proportions), size = 3) +
    expand_limits(x = c(.99,1), y = c(.9,1)) +
    labs(x = "Precision", y = "Recall") +
    theme_bw() +
    #coord_equal() +
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
    ggsave(str_c(supp, "precision_recall.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)

  manuscript_numbers[["03_somatic"]][["precision_at_phase_proportion"]] <- prec[which(proportions == phase_proportion)]
  manuscript_numbers[["03_somatic"]][["recall_at_phase_proportion"]] <- rec[which(proportions == phase_proportion)]
  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased"]] <- phasing_variants_mapq20_tbl %>% filter(phased == "Phased") %>% nrow()
  manuscript_numbers[["03_somatic"]][["pct_somatic_mutations_with_enough_coverage_phased"]] <- 100*manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased"]]/manuscript_numbers[["03_somatic"]][["n_somatic_mutations_from_10Xmapping_enough_coverage"]]

  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_0/1"]] <- phasing_variants_mapq20_tbl %>% filter(enough_coverage, Genotype == "0/1") %>% nrow()
  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_0/1_phased"]] <- phasing_variants_mapq20_tbl %>% filter(Genotype == "0/1", phased == "Phased") %>% nrow()
  manuscript_numbers[["03_somatic"]][["pct_somatic_mutations_with_enough_coverage_0/1_phased"]] <- 100*manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_0/1_phased"]]/manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_0/1"]]

  # compare with barcode method
  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased_table"]] <- phasing_variants_mapq20_tbl %>% select(phased_by_linked_alleles, phased_by_barcodes) %>% table()

  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_both_phased"]] <- sum(manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased_table"]][1:2, 1:2])
  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_both_phased_concordant"]] <- sum(diag(manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased_table"]])[1:2])
  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_both_phased_concordant_pct"]] <- 100*manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_both_phased_concordant"]]/manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_both_phased"]]

  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_added_by_bc"]] <- sum(manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased_table"]][3:4,1:2])
  manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_added_by_la"]] <- sum(manuscript_numbers[["03_somatic"]][["n_somatic_mutations_with_enough_coverage_phased_table"]][1:2,3:4])



  rm(proportions, tp, tn, fp, fn, prec, rec, i)
}

# Concordance with longranger
{
  phasing_variants_mapq20_tbl %>%
    #filter(phased_by != "NC") %>%
    #filter(!is.na(Genotype)) %>%
    replace_na(list(pct_ALT_on_H1 = 0.5)) %>%
    replace_na(list(Genotype = "Missing")) %>%
    mutate(phasing_pair = str_c(phased_by_linked_alleles, phased_by_barcodes, sep = "/")) %>%
    select(Genotype, phasing_pair, phased_by) %>%
    mutate(Genotype = factor(Genotype),
           phasing_pair = factor(phasing_pair),
           phased_by = factor(phased_by,
                              levels = c("Both (agree)", "BC", "LA", "LR", "Both (conflict)", "NC", "NP"),
                              labels = c("Both (agree)", "Barcodes", "Linked Alleles", "Linked Alleles", "Both (conflict)", "Not Enough Coverage", "Not Phased"),
                              ordered = TRUE)) %>%
    group_by(Genotype, phasing_pair, phased_by) %>%
    count(Genotype, phasing_pair, phased_by) %>%
    mutate(my_color = case_when(phasing_pair == "NC/NC" & Genotype == "Missing" ~ "white",
                                TRUE ~ "black")) %>%
    ggplot(aes(x = Genotype, y = fct_rev(phasing_pair), fill = n)) +
    geom_tile(color = "black", show.legend = FALSE) +
    geom_text(aes(label = n, color = my_color), size = 8/ggplot2:::.pt) +
    facet_wrap(~ phased_by, ncol = 1, scales = "free_y") +
    labs(x = "Genotype (Long Ranger)", y = "Somatic Haplotype Phasing Call (Linked Alleles / Barcodes)") +
    scale_fill_gradient2(low = "white", high = "#969696") +
    scale_color_identity() +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "longranger_concordance.pdf"),
           width = 3.5, height = 7.25, useDingbats = FALSE)

  n_nc <- phasing_variants_mapq20_tbl %>%
    filter(phased_by_linked_alleles == "NC" & phased_by_barcodes == "NC") %>%
    nrow()

  phasing_variants_mapq20_tbl %>%
    select(phased_by_linked_alleles, phased_by_barcodes) %>%
    group_by(phased_by_barcodes, phased_by_linked_alleles) %>%
    summarize(count = n()) %>%
    filter(!(phased_by_linked_alleles == "NC" & phased_by_barcodes == "NC")) %>%
    mutate(my_color = case_when(phased_by_linked_alleles == phased_by_barcodes ~ "white",
                                TRUE ~ "black")) %>%
    ggplot(aes(x = phased_by_linked_alleles,
               y = phased_by_barcodes,
               fill = count)) +
    geom_tile(color = "black", show.legend = FALSE) +
    geom_text(aes(label = count, color = my_color), size = 8/ggplot2:::.pt) +
    annotate("text", x = "NC", y = "NC", label = str_c("(", n_nc, ")"), size = 8/ggplot2:::.pt) +
    labs(x = "Linked Alleles Method", y = "Barcodes Method") +
    scale_fill_gradient2(low = "white", high = "#969696") +
    #scale_fill_gradient2(low = "white", high = "#000000") +
    scale_color_identity() +
    coord_equal() +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "methods_phasing.pdf"),
           width = 2.25, height = 2.25, useDingbats = FALSE)

  rm(n_nc)
}

# Phased somatic mutations on haplotypes
{
  sample_id <- "59114_4"
  chromosome_of_interest <- "chr1"
  min_length = 1e3

  plot_df <- phase_sets_tbl %>%
    filter(sample == sample_id) %>%
    filter(length_variants > min_length & !is.na(length_variants)) %>%
    filter(chromosome == chromosome_of_interest) %>%
    droplevels() %>%
    mutate(alternating_01 = row_number() %% 2)

  plot_df1 <- plot_df %>% mutate(y_min = 1.1, y_max = 1.9,
                                 color_vector = case_when(alternating_01 == 0 ~ my_color_100,
                                                          TRUE ~ my_color_50))
  plot_df2 <- plot_df %>% mutate(y_min = 0.1, y_max = 0.9,
                                 color_vector = case_when(alternating_01 == 0 ~ my_color_100,
                                                          TRUE ~ my_color_50))

  plot_df_together <- bind_rows(plot_df1, plot_df2) %>% arrange(ps_id)

  mut_df <- phasing_variants_mapq20_tbl %>%
    filter(sample == sample_id, Chromosome == chromosome_of_interest) %>%
    mutate(haplotype_of_variant_shape = case_when(
      !enough_coverage ~ "Not Enough Coverage",
      phased == "Phased" & phased_by == "Both (agree)" & phased_by_linked_alleles == "H1" ~ "Haplotype 1",
      phased == "Phased" & phased_by == "LA" & phased_by_linked_alleles == "H1" ~ "Haplotype 1",
      phased == "Phased" & phased_by == "BC" & phased_by_barcodes == "H1" ~ "Haplotype 1",
      phased == "Phased" & phased_by == "Both (agree)" & phased_by_linked_alleles == "H2" ~ "Haplotype 2",
      phased == "Phased" & phased_by == "LA" & phased_by_linked_alleles == "H2" ~ "Haplotype 2",
      phased == "Phased" & phased_by == "BC" & phased_by_barcodes == "H2" ~ "Haplotype 2",
      phased_by == "Both (conflict)" ~ "Not Phased",
      TRUE ~ "Not Phased")) %>%
    mutate(haplotype_of_variant_yaxis = case_when(haplotype_of_variant_shape == "Haplotype 1" ~ 0.1,
                                                  haplotype_of_variant_shape == "Haplotype 2" ~ 1.9,
                                                  TRUE ~ 1)) %>%
    filter(haplotype_of_variant_shape != "Not Enough Coverage")

  ggplot(plot_df_together) +
    geom_segment(data = chromosome_length_tbl %>%
                   filter(contig == chromosome_of_interest),
                 aes(y = 1.5, yend = 1.45, x = 0, xend = size/1e6)) +
    geom_segment(data = chromosome_length_tbl %>%
                   filter(contig == chromosome_of_interest),
                 aes(y = 0.5, yend = 0.55, x = 0, xend = size/1e6)) +
    geom_rect(aes(xmin = first_variant_pos/1e6,
                  xmax = last_variant_pos/1e6,
                  ymin = y_min,
                  ymax = y_max),
              fill = "#ffffff",
              show.legend = FALSE) +
    geom_rect(aes(xmin = first_variant_pos/1e6,
                  xmax = last_variant_pos/1e6,
                  ymin = y_min,
                  ymax = y_max,
                  fill = color_vector),
              show.legend = FALSE) +
    geom_point(data = mut_df, aes(x = Position/1e6,
                                  y = haplotype_of_variant_yaxis,
                                  shape = factor(haplotype_of_variant_shape)),
               fill = "#000000",
               alpha = 0.5,
               size = 2,
               show.legend = FALSE) +
    #geom_text(aes(x = -Inf, y = 2.5, label = "Phased Somatic Mutation"), size = 3) +
    labs(x = str_c(chromosome_of_interest, " Position (Mb)"),
         y = NULL,
         shape = "Somatic Variant\nHaplotype") +
    scale_y_continuous(expand = c(0.05, 0.05), limits = c(0, 2.5)) +
    scale_x_continuous(expand = c(0.03, 0.03)) +
    #scale_fill_manual(values = c("#a6cee3", "#fb9a99", "#1f78b4", "#e31a1c")) +
    scale_fill_identity() +
    scale_shape_manual(values = c(17, 25, 16)) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
  ggsave(str_c(main, "somatic_mutations_on_phase_sets.pdf"),
         width = 4.75, height = 1.125, useDingbats = FALSE)

  rm(plot_df, plot_df1, plot_df2, plot_df_together, mut_df,
     min_length, sample_id, chromosome_of_interest)

}

# variants per phase set
{

  phasing_variants_grouped <- phasing_variants_mapq20_tbl %>%
    filter(enough_coverage, Phase_Set_Length >= 1e3) %>%
    group_by(sample, Phase_Set, Phase_Set_Length) %>%
    summarize(n_somatic_variants = n(),
              n_phased = sum(phased == "Phased")) %>%
    mutate(somatic_mutations_per_Mb = 1e6 * n_somatic_variants/Phase_Set_Length,
           proportion_variants_phased = n_phased / n_somatic_variants,
           n_pairs_phased = choose(n_phased, 2))

  manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb"]] <- phasing_variants_grouped %>% filter(Phase_Set_Length >= 1e3) %>% nrow()
  manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb_0_somatic_mutation_pairs"]] <- phasing_variants_grouped %>% filter(Phase_Set_Length >= 1e3) %>% filter(n_pairs_phased == 0) %>% nrow()
  manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb_0_somatic_mutations"]] <- phasing_variants_grouped %>% filter(Phase_Set_Length >= 1e3) %>% filter(n_phased == 0) %>% nrow()
  manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb_1_somatic_mutations"]] <- phasing_variants_grouped %>% filter(Phase_Set_Length >= 1e3) %>% filter(n_phased == 1) %>% nrow()

  manuscript_numbers[["03_somatic"]][["pct_phase_sets_1kb_0_somatic_mutation_pairs"]] <- 100*manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb_0_somatic_mutation_pairs"]]/manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb"]]
  manuscript_numbers[["03_somatic"]][["pct_phase_sets_1kb_0_somatic_mutations"]] <- 100*manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb_0_somatic_mutations"]]/manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb"]]
  manuscript_numbers[["03_somatic"]][["pct_phase_sets_1kb_1_somatic_mutation"]] <- 100*manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb_1_somatic_mutations"]]/manuscript_numbers[["03_somatic"]][["n_phase_sets_1kb"]]

  plot_data <- phasing_variants_grouped %>%
    filter(Phase_Set_Length >= 1e3, n_pairs_phased >= 1) %>%
    mutate(n_pairs_color = case_when(n_pairs_phased == 1 ~ "1 pair",
                                     n_pairs_phased <= 3 ~ "<= 3",
                                     n_pairs_phased <= 10 ~ "<= 10",
                                     n_pairs_phased <= 100 ~ "<= 100",
                                     TRUE ~ "> 100")) %>%
    mutate(n_pairs_color = factor(n_pairs_color,
                                  levels = c("0", "1 pair", "<= 3", "<= 10", "<= 100", "> 100"),
                                  ordered = TRUE)) %>%
    arrange(n_pairs_color)

  n_phase_sets <- plot_data %>% nrow()
  manuscript_numbers[["03_somatic"]][["n_phase_sets_gt_1kb_1_pair_by_n_pairs"]] <- plot_data %>% pull(n_pairs_color) %>% table()/n_phase_sets
  n_phase_sets_1.0 <- plot_data %>% filter(proportion_variants_phased == 1) %>% nrow()
  n_phase_sets_0.75 <- plot_data %>% filter(proportion_variants_phased >= 0.75) %>% nrow()
  manuscript_numbers[["03_somatic"]][["n_phase_sets_gt_1kb_1_pair"]] <- n_phase_sets
  manuscript_numbers[["03_somatic"]][["n_phase_sets_gt_1kb_1_pair_all_phased"]] <- n_phase_sets_1.0
  manuscript_numbers[["03_somatic"]][["n_phase_sets_gt_1kb_1_pair_0.75_phased"]] <- n_phase_sets_0.75
  manuscript_numbers[["03_somatic"]][["pct_phase_sets_all_variants_phased"]] <- 100*n_phase_sets_1.0/n_phase_sets
  manuscript_numbers[["03_somatic"]][["pct_phase_sets_all_variants_phased"]] <- 100*n_phase_sets_1.0/n_phase_sets
  manuscript_numbers[["03_somatic"]][["pct_phase_sets_0.75_variants_phased"]] <- 100*n_phase_sets_0.75/n_phase_sets

  min_log2 <- plot_data %>% pull(somatic_mutations_per_Mb) %>%
    log2() %>% min() %>% round(digits = 0)
  max_log2 <- plot_data %>% pull(somatic_mutations_per_Mb) %>%
    log2() %>% max() %>% round(digits = 0)
  my_breaks <- seq(from = min_log2, to = max_log2, by = 1)

  manuscript_numbers[["03_somatic"]][["min_log2_somatic_mutations_per_Mb"]] <- min_log2
  manuscript_numbers[["03_somatic"]][["max_log2_somatic_mutations_per_Mb"]] <- max_log2
  manuscript_numbers[["03_somatic"]][["min_somatic_mutations_per_Mb"]] <- 2^min_log2
  manuscript_numbers[["03_somatic"]][["max_somatic_mutations_per_Mb"]] <- 2^max_log2
  manuscript_numbers[["03_somatic"]][["median_somatic_mutations_per_Mb"]] <- plot_data %>% pull(somatic_mutations_per_Mb) %>% summary()

  p <- ggplot(data = plot_data,
              aes(x = log2(somatic_mutations_per_Mb),
                  y = proportion_variants_phased,
                  color = n_pairs_color)) +
    geom_point(shape = 16, alpha = 0.75) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = my_breaks) +
    labs(x = "Somatic Mutations per Mb (log2)",
         y = "Proportion of Mutations Phased",
         color = "Pairs of Phased\nSomatic Mutations",
         size = "Pairs of Phased\nSomatic Mutations") +
    scale_color_brewer(palette = "YlGnBu", drop = FALSE) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"))

  q_with_legend <- ggExtra::ggMarginal(p, type = "histogram", groupFill = TRUE, color = NA)
  q <- ggExtra::ggMarginal(p + guides(color = FALSE, size = FALSE),
                           type = "histogram",
                           groupFill = TRUE,
                           color = NA)

  ggsave(str_c(main, "somatic_mutations_per_Mb.with_legend.pdf"), q_with_legend,
         width = 4.75, height = 1.5*2.25, useDingbats = FALSE)

  ggsave(str_c(main, "somatic_mutations_per_Mb.pdf"), q,
         width = 4.75, height = 1.5*2.25, useDingbats = FALSE)

  rm(phasing_variants_grouped, plot_data, p, q, q_with_legend,
     n_phase_sets, n_phase_sets_1.0, n_phase_sets_0.75,
     min_log2, max_log2, my_breaks)
}

# looking at driver mutations
{
  drivers <- driver_mutations_tbl %>%
    rowwise() %>%
    mutate(p1 = str_sub(str_split(protein, pattern = "\\.", simplify = TRUE)[2], 1, 1),
           p2 = str_sub(str_split(protein, pattern = "\\.", simplify = TRUE)[2], -1),
           synonymous = (p1 == p2),
           Variant = str_c(chr, pos, ref, alt, sep = ":"))

  tumor_vafs <- driver_mutations_vaf_tbl %>%
    filter(timepoint != "Normal") %>%
    mutate(variant_key = str_c(chr, pos, ref, alt, sep = ":"))
  normal_vafs <- driver_mutations_vaf_tbl %>%
    filter(timepoint == "Normal")  %>%
    mutate(variant_key = str_c(chr, pos, ref, alt, sep = ":")) %>%
    mutate(normal_vaf = vaf,
           normal_alleles = alleles)

  plot_df <- sombx_driver_mapq20_tbl %>%
    mutate(Variant = str_c(chromosome, position, ref, alt, sep = ":")) %>%
    left_join(drivers, by = c("Variant", "ref", "alt")) %>%
    filter(!synonymous) %>%
    rowwise() %>%
    mutate(n_ref_h1_bx = case_when(is.na(ref_barcodes_H1) ~ 0,
                                   TRUE ~ as.double(length(str_split(ref_barcodes_H1, ";", simplify = TRUE)))),
           n_ref_h2_bx = case_when(is.na(ref_barcodes_H2) ~ 0,
                                   TRUE ~ as.double(length(str_split(ref_barcodes_H2, ";", simplify = TRUE)))),
           n_ref_none_bx = case_when(is.na(ref_barcodes_None) ~ 0,
                                   TRUE ~ as.double(length(str_split(ref_barcodes_None, ";", simplify = TRUE)))),
           n_alt_h1_bx = case_when(is.na(alt_barcodes_H1) ~ 0,
                                   TRUE ~ as.double(length(str_split(alt_barcodes_H1, ";", simplify = TRUE)))),
           n_alt_h2_bx = case_when(is.na(alt_barcodes_H2) ~ 0,
                                   TRUE ~ as.double(length(str_split(alt_barcodes_H2, ";", simplify = TRUE)))),
           n_alt_none_bx = case_when(is.na(alt_barcodes_None) ~ 0,
                                   TRUE ~ as.double(length(str_split(alt_barcodes_None, ";", simplify = TRUE)))),
           total_ref = sum(n_ref_h1_bx, n_ref_h2_bx, n_ref_none_bx),
           total_alt = sum(n_alt_h1_bx, n_alt_h2_bx, n_alt_none_bx),
           barcode_vaf = total_alt/(total_alt + total_ref)) %>%
    left_join(tumor_vafs %>% select(variant_key, sample, vaf, alleles),
              by = c("variant_key", "sample")) %>%
    left_join(normal_vafs %>% select(variant_key, patient, normal_vaf, normal_alleles),
              by = c("variant_key", "patient")) %>%
    filter(!is.na(normal_vaf),
           normal_vaf <= 0.02,
           vaf >= 0.05,
           nchar(alleles) >= 10,
           n_alt_h1_bx + n_alt_h2_bx > 0) %>%
    left_join(phasing_variants_driver_mapq20_tbl %>%
                select(-c("patient", "timepoint", "vcf_column", "sorted",
                          "cnv_maf_status", "display_name", "sample_n",
                          "my_color_100", "my_color_75", "my_color_50",
                          "my_color_25", "my_shape")),
              by = c("variant_key" = "Variant", "sample")) %>%
    mutate(variant_name = str_c(gene, protein, sep = " ")) %>%
    filter(enough_coverage) %>%
    mutate(phased_by_SH = case_when(phased == "Phased" & phased_by == "Both (agree)" & phased_by_linked_alleles == "H1" ~ "H1",
                                    phased == "Phased" & phased_by == "LA" & phased_by_linked_alleles == "H1" ~ "H1",
                                    phased == "Phased" & phased_by == "BC" & phased_by_barcodes == "H1" ~ "H1",
                                    phased == "Phased" & phased_by == "Both (agree)" & phased_by_linked_alleles == "H2" ~ "H2",
                                    phased == "Phased" & phased_by == "LA" & phased_by_linked_alleles == "H2" ~ "H2",
                                    phased == "Phased" & phased_by == "BC" & phased_by_barcodes == "H2" ~ "H2",
                                    phased_by == "Both (conflict)" ~ "NP",
                                    TRUE ~ "NP")) %>%
    replace_na(list(Genotype = "NA")) %>%
    filter(gene %in% important_mutations_tbl$gene) %>%
    mutate(plot_category = case_when(phased_by_SH == "NP" & Genotype == "0/1" ~ "Not phased",
                                     phased_by_SH == "NP" & Genotype == "1|0" ~ "Not phased",
                                     phased_by_SH == "NP" & Genotype == "0|1" ~ "Not phased",
                                     phased_by_SH == "NP" & Genotype == "NA" ~ "Not phased",
                                     phased_by_SH == "H1" & Genotype == "0/1" ~ "Phase added",
                                     phased_by_SH == "H1" & Genotype == "1|0" ~ "Phase consistent",
                                     phased_by_SH == "H1" & Genotype == "0|1" ~ "Phase inconsistent",
                                     phased_by_SH == "H1" & Genotype == "NA" ~ "Phased non-call",
                                     phased_by_SH == "H2" & Genotype == "0/1" ~ "Phase added",
                                     phased_by_SH == "H2" & Genotype == "1|0" ~ "Phase inconsistent",
                                     phased_by_SH == "H2" & Genotype == "0|1" ~ "Phase consistent",
                                     phased_by_SH == "H2" & Genotype == "NA" ~ "Phased non-call")) %>%
    mutate(plot_category = factor(plot_category,
                                  levels = c("Phase added", "Phase consistent",
                                             "Phase inconsistent", "Phased non-call",
                                             "Not phased"), ordered = TRUE))

  correlation_value <- round(as.numeric(cor.test(plot_df$vaf, plot_df$barcode_vaf)$estimate), 2)

  ggplot(plot_df, aes(x = vaf, y = barcode_vaf)) +
    geom_abline(lty = 2) +
    geom_smooth(method = "lm") +
    geom_point() +
    coord_equal() +
    expand_limits(x = c(0,1), y = c(0,1)) +
    annotate("text", x = 1, y = 0, label = correlation_value, size = 2) +
    labs(x = "Variant Allele Frequency (calculated from reads)",
         y = "Variant Barcode Frequency (calculated from barcodes)") +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "vaf_correlation.pdf"),
           height = 3.5, width = 3.5, useDingbats = FALSE)

  p <- ggplot(plot_df, aes(y = variant_name, x = display_name)) +
    geom_point(aes(shape = plot_category), fill = "black", size = 3, stroke = 0.75) +
    scale_shape_manual(values = c(9, 23, 13, 5, 4), drop = FALSE) +
    labs(shape = "GT", color = "Hap") +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 270, size = 8, hjust = 0, vjust = 0.5),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines"),
          legend.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.key = element_blank())

  ggsave(str_c(main, "mm_mutations.with_legend.pdf"), p,
         height = 3.375, width = 2.25, useDingbats = FALSE)
  ggsave(str_c(main, "mm_mutations.pdf"), p + guides(color = FALSE, shape = FALSE),
         height = 3.375, width = 2.25, useDingbats = FALSE)

  manuscript_numbers[["03_somatic"]][["ATR_pct_ALT_on_H1"]] <- plot_df %>% filter(gene == "ATR") %>% select(pct_ALT_on_H1, pct_ALT_on_H2)

  rm(drivers, plot_df, correlation_value, tumor_vafs, normal_vafs, p)

}

rm(main, supp)
