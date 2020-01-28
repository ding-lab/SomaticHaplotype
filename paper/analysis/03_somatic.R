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

# Precision/recall balance threshold to maximize longranger concordance
{
  proportions <- seq(0, 0.1, by = 0.01)
  tp <- rep(NA, length(proportions))
  tn <- rep(NA, length(proportions))
  fp <- rep(NA, length(proportions))
  fn <- rep(NA, length(proportions))

  for (i in 1:length(proportions)) {
    tp[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(n_ALT_H1 + n_ALT_H2 >= 10) %>%
      filter(pct_ALT_on_H1 <= proportions[i] & Genotype == "0|1") %>%
      nrow()

    fp[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(n_ALT_H1 + n_ALT_H2 >= 10) %>%
      filter(pct_ALT_on_H1 <= proportions[i] & Genotype == "1|0") %>%
      nrow()

    fn[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(n_ALT_H1 + n_ALT_H2 >= 10) %>%
      filter(pct_ALT_on_H1 > proportions[i] & Genotype == "0|1") %>%
      nrow()

    tn[i] <- phasing_variants_mapq20_tbl %>%
      filter(Genotype %in% c("0|1", "1|0")) %>%
      filter(n_ALT_H1 + n_ALT_H2 >= 10) %>%
      filter(pct_ALT_on_H1 > proportions[i] & Genotype == "1|0") %>%
      nrow()
  }

  prec = tp/(tp + fp)
  rec = tp/(tp + fn)

  tpr = tp/(tp + fn)
  fpr = fp/(fp + tn)

  ggplot(tibble(proportions, prec, rec), aes(x = prec, y = rec)) +
    geom_point() +
    geom_label_repel(aes(label = proportions), size = 3) +
    expand_limits(x = 1, y = 1) +
    labs(x = "Precision", y = "Recall") +
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
    ggsave(str_c(supp, "precision_recall.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)

  rm(proportions, tp, tn, fp, fn, prec, rec, tpr, fpr, i)
}

# Concordance with longranger
{
  H1_proportion <- 0.09

  plot_df <- phasing_variants_mapq20_tbl %>%
    filter(Genotype %in% c("0|1", "1|0")) %>%
    filter(!is.na(pct_ALT_on_H1)) %>%
    filter(n_ALT_H1 + n_ALT_H2 >= 10) %>%
    mutate(color_gt = case_when(
      pct_ALT_on_H1 <= H1_proportion & Genotype == "0|1" ~ "Agree",
      pct_ALT_on_H1 >= 1 - H1_proportion & Genotype == "1|0" ~ "Agree",
      pct_ALT_on_H1 <= H1_proportion & Genotype == "1|0" ~ "Disagree",
      pct_ALT_on_H1 >= 1 - H1_proportion & Genotype == "0|1" ~ "Disagree",
      TRUE ~ "No call"
    ))

  n_total <- plot_df %>% nrow()
  n_agree <- plot_df %>% filter(color_gt == "Agree") %>% nrow()
  n_disagree <- plot_df %>% filter(color_gt == "Disagree") %>% nrow()
  concordance <- n_agree/n_total
  discordance <- n_disagree/n_total
  print(c("Concordance: ", 100*concordance, n_agree, n_total))
  print(c("Discordance: ", 100*discordance, n_disagree, n_total))

  plot_df %>%
    ggplot(aes(x = Genotype, y = pct_ALT_on_H1*100)) +
    geom_jitter(aes(color = color_gt),
                shape = 16, height = 0, width = 0.25,
                alpha = 0.25, show.legend = FALSE) +
    geom_violin(fill = "white", draw_quantiles = 0.5) +
    geom_hline(yintercept = H1_proportion*100, linetype = 2) +
    geom_hline(yintercept = (1 - H1_proportion)*100, linetype = 2) +
    scale_y_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 10),
                       expand = c(0.01, 0.01)) +
    scale_color_manual(values = c("#3182bd", "#e31a1c", "#bdbdbd")) +
    annotate("text", x = "0|1", y = 93.5, label = "Phased as 1|0", size = 2) +
    annotate("text", x = "1|0", y = 6.5, label = "Phased as 0|1", size = 2) +
    labs(x = "Somatic Mutation Genotype (longranger)",
         y = "Linked Variants Assigned to H1 (%)") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "longranger_concordance.pdf"),
           width = 2.25, height = 2.25, useDingbats = FALSE)

  rm(plot_df, n_total, n_agree, n_disagree, concordance, discordance)
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
          n_ALT_H1 + n_ALT_H2 < 10 ~ "Not Enough Coverage",
          pct_ALT_on_H1 >= 0.91 ~ "Haplotype 1",
          pct_ALT_on_H2 >= 0.91 ~ "Haplotype 2",
          TRUE ~ "Not Phased")) %>%
        mutate(haplotype_of_variant_yaxis = case_when(
          n_ALT_H1 + n_ALT_H2 < 10 ~ 1,
          pct_ALT_on_H1 > 0.91 ~ 0.1,
          pct_ALT_on_H2 > 0.91 ~ 1.9,
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
    filter(n_ALT_H1 + n_ALT_H2 >= 10, Phase_Set_Length > 1e3) %>%
    group_by(sample, Phase_Set, Phase_Set_Length) %>%
    summarize(n_somatic_variants = n(),
              n_phased = sum(pct_ALT_on_H1 >= 0.91 | pct_ALT_on_H2 >= 0.91, na.rm = TRUE)) %>%
    mutate(somatic_mutations_per_Mb = 1e6 * n_somatic_variants/Phase_Set_Length,
           proportion_variants_phased = n_phased / n_somatic_variants,
           n_pairs_phased = choose(n_phased, 2))

  stat_data <- somatic_per_phase_set_mapq20_tbl %>%
    left_join(phasing_variants_grouped,
              by = c("sample", "ps_id" = "Phase_Set",
                     "length_variants" = "Phase_Set_Length")) %>%
    replace_na(list(somatic_mutations_per_Mb = 0)) %>%
    filter(length_variants >= 1e3)

  n_pb <- stat_data %>% nrow()
  n_pb_1_sm <- stat_data %>% filter(n_somatic_variants.x > 0) %>% nrow()
  tot_len <- stat_data %>% pull(length_variants) %>% sum()
  tot_sm <- stat_data %>% pull(n_somatic_variants.x) %>% sum()
  pb_1_sm <- 100*n_pb_1_sm/n_pb
  sm_mb <- 1e6*tot_sm/tot_len
  print(c("% phase blocks with > 0 somatic mutations:", pb_1_sm))
  print(c("Somatic mutations per Mb within phase blocks:", sm_mb))

  plot_data <- phasing_variants_grouped %>%
    filter(Phase_Set_Length >= 1e3, n_phased > 1) %>%
    mutate(n_pairs_color = case_when(n_pairs_phased == 1 ~ "1 pair",
                                     n_pairs_phased <= 3 ~ "< 3",
                                     n_pairs_phased <= 10 ~ "< 10",
                                     n_pairs_phased <= 100 ~ "< 100",
                                     TRUE ~ "> 100")) %>%
    mutate(n_pairs_color = factor(n_pairs_color,
                                  levels = c("0", "1 pair", "< 3", "< 10", "< 100", "> 100"),
                                  ordered = TRUE)) %>%
    arrange(n_pairs_color)

  print("With WGS only:")
  plot_data %>% pull(n_pairs_color) %>% table() %>% print()

  n_pb <- plot_data %>% nrow()
  n_pb_1.0 <- plot_data %>% filter(proportion_variants_phased == 1) %>% nrow()
  n_pb_0.75 <- plot_data %>% filter(proportion_variants_phased >= 0.75) %>% nrow()
  print(c("Percentage of phase blocks with all variants phased: ", 100*n_pb_1.0/n_pb))
  print(c("Percentage of phase blocks with half variants phased: ", 100*n_pb_0.75/n_pb))

  min_log2 <- plot_data %>% pull(somatic_mutations_per_Mb) %>%
    log2() %>% min() %>% round(digits = 0)
  max_log2 <- plot_data %>% pull(somatic_mutations_per_Mb) %>%
    log2() %>% max() %>% round(digits = 0)
  my_breaks <- seq(from = min_log2, to = max_log2, by = 1)

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

  rm(phasing_variants_grouped, plot_data, stat_data, p, q, q_with_legend,
     n_pb, n_pb_0.75, n_pb_1_sm, n_pb_1.0, pb_1_sm, sm_mb, tot_len, tot_sm,
     min_log2, max_log2, my_breaks)
}

# looking at important mutations
{
  important_ones <- important_mutations_tbl %>%
    rowwise() %>%
    mutate(p1 = str_sub(str_split(protein, pattern = "\\.", simplify = TRUE)[2], 1, 1),
           p2 = str_sub(str_split(protein, pattern = "\\.", simplify = TRUE)[2], -1),
           synonymous = (p1 == p2),
           Variant = str_c(chr, pos, ref, alt, sep = ":"))

  tumor_vafs <- important_mutations_vaf_tbl %>%
    filter(timepoint != "Normal") %>%
    mutate(variant_key = str_c(chr, pos, ref, alt, sep = ":"))
  normal_vafs <- important_mutations_vaf_tbl %>%
    filter(timepoint == "Normal")  %>%
    mutate(variant_key = str_c(chr, pos, ref, alt, sep = ":")) %>%
    mutate(normal_vaf = vaf,
           normal_alleles = alleles)

  plot_df <- sombx_all_mapq20_tbl %>%
    mutate(Variant = str_c(chromosome, position, ref, alt, sep = ":")) %>%
    left_join(important_ones, by = c("Variant", "ref", "alt")) %>%
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
           vaf > 0.05,
           nchar(alleles) > 10,
           n_alt_h1_bx + n_alt_h2_bx > 0) %>%
    left_join(phasing_variants_all_mapq20_tbl %>%
                select(-c("patient", "timepoint", "vcf_column", "sorted",
                          "cnv_maf_status", "display_name", "sample_n",
                          "my_color_100", "my_color_75", "my_color_50",
                          "my_color_25", "my_shape")),
              by = c("variant_key" = "Variant", "sample")) %>%
    mutate(variant_name = str_c(gene, protein, sep = " ")) %>%
    mutate(phased_by_SH = pct_ALT_on_H1 <= H1_proportion) %>%
    replace_na(list(Genotype = "Missing"))

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

  max_vaf <- plot_df %>% pull(barcode_vaf) %>% max() %>% plyr::round_any(.1, f = ceiling)
  ggplot(plot_df, aes(y = variant_name, x = barcode_vaf)) +
    geom_segment(aes(yend = variant_name, xend = 0), show.legend = FALSE) +
    geom_point(aes(shape = Genotype, color = phased_by_SH), show.legend = FALSE) +
    geom_text(aes(label = variant_name), size = 2, hjust = 0, nudge_x = 0.025) +
    facet_wrap(~display_name, ncol = 1, scales = "free_y", strip.position = "left") +
    labs(x = "Variant Allele Frequency (VAF)", y = NULL) +
    expand_limits(x = c(0, max_vaf)) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.text.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(main, "mm_mutation.pdf"),
           height = 4.75, width = 2.25, useDingbats = FALSE)

  rm(plot_df, correlation_value, max_vaf, H1_proportion)

}


rm(main, supp)
