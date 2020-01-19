################################################################################
# Closer look at allele relationships
################################################################################

main = "figures/04_alleles/main/"
supp = "figures/04_alleles/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# distance between mutations
{
  maximum_gene_length <- protein_coding_genes_tbl %>%
    mutate(gene_length = end - start) %>% pull(gene_length) %>% max()

  maximum_distance_10bx_mutations <- variant_pairs_mapq20_tbl %>%
    filter(weighted_log2.copyRatio > -.25,
           weighted_log2.copyRatio < 0.2) %>%
    filter(n_overlapping_barcodes >= 10) %>%
    pull(distance_between_variants) %>% max()

  variant_pairs_mapq20_tbl %>%
    filter(weighted_log2.copyRatio > -.25,
           weighted_log2.copyRatio < 0.2) %>%
    filter(distance_between_variants <= maximum_distance_10bx_mutations,
           distance_between_variants <= 10) %>%
    mutate(overlap_category = case_when(n_overlapping_barcodes == 0 ~ "Zero",
                                        n_overlapping_barcodes < 10 ~ "< 10",
                                        n_overlapping_barcodes < 20 ~ "< 20",
                                        n_overlapping_barcodes < 50 ~ "< 50",
                                        n_overlapping_barcodes < 100 ~ "< 100",
                                        TRUE ~ "100+")) %>%
    mutate(overlap_category = fct_rev(factor(overlap_category,
                                             levels = c("Zero", "< 10", "< 20",
                                                        "< 50", "< 100", "100+"),
                                             ordered = TRUE))) %>%
    ggplot(aes(x = distance_between_variants)) +
    geom_histogram(aes(fill = overlap_category), boundary = 1, binwidth = 1) +
    scale_fill_viridis(discrete = TRUE, option = "C") +
    theme_bw()
}

# coverage of two sites
{

}

# automated allele combination stats
{
  somatic_barcodes_mod <- somatic_barcodes %>% rowwise() %>%
    mutate(
      n_ref_barcodes_H1 = case_when(
        !is.na(ref_barcodes_H1) ~
          length(str_split(ref_barcodes_H1, ";", simplify = TRUE))),
      n_ref_barcodes_H2 = case_when(
        !is.na(ref_barcodes_H2) ~
          length(str_split(ref_barcodes_H2, ";", simplify = TRUE))),
      n_ref_barcodes_None = case_when(
        !is.na(ref_barcodes_None) ~
          length(str_split(ref_barcodes_None, ";", simplify = TRUE))),
      n_alt_barcodes_H1 = case_when(
        !is.na(alt_barcodes_H1) ~
          length(str_split(alt_barcodes_H1, ";", simplify = TRUE))),
      n_alt_barcodes_H2 = case_when(
        !is.na(alt_barcodes_H2) ~
          length(str_split(alt_barcodes_H2, ";", simplify = TRUE))),
      n_alt_barcodes_None = case_when(
        !is.na(alt_barcodes_None) ~
          length(str_split(alt_barcodes_None, ";", simplify = TRUE))),
      n_alt_barcodes = sum(c(n_alt_barcodes_H1, n_alt_barcodes_H2, n_alt_barcodes_None), na.rm = TRUE),
      n_total_barcodes = sum(c(n_alt_barcodes_H1, n_alt_barcodes_H2, n_alt_barcodes_None, n_ref_barcodes_H1, n_ref_barcodes_H2, n_ref_barcodes_None), na.rm = TRUE),
      barcode_vaf = case_when(n_alt_barcodes == 0 ~ 0,
                              n_total_barcodes == 0 ~ 0,
                              TRUE ~ n_alt_barcodes/n_total_barcodes)) %>% ungroup()

  ################################################################################

  somatic_barcodes_n <- somatic_barcodes_mod %>%
    select(Sample, variant_key, n_alt_barcodes, n_total_barcodes, barcode_vaf)

  variant_pairs_coverage <- variant_pairs %>%
    left_join(somatic_barcodes_n, by = c("Sample", "Variant1" = "variant_key")) %>%
    left_join(somatic_barcodes_n, by = c("Sample", "Variant2" = "variant_key"))

  n_variant_pairs <- variant_pairs_coverage %>% nrow()
  n_variant_pairs_both_covered <- variant_pairs_coverage %>%
    filter(n_alt_barcodes.x > 3, n_alt_barcodes.y > 3,
           n_total_barcodes.x > 10, n_total_barcodes.y > 10) %>% nrow()
  n_variant_pairs_share_1_barcode <- variant_pairs_coverage %>%
    filter(n_alt_barcodes.x > 3, n_alt_barcodes.y > 3,
           n_total_barcodes.x > 10, n_total_barcodes.y > 10) %>%
    filter(n_bx_overlap_00 != 0 | n_bx_overlap_01 != 0 | n_bx_overlap_10 != 0 | n_bx_overlap_11 != 0) %>% nrow()
  n_variant_pairs_share_1_barcode_variant <- variant_pairs_coverage %>%
    filter(n_alt_barcodes.x > 3, n_alt_barcodes.y > 3,
           n_total_barcodes.x > 10, n_total_barcodes.y > 10) %>%
    filter(n_bx_overlap_00 != 0 | n_bx_overlap_01 != 0 | n_bx_overlap_10 != 0 | n_bx_overlap_11 != 0) %>%
    filter(n_bx_overlap_01 != 0 | n_bx_overlap_10 != 0 | n_bx_overlap_11 != 0) %>% nrow()

  print(c("Number of variant pairs: ", n_variant_pairs))
  print(c("Number of variant pairs, both with coverage: ", n_variant_pairs_both_covered))
  print(c("Number of variant pairs, both with coverage, share >= 1 barcode: ", n_variant_pairs_share_1_barcode))
  print(c("Number of variant pairs, both with covearge, shard >= 1 barcode at variant allele: ", n_variant_pairs_share_1_barcode_variant))

  variant_pairs_coverage %>%
    filter(n_alt_barcodes.x > 3, n_alt_barcodes.y > 3,
           n_total_barcodes.x > 10, n_total_barcodes.y > 10) %>%
    filter(n_bx_overlap_00 != 0 | n_bx_overlap_01 != 0 | n_bx_overlap_10 != 0 | n_bx_overlap_11 != 0) %>%
    filter(n_bx_overlap_01 != 0 | n_bx_overlap_10 != 0 | n_bx_overlap_11 != 0) %>%
    rowwise() %>%
    mutate(combo = str_c(as.numeric(n_bx_overlap_00 != 0),
                         as.numeric(n_bx_overlap_01 != 0),
                         as.numeric(n_bx_overlap_10 != 0),
                         as.numeric(n_bx_overlap_11 != 0))) %>%
    ungroup() %>% pull(combo) %>% table() %>% sort()
}

# allele combinations (handmade with known figures)
{
  plot_df <- tribble(~level1, ~level2, ~value,
                     "Both\ncovered",  "Yes",  6,
                     "Both\ncovered",   "No", 94,
                     "Share\nbarcode", "Yes", 22,
                     "Share\nbarcode",  "No", 78,
                     "Share\nvariant", "Yes", 67,
                     "Share\nvariant",  "No", 33)

  p <- ggplot(data = plot_df, aes(x = level1, y = value, fill = level2)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(data = plot_df %>% filter(level2 == "Yes"), aes(label = value),
              vjust = 1, nudge_y = -1,
              color = "white",
              size = 2.5) +
    geom_text(data = plot_df %>% filter(level2 == "No"), aes(label = value, y = 100),
              vjust = 1, nudge_y = -1,
              color = "white",
              size = 2.5) +
    geom_segment(x = 1.25, y = 6, xend = 1.75, yend = 100, lty = 2, color = "#bdbdbd", lwd = 0.3) +
    geom_segment(x = 1.25, y = 0, xend = 1.75, yend = 0, lty = 2, color = "#bdbdbd", lwd = 0.3) +
    geom_segment(x = 2.25, y = 22, xend = 2.75, yend = 100, lty = 2, color = "#bdbdbd", lwd = 0.3) +
    geom_segment(x = 2.25, y = 0, xend = 2.75, yend = 0, lty = 2, color = "#bdbdbd", lwd = 0.3) +
    geom_segment(x = 3.25, y = 67, xend = 3.75, yend = 100, lty = 2, color = "#bdbdbd", lwd = 0.3) +
    geom_segment(x = 3.25, y = 0, xend = 3.75, yend = 0, lty = 2, color = "#bdbdbd", lwd = 0.3) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_manual(values = c("#80cdc1", "#01665e")) +
    scale_x_discrete(expand = c(0,1)) +
    labs(x = NULL, y = "Somatic Mutation Pairs (%)", fill = NULL)

  q <- p + guides(fill = FALSE)

  ggsave("n_variant_pairs.pdf", p, width = 5, height = 5, useDingbats = FALSE)
  ggsave("n_variant_pairs.no_legend.pdf", q, width = 5, height = 5, useDingbats = FALSE)

  ################################################################################

  text_df <- tribble(~category, ~total,
                     7,    969,
                     6,    274,
                     5,    252,
                     4,    210,
                     3,     97,
                     2,     90,
                     1,     96)

  plot_df <- tribble(~alleles, ~category, ~filled,
                     1, 7, 1,
                     2, 7, 0,
                     3, 7, 0,
                     4, 7, 1,
                     1, 6, 1,
                     2, 6, 0,
                     3, 6, 1,
                     4, 6, 0,
                     1, 5, 1,
                     2, 5, 1,
                     3, 5, 0,
                     4, 5, 0,
                     1, 4, 1,
                     2, 4, 1,
                     3, 4, 1,
                     4, 4, 0,
                     1, 3, 1,
                     2, 3, 0,
                     3, 3, 1,
                     4, 3, 1,
                     1, 2, 1,
                     2, 2, 1,
                     3, 2, 0,
                     4, 2, 1,
                     1, 1, 2,
                     2, 1, 2,
                     3, 1, 2,
                     4, 1, 2)

  p <- ggplot() +
    geom_tile(data = plot_df,
              aes(x = alleles, y = category, fill = as.factor(filled)),
              color = NA) +
    geom_text(data = text_df, aes(y = category, x = 5, label = total)) +
    geom_vline(xintercept = seq(1,3) + 0.5, lty = 1, lwd = 0.5, color = "#ffffff") +
    geom_hline(yintercept = seq(1,6) + 0.5, lty = 1, lwd = 2, color = "#ffffff") +
    coord_equal() +
    scale_fill_manual(values = c("#80cdc1", "#01665e", "#bdbdbd")) +
    scale_x_continuous(limits = c(0.5,5.5),
                       breaks = seq(1,5),
                       labels = c("REF/REF", "REF/ALT", "ALT/REF", "ALT/ALT", "Total"),
                       position = "top") +
    scale_y_continuous(expand = c(0,0), position = "right") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank()) +
    labs(y = "Observed Combinations of Linked Somatic Mutations", x = NULL) +
    guides(fill = FALSE)

  ggsave("allele_combinations.pdf", p, height = 5, width = 5, useDingbats = FALSE)

}

# simplified variant pairs
{
  #bcv <- read.delim("27522_1.NRAS.barcodes_variants.tsv")
  bcv <- read.delim("~/Desktop/27522_1.no_NRAS_Q61.barcodes_variants.tsv")
  variants_of_interest <- c("chr1:114716124:C:G", "chr1:114713909:G:T")

  bcv_filtered <- bcv %>%
    filter(Allele != "No Coverage" &
             ((Filter %in% c("PASS", "10X_PHASING_INCONSISTENT") &
                 Genotype != "1|1") |
                Variant %in% variants_of_interest))

  variants_more_than_1 <- bcv_filtered %>%
    group_by(Variant) %>%
    summarize(count = n()) %>%
    filter(count > 1) %>%
    mutate(Variant_rank = seq(1:n()))

  barcodes_more_than_1 <- bcv_filtered %>%
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
    arrange(desc(proportion_H1), proportion_H2) %>%
    select(Barcode) %>%
    mutate(Barcode_rank = seq(1:n())) %>%
    mutate(Barcode_ordered = factor(Barcode_rank, labels = Barcode))

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

  plot_data %>% left_join(horizontal_lines, by = "Barcode") %>%
    ggplot(aes(x = Variant_rank, y = Barcode_rank,
               label = Allele, fill = Haplotype)) +
    geom_segment(aes(y = y, yend = yend, x = x, xend = xend)) +
    geom_vline(xintercept = seq(1:n_variants), linetype = 1, alpha = 0.25) +
    geom_hline(yintercept = seq(1:n_barcodes), linetype = 1, alpha = 0.25) +
    geom_tile() +
    geom_text() +
    coord_equal(ratio = 1) +
    labs(x = "Variant", y = "Barcode", fill = "Haplotype") +
    theme_bw(base_size = 20) +
    scale_x_continuous(breaks = seq(1:n_variants),
                       labels = plot_data %>% pull(Variant) %>% sort() %>% unique()) +
    scale_y_continuous(breaks = seq(1:n_barcodes),
                       labels = plot_data %>% arrange(Barcode_rank) %>% pull(Barcode) %>% unique()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(panel.grid.major.y = element_blank()) +
    theme(panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.x = element_blank())

  ggsave("~/Desktop/27522_1.no_NRAS_Q61.barcode_variants.pdf", width = 30, height = 20, useDingbats = FALSE)

}

# other one
{
  library(tidyverse)

  bcv <- read.delim("27522_1.NRAS.barcodes_variants.tsv")

  variants_of_interest <- c("chr1:114713909:G:T", "chr1:114716124:C:G")

  bcv_filtered <- bcv %>%
    filter(Allele != "No Coverage" &
             ((Filter %in% c("PASS", "10X_PHASING_INCONSISTENT") &
                 Genotype != "1|1") |
                Variant %in% variants_of_interest))

  variants_more_than_1 <- bcv_filtered %>%
    group_by(Variant) %>%
    summarize(count = n()) %>%
    filter(count > 1) %>%
    mutate(Variant_rank = seq(1:n())) %>%
    select(Variant, Variant_rank)

  n_variants <- variants_more_than_1 %>% nrow()

  barcodes_more_than_1 <- bcv_filtered %>%
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
    arrange(desc(proportion_H1), proportion_H2) %>%
    select(Barcode) %>%
    mutate(Barcode_rank = seq(1:n())) %>%
    mutate(Barcode_ordered = factor(Barcode_rank, labels = Barcode))

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

  plot_data <- plot_data %>% left_join(horizontal_lines, by = "Barcode")

  ################################################################################

  barcodes_support_variant1_0 <- plot_data %>%
    filter(Variant == variants_of_interest[1], Allele == 0) %>% pull(Barcode)

  barcodes_support_variant1_1 <- plot_data %>%
    filter(Variant == variants_of_interest[1], Allele == 1) %>% pull(Barcode)

  barcodes_support_variant2_0 <- plot_data %>%
    filter(Variant == variants_of_interest[2], Allele == 0) %>% pull(Barcode)

  barcodes_support_variant2_1 <- plot_data %>%
    filter(Variant == variants_of_interest[2], Allele == 1) %>% pull(Barcode)

  barcodes_support_both_00 <- intersect(barcodes_support_variant1_0,
                                        barcodes_support_variant2_0)

  barcodes_support_both_01 <- intersect(barcodes_support_variant1_0,
                                        barcodes_support_variant2_1)

  barcodes_support_both_10 <- intersect(barcodes_support_variant1_1,
                                        barcodes_support_variant2_0)

  barcodes_support_both_11 <- intersect(barcodes_support_variant1_1,
                                        barcodes_support_variant2_1)

  barcodes_support_one_0x <- barcodes_support_variant1_0[
    !(barcodes_support_variant1_0 %in% barcodes_support_variant2_0) &
      !(barcodes_support_variant1_0 %in% barcodes_support_variant2_1)]

  barcodes_support_one_1x <- barcodes_support_variant1_1[
    !(barcodes_support_variant1_1 %in% barcodes_support_variant2_0) &
      !(barcodes_support_variant1_1 %in% barcodes_support_variant2_1)]

  barcodes_support_one_x0 <- barcodes_support_variant2_0[
    !(barcodes_support_variant2_0 %in% barcodes_support_variant1_0) &
      !(barcodes_support_variant2_0 %in% barcodes_support_variant1_1)]

  barcodes_support_one_x1 <- barcodes_support_variant2_1[
    !(barcodes_support_variant2_1 %in% barcodes_support_variant1_0) &
      !(barcodes_support_variant2_1 %in% barcodes_support_variant1_1)]

  max_hap_vector <- NULL
  for (bx in barcodes_more_than_1) {
    haps <- plot_data %>% filter(Barcode == bx) %>% pull(Haplotype) %>% unique()
    if ("H1" %in% haps & !("H2" %in% haps)) {
      max_hap <- "H1"
    } else if ("H2" %in% haps & !("H1" %in% haps)) {
      max_hap <- "H2"
    } else if ("Not Phased Heterozygote" %in% haps & length(haps) == 1) {
      max_hap <- "Not Phased Heterozygote"
    } else {
      max_hap <- plot_data %>% filter(Barcode == bx) %>%
        pull(Haplotype) %>% table() %>% which.max() %>% names()
    }
    max_hap_vector <- c(max_hap_vector, max_hap)
  }
  max_hap_tbl <- as.tibble(data.frame(Barcode = barcodes_more_than_1,
                                      max_Haplotype = max_hap_vector))

  plot_data <- plot_data %>%
    mutate(mutation_pair = case_when(
      Barcode %in% barcodes_support_both_00 ~ "REF/REF",
      Barcode %in% barcodes_support_both_01 ~ "REF/ALT",
      Barcode %in% barcodes_support_both_10 ~ "ALT/REF",
      Barcode %in% barcodes_support_both_11 ~ "ALT/ALT",
      Barcode %in% barcodes_support_one_0x ~ "REF/NC",
      Barcode %in% barcodes_support_one_1x ~ "ALT/NC",
      Barcode %in% barcodes_support_one_x0 ~ "NC/REF",
      Barcode %in% barcodes_support_one_x1 ~ "NC/ALT")) %>%
    left_join(max_hap_tbl, by = "Barcode")

  mutation_pair_rank <- tribble(~mutation_pair, ~mutation_pair_rank,
                                "NC/ALT",  1,
                                "ALT/NC",  2,
                                "NC/REF",  3,
                                "REF/NC",  4,
                                "ALT/ALT", 5,
                                "ALT/REF", 6,
                                "REF/ALT", 7,
                                "REF/REF", 8)

  plot_data <- plot_data %>%
    left_join(mutation_pair_rank, by = "mutation_pair") %>%
    group_by(Variant, Variant_rank, Allele, max_Haplotype, mutation_pair, mutation_pair_rank) %>%
    summarize(count = n(), x = min(x), xend = max(xend)) %>%
    mutate(hap_allele = str_c(max_Haplotype, Allele)) %>%
    mutate(Haplotype_Allele = case_when(
      hap_allele == "H10" ~ "Haplotype 1, REF",
      hap_allele == "H11" ~ "Haplotype 1, ALT",
      hap_allele == "H20" ~ "Haplotype 2, REF",
      hap_allele == "H21" ~ "Haplotype 2, ALT",
      hap_allele == "Not Phased Heterozygote0" ~ "Not Phased, REF",
      hap_allele == "Not Phased Heterozygote1" ~ "Not Phased, ALT")) %>%
    mutate(Haplotype_Allele_factor = factor(Haplotype_Allele, levels = c("Haplotype 1, REF", "Haplotype 1, ALT", "Haplotype 2, REF", "Haplotype 2, ALT", "Not Phased, REF", "Not Phased, ALT"), ordered = TRUE)) %>%
    mutate(max_Haplotype_full = case_when(
      max_Haplotype == "H1" ~ "Haplotype 1",
      max_Haplotype == "H2" ~ "Haplotype 2",
      max_Haplotype == "Not Phased Heterozygote" ~ "Haplotype 2")) %>%
    mutate(max_Haplotype_full_2 = case_when(
      max_Haplotype == "H1" ~ "Haplotype 1",
      max_Haplotype == "H2" ~ "Haplotype 2",
      max_Haplotype == "Not Phased Heterozygote" ~ "Haplotype 11")) %>%
    ungroup()

  variant_labels <- variants_more_than_1 %>% mutate(variant_labels = case_when(
    as.character(Variant) == "chr1:114716124:C:G" ~ "G13",
    as.character(Variant) == "chr1:114713909:G:T" ~ "Q61",
    TRUE ~ "" )) %>% pull(variant_labels)

  plot_data %>%
    ggplot(aes(x = Variant_rank, y = mutation_pair_rank, label = count, fill = Haplotype_Allele_factor)) +
    geom_segment(aes(y = mutation_pair_rank, yend = mutation_pair_rank, x = x, xend = xend)) +
    geom_vline(xintercept = seq(1:n_variants), linetype = 1, alpha = 0.25) +
    geom_hline(yintercept = seq(1:8), linetype = 1, alpha = 0.25) +
    geom_vline(xintercept = plot_data %>% filter(Variant %in% variants_of_interest) %>% pull(Variant_rank) %>% unique(), size = 1.5) +
    geom_tile() +
    geom_text() +
    coord_equal(ratio = 1) +
    labs(x = NULL, y = "NRAS Q61/G13 Mutation Combination", fill = "Haplotype, Allele") +
    theme_bw() +
    scale_x_continuous(breaks = seq(1:n_variants),
                       labels = c(variant_labels),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(1:8),
                       labels = mutation_pair_rank$mutation_pair,
                       expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5)) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_blank()) +
    facet_wrap(~ max_Haplotype_full, ncol = 1, strip.position = c("top")) +
    scale_fill_brewer(palette = "Paired")

  ggsave("27522_1.NRAS.condensed_haplotypes.pdf", width = 20, height = 10, useDingbats = FALSE)

  ########

  plot_data %>%
    ggplot(aes(x = Variant_rank, y = mutation_pair_rank, label = count, fill = Haplotype_Allele_factor)) +
    geom_segment(aes(y = mutation_pair_rank, yend = mutation_pair_rank, x = x, xend = xend)) +
    geom_vline(xintercept = seq(1:n_variants), linetype = 1, alpha = 0.25) +
    geom_hline(yintercept = seq(1:8), linetype = 1, alpha = 0.25) +
    geom_vline(xintercept = plot_data %>% filter(Variant %in% variants_of_interest) %>% pull(Variant_rank) %>% unique(), size = 1.5) +
    geom_tile() +
    geom_text() +
    coord_equal(ratio = 1) +
    labs(x = NULL, y = "NRAS Q61/G13 Mutation Combination", fill = "Haplotype, Allele") +
    theme_bw() +
    scale_x_continuous(breaks = seq(1:n_variants),
                       labels = c(variant_labels),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(1:8),
                       labels = mutation_pair_rank$mutation_pair,
                       expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5)) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_blank()) +
    facet_wrap(~ max_Haplotype_full_2, ncol = 1, strip.position = c("top")) +
    scale_fill_brewer(palette = "Paired")

  ggsave("27522_1.NRAS.condensed_haplotype.three_hap.pdf", width = 20, height = 10, useDingbats = FALSE)

}

# 27522_1 NRAS G13 Q61 figure
{
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
