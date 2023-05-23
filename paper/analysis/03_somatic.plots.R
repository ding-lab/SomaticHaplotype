################################################################################
# Somatic phasing
# Concordance with longranger
# Mutations on phase sets
# Mutation pairs per phase set
################################################################################

library(tidyverse)
library(ggrepel)

data_dir = file.path("data_for_plots/03_somatic")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

main = "figures/03_somatic_phasing/main/"
supp = "figures/03_somatic_phasing/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# Precision/recall balance threshold to maximize longranger concordance
{
  precision_recall_plot_df <- read_tsv(file = file.path(data_dir, "precision_recall_plot_df.tsv"),
                                       show_col_types = FALSE)

  ggplot(data = precision_recall_plot_df,
         aes(x = prec, y = rec)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "precision_recall.pdf"),
         width = 3.5, height = 3.5, useDingbats = FALSE)

  rm(precision_recall_plot_df)
}

# Concordance with longranger
{

  longranger_concordance_plot_df <- read_tsv(file = file.path(data_dir, "longranger_concordance_plot_df.tsv"),
                                             show_col_types = FALSE)

  ggplot(data = longranger_concordance_plot_df,
         aes(x = Genotype, y = fct_rev(phasing_pair), fill = n)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "longranger_concordance.pdf"),
         width = 3.5, height = 7.25, useDingbats = FALSE)

  rm(longranger_concordance_plot_df)

  n_nc <- read_tsv(file = file.path(data_dir, "n_nc.tsv"),
                   show_col_types = FALSE)$n_nc

  methods_phasing_plot_df <- read_tsv(file = file.path(data_dir, "methods_phasing_plot_df.tsv"),
                                      show_col_types = FALSE)

  ggplot(data = methods_phasing_plot_df,
         aes(x = phased_by_linked_alleles,
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(main, "methods_phasing.pdf"),
         width = 2.25, height = 2.25, useDingbats = FALSE)

  rm(n_nc, methods_phasing_plot_df)
}

# Phased somatic mutations on haplotypes
{
  sample_id <- "59114_4"
  chromosome_of_interest <- "chr1"
  min_length = 1e3

  chromosome_length_tbl <- read_tsv(file = file.path("data_for_plots/02_phase_sets/plot_chr.tsv"),
                                    show_col_types = FALSE)

  somatic_phase_sets_plot_df <- read_tsv(file = file.path(data_dir, "somatic_phase_sets_plot_df.tsv"),
                                         show_col_types = FALSE)

  somatic_mutation_plot_df <- read_tsv(file = file.path(data_dir, "somatic_mutation_plot_df.tsv"),
                                       show_col_types = FALSE)

  ggplot(somatic_phase_sets_plot_df) +
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
    geom_point(data = somatic_mutation_plot_df,
               aes(x = Position/1e6,
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(main, "somatic_mutations_on_phase_sets.pdf"),
           width = 4.75, height = 1.125, useDingbats = FALSE)

  rm(somatic_phase_sets_plot_df,
     somatic_mutation_plot_df,
     min_length, sample_id, chromosome_of_interest)

}

# variants per phase set
{

  somatic_mutations_per_Mb_plot_df <- read_tsv(file = file.path(data_dir, "somatic_mutations_per_Mb_plot_df.tsv"),
                                               show_col_types = FALSE)

  min_log2 <- somatic_mutations_per_Mb_plot_df %>% pull(somatic_mutations_per_Mb) %>%
    log2() %>% min() %>% round(digits = 0)
  max_log2 <- somatic_mutations_per_Mb_plot_df %>% pull(somatic_mutations_per_Mb) %>%
    log2() %>% max() %>% round(digits = 0)
  my_breaks <- seq(from = min_log2, to = max_log2, by = 1)

  p <- ggplot(data = somatic_mutations_per_Mb_plot_df,
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
  q <- ggExtra::ggMarginal(p + guides(color = FALSE, size = FALSE, scale = "none"),
                           type = "histogram",
                           groupFill = TRUE,
                           color = NA)

  ggsave(str_c(main, "somatic_mutations_per_Mb.with_legend.pdf"), q_with_legend,
         width = 4.75, height = 1.5*2.25, useDingbats = FALSE)

  ggsave(str_c(main, "somatic_mutations_per_Mb.pdf"), q,
         width = 4.75, height = 1.5*2.25, useDingbats = FALSE)

  rm(somatic_mutations_per_Mb_plot_df,
     p, q, q_with_legend,
     min_log2, max_log2, my_breaks)
}

# looking at driver mutations
{
  vaf_correlation_plot_df <- read_tsv(file = file.path(data_dir, "vaf_correlation_plot_df.tsv"),
                                      show_col_types = FALSE)

  correlation_value <- round(as.numeric(cor.test(vaf_correlation_plot_df$vaf,
                                                 vaf_correlation_plot_df$barcode_vaf)$estimate), 2)

  ggplot(vaf_correlation_plot_df, aes(x = vaf, y = barcode_vaf)) +
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
          plot.margin = unit(c(0,0,0,0), "lines"))

  ggsave(str_c(supp, "vaf_correlation.pdf"),
         height = 3.5, width = 3.5, useDingbats = FALSE)

  p <- ggplot(vaf_correlation_plot_df,
              aes(y = variant_name, x = display_name)) +
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

  rm(vaf_correlation_plot_df, correlation_value, p)

}

# Can we correct for copy number using linked alleles
if (FALSE) {

  for (this_sample in cnv_tbl %>% pull(sample) %>% unique()) {
    print(this_sample)
    # Somatic variants phased with CNV data
    x <- phasing_variants_driver_mapq20_tbl %>%
      filter(phased_by_linked_alleles %in% c("H1", "H2"),
             cnv_maf_status == TRUE) %>%
      select(Variant, Chromosome, Position, Reference, Alternate,
             Phase_Set, sample, phased_by_linked_alleles) %>%
      filter(sample == this_sample)

    # Multiplication factors
    y <- barcodes_variants_driver_mapq20_tbl[[this_sample]] %>%
      filter(Phased_Heterozygote == "True",
             Variant != Somatic_Variant) %>%
      group_by(Somatic_Variant, Variant) %>%
      summarize(n_H1 = sum(Haplotype == "H1"),
                n_H2 = sum(Haplotype == "H2"),
                n_total = n_H1 + n_H2,
                mult_H1 = 0.5*n_total/n_H1,
                mult_H2 = 0.5*n_total/n_H2) %>%
      filter(n_H1 >= 5, n_H2 >= 5) %>%
      ungroup() %>% group_by(Somatic_Variant) %>%
      summarize(total_variants = n(),
                mult_H1_mean = mean(mult_H1), mult_H1_sd = sd(mult_H1),
                mult_H2_mean = mean(mult_H2), mult_H2_sd = sd(mult_H2)) %>%
      filter(total_variants >= 10)

    get_cnv_for_position <- function(cnv, my_chr, my_pos, my_sample){
      return_value <- cnv %>%
        filter(sample == my_sample,
               chrom == my_chr,
               start <= my_pos,
               end >= my_pos) %>%
        pull(log2.copyRatio)

      if (length(return_value) == 0) {
        return_value <- NA
      }
      return(return_value)

    }

    x %>%
      left_join(y, by = c("Variant" = "Somatic_Variant")) %>%
      filter(!is.na(total_variants)) %>%
      rowwise() %>%
      mutate(log2.copyRatio = get_cnv_for_position(cnv_tbl,
                                                   my_chr = Chromosome,
                                                   my_pos = Position,
                                                   my_sample = sample)) %>%
      ungroup() %>%
      View()

    rm(x, y)

  }

}

rm(main, supp, data_dir, chromosome_length_tbl)
