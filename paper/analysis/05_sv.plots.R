################################################################################
# SVs in our data
################################################################################

library(tidyverse)

data_dir = file.path("data_for_plots/05_sv")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

main = "figures/05_sv/main/"
supp = "figures/05_sv/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# plot barcodes
if (TRUE) {

  sample_ids <- list.files(data_dir) %>%
    str_split(pattern = "\\.") %>%
    purrr::map_chr(function(x) x[2])

  sv_barcodes_tbl <- purrr::map(sample_ids,
                                function(x) read_tsv(file = file.path(data_dir, str_c("sv_barcodes_plot_df.", x, ".tsv")),
                                                     show_col_types = FALSE) %>%
                                  mutate(chrom = factor(chrom,
                                                        levels = str_c("chr", seq(1:22)),
                                                        ordered = TRUE)) %>%
                                  mutate(read_haplotype = factor(read_haplotype,
                                                                 levels = c(0,1,2),
                                                                 ordered = TRUE)) %>%
                                  mutate(supports_trans = factor(supports_trans,
                                                                 levels = c(0,1),
                                                                 ordered = TRUE)))

  names(sv_barcodes_tbl) <- sample_ids

  large_barcode_plot_translocation_only <- function(barcodes_tbl, sample_id,
                                                    breakpoint_pair,
                                                    display_name,
                                                    translocation_name,
                                                    filename,
                                                    my_height, my_width,
                                                    downsample = 1,
                                                    keep_bc = NULL,
                                                    min_max = FALSE,
                                                    min_chr = 2){

    breakpoints_list <- str_split(str_split(breakpoint_pair,
                                            pattern = ",",
                                            simplify = TRUE),
                                  pattern = ":") %>% unique()

    if (downsample < 1 & downsample > 0) {
      number_of_rows <- barcodes_tbl[[sample_id]] %>% nrow()
      bc_tbl <- barcodes_tbl[[sample_id]] %>%
        filter(as.logical(rbernoulli(number_of_rows, p = downsample))) %>%
        filter(breakpoints %in% breakpoint_pair)
    } else {
      bc_tbl <- barcodes_tbl[[sample_id]] %>%
        filter(breakpoints %in% breakpoint_pair) %>%
        filter(barcode %in% keep_bc)
    }

    if (is.null(keep_bc)) {
      consistent_barcodes <- bc_tbl %>%
        filter(read_haplotype %in% c(1, 2)) %>%
        group_by(chrom, barcode) %>%
        summarize(n_haps = length(unique(read_haplotype)),
                  .groups = "drop") %>%
        filter(n_haps == 1) %>%
        pull(barcode) %>%
        unique()
    } else {
      consistent_barcodes <- keep_bc
    }

    if (min_max) {
      min_max_read_positions <- bc_tbl %>%
        filter(supports_trans == 1,
               read_haplotype %in% c(1, 2),
               barcode %in% consistent_barcodes) %>%
        group_by(breakpoints, chrom, position) %>%
        summarize(min = min(read_position), max = max(read_position),
                  .groups = "drop")
    } else {
      min_max_read_positions <- bc_tbl %>%
        group_by(breakpoints, chrom, position) %>%
        summarize(min = 0, max = Inf,
                  .groups = "drop")
    }

    if (is.null(keep_bc)) {
      barcodes_both_chrom <- bc_tbl %>%
        left_join(min_max_read_positions,
                  by = c("breakpoints", "chrom", "position")) %>%
        filter(read_position >= min, read_position <= max) %>%
        filter(barcode %in% consistent_barcodes) %>%
        group_by(barcode) %>%
        summarize(both_chrom = length(unique(chrom)) == min_chr,
                  gemtools_bc = any(supports_trans == 1),
                  .groups = "drop") %>%
        filter(both_chrom == TRUE) %>% # | gemtools_bc == TRUE) %>%
        pull(barcode) %>% unique()
    } else {
      barcodes_both_chrom <- keep_bc
    }

    plot_df <- bc_tbl %>%
      filter(barcode %in% barcodes_both_chrom) %>%
      #filter(barcode %in% consistent_barcodes) %>%
      left_join(min_max_read_positions,
                by = c("breakpoints", "chrom", "position")) %>%
      filter(read_position >= min, read_position <= max) %>%
      mutate(both_chrom_name = str_c("Translocation ", translocation_name)) %>%
      arrange(chrom, desc(read_position))

    barcode_order <- plot_df %>% pull(barcode) %>% unique()

    plot_df <- plot_df %>%
      mutate(barcode = factor(barcode, levels = barcode_order, ordered = TRUE))

    p <- ggplot(plot_df, aes(x = read_position/1e6, y = barcode)) +
      geom_point(aes(color = read_haplotype), size = 0.5, shape = 16, alpha = 1, show.legend = FALSE) +
      scale_color_manual(values = c("#bdbdbd", "#ae8dc1", "#7fbf7b"), drop = FALSE)

    for (i in 1:length(breakpoints_list)) {
      p <- p + geom_vline(data = plot_df %>% filter(chrom == breakpoints_list[[i]][1]),
                          aes_string(xintercept = as.numeric(breakpoints_list[[i]][2])/1e6),
                          lty = 2)
    }

    p <- p +
      facet_grid(both_chrom_name ~ chrom, scales = "free") +
      labs(x = "Read Mapping Position (Mb)", y = "Barcode (one per row)") +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            panel.background = element_blank(),
            panel.grid = element_line(size = 0.5),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.placement = "right",
            strip.text = element_text(size = 8))

    ggsave(str_c(supp, filename, ".pdf"), p,
           width = my_width, height = my_height, useDingbats = FALSE)
  }

  set.seed(1)
  large_barcode_plot_translocation_only(sv_barcodes_tbl, "27522_1", "chr14:105858088,chr4:1871962", "27522 (P)", "t(4;14)", "27522_1.t414", 3.5, 5.25, downsample = 0.1, min_max = TRUE)

  keep_t4617_bc <- sv_barcodes_tbl[["27522_1"]] %>%
    filter(breakpoints %in% c("chr17:76201083,chr6:37194888", "chr4:152580898,chr6:37208806")) %>%
    select(breakpoints, chrom, barcode) %>%
    unique() %>%
    group_by(barcode) %>%
    summarize(n_chrom = length(unique(chrom))) %>%
    filter(n_chrom == 3) %>%
    pull(barcode) %>%
    unique()

  set.seed(1)
  large_barcode_plot_translocation_only(sv_barcodes_tbl, "27522_1", c("chr17:76201083,chr6:37194888", "chr4:152580898,chr6:37208806"), "27522 (P)", "t(4;6) + t(6;17)", "27522_1.t4617", 3, 7.25, keep_bc = keep_t4617_bc, min_chr = 3)

  rm(keep_t4617_bc)

  set.seed(1)
  large_barcode_plot_translocation_only(sv_barcodes_tbl, "27522_1", c("chr17:76201083,chr6:37194888"), "27522 (P)", "t(6;17)", "27522_1.t617", 3, 7.25, downsample = 0.1, min_max = TRUE)

  set.seed(1)
  large_barcode_plot_translocation_only(sv_barcodes_tbl, "27522_1", c("chr4:152580898,chr6:37208806"), "27522 (P)", "t(4;6)", "27522_1.t46", 3, 7.25, downsample = 0.1, min_max = TRUE)

  set.seed(1)
  large_barcode_plot_translocation_only(sv_barcodes_tbl, "77570", c("chr11:69404963,chr14:105741942", "chr11:69404963,chr14:106269142"), "77570 (P)", "t(11;14)", "77570.t1114", 3.5, 5.25, downsample = 0.1, min_max = TRUE)

}

# MMRF SV
if (TRUE) {

  mmrf_tbl <- tribble(~SV, ~target, ~Barwick, ~Barwick_label, ~Manier,
                      "t(11;14)", "CCND1", 16.0, NA, "15-20%",
                      "t(4;14)", "WHSC1", 11.0, "11.0%", "15%",
                      "t(14;16)", "MAF", 3.3, "3.3%", "5%",
                      "t(6;14)", "CCND3", 1.1, "1.1%", "1-2%",
                      "t(14;20)", "MAFB", 1.0, "1.0%", "~1%") %>%
    mutate(Barwick_axis = str_c(SV, target, sep = "\n"))

  max_label = "16.0%"
  max_pct = 16
  max_t = "t(11;14)\nCCND1"

  ggplot(data = mmrf_tbl,
         aes(x = Barwick, y = fct_reorder(Barwick_axis, Barwick))) +
    geom_col() +
    geom_text(aes(label = Barwick_label), hjust = 0, nudge_x = 0.1, size = 8/ggplot2:::.pt) +
    annotate("text", x = max_pct - 0.1, y = max_t, label = max_label, hjust = 1, size = 8/ggplot2:::.pt, color = "white") +
    annotate("text", x = Inf, y = -Inf, label = "Reported by\nBarwick, et al.", hjust = 1, vjust = 0, size = 8/ggplot2:::.pt) +
    labs(title = "Multiple Myeloma Translocations", x = "Percentage of MMRF with Translocation", y = NULL) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, mmrf_tbl %>% pull(Barwick) %>% max() + 1)) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 10),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines")) #+
    #ggsave(str_c(supp, "barwick_translocations.pdf"), width = 2, height = 2, useDingbats = FALSE)

  rm(mmrf_tbl, max_label, max_pct, max_t)
}

# Our data SVs
if (TRUE) {

  barcodes_tbl <- tribble(~sample, ~event, ~bp1, ~bp2, ~ps_id, ~H1, ~H2, ~max,
                          "27522 (P)", "t(4;14)", "chr4:1871962", "chr14:105858088", "chr14:105634857", 0, 17, 96,
                          "27522 (P)", "t(4;14)", "chr4:1871962", "chr14:105858088", "chr14:105858995", 0, 96, 96,
                          "27522 (P)", "t(4;14)", "chr4:1871962", "chr14:105858088", "chr4:228869", 1, 82, 96,
                          "27522 (Rel)", "t(4;14)", "chr4:1871962", "chr14:105858093", "chr14:105699569", 55, 0, 58,
                          "27522 (Rel)", "t(4;14)", "chr4:1871962", "chr14:105858093", "chr14:105858970", 0, 58, 58,
                          "27522 (Rel)", "t(4;14)", "chr4:1871962", "chr14:105858093", "chr4:222307", 0, 13, 58,
                          "77570 (Rel)", "t(11;14)", "chr11:69404963", "chr14:106269142", "(1) chr11:64789730", 12, 0, 37,
                          "77570 (Rel)", "t(11;14)", "chr11:69404963", "chr14:106269142", "(1) chr14:105596163", 12, 0, 37,
                          "77570 (Rel)", "t(11;14)", "chr11:69404963", "chr14:105741942", "(2) chr11:64789730", 37, 0, 37,
                          "77570 (Rel)", "t(11;14)", "chr11:69404963", "chr14:105741942", "(2) chr14:105596163", 22, 0,  37) %>%
    gather(H1, H2, key = "Haplotype", value = "Count") %>%
    mutate(strip_label = str_c(sample, event, sep = " "))

  barcodes_tbl %>%
    ggplot(aes(x = Count, y = ps_id, fill = fct_rev(Haplotype))) +
    geom_col(position = "dodge", width = 0.6, show.legend = FALSE) +
    geom_text(aes(x = Inf, label = ps_id), hjust = 1, size = 8/ggplot2:::.pt) +
    geom_text(aes(x = 0, label = Haplotype), hjust = 0, size = 8/ggplot2:::.pt, position = position_dodge(width = 0.6)) +
    facet_wrap(~strip_label, nrow = 1, scales = "free") +
    scale_fill_manual(values = c("#7fbf7b", "#ae8dc1")) +
    labs(x = "Number of SV-supporting Barcodes Assigned to Each Haplotype", y = "Phase Sets") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.text.y = element_blank(),
          #axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.title = element_text(size = 8),
          panel.background = element_blank(),
          #panel.border = element_blank(),
          panel.grid = element_line(size = 0.5),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines")) #+
    #ggsave(str_c(supp, "sample_translocations.pdf"), width = 5, height = 2, useDingbats = FALSE)

  rm(barcodes_tbl)

}

rm(main, supp, data_dir, sample_ids, sv_barcodes_tbl)
