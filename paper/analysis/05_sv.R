################################################################################
# SVs in our data
################################################################################

main = "figures/05_sv/main/"
supp = "figures/05_sv/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

# plot barcodes
if (TRUE) {

  large_barcode_plot <- function(barcodes_tbl, sample_id, breakpoint_pair, display_name, translocation_name, filename, my_height, my_width, gemtools_support_only = TRUE, min_max = TRUE){

    breakpoints_list <- str_split(str_split(breakpoint_pair,
                                            pattern = ",",
                                            simplify = TRUE),
                                  pattern = ":") %>% unique()

    bc_tbl <- barcodes_tbl[[sample_id]] %>%
      filter(breakpoints %in% breakpoint_pair)

    consistent_barcodes <- bc_tbl %>%
      filter(read_haplotype %in% c(1, 2)) %>%
      group_by(chrom, barcode) %>%
      summarize(n_haps = length(unique(read_haplotype))) %>%
      filter(n_haps == 1) %>%
      pull(barcode) %>%
      unique()

    if (min_max) {
      min_max_read_positions <- bc_tbl %>%
        filter(supports_trans == 1,
               read_haplotype %in% c(1, 2),
               barcode %in% consistent_barcodes) %>%
        group_by(breakpoints, chrom, position) %>%
        summarize(min = min(read_position), max = max(read_position))
    } else {
      min_max_read_positions <- bc_tbl %>%
        group_by(breakpoints, chrom, position) %>%
        mutate(min = 0, max = Inf)
    }

    if (gemtools_support_only) {
      barcodes_both_chrom <- bc_tbl %>%
        filter(supports_trans == 1) %>%
        pull(barcode) %>% unique()
    } else {
      barcodes_both_chrom <- bc_tbl %>%
        group_by(barcode) %>%
        summarize(both_chrom = length(unique(chrom)) > 1) %>%
        filter(both_chrom == TRUE) %>%
        pull(barcode) %>% unique()
    }

    barcodes_dominant_haplotype <- bc_tbl %>%
      filter(read_haplotype != 0) %>%
      group_by(barcode) %>%
      summarize(dominant_haplotype = median(as.numeric(read_haplotype) - 1))
      #summarize(dominant_haplotype = median(as.numeric(read_haplotype)))

    barcodes_both_sides <- bc_tbl %>%
      filter(!(barcode %in% barcodes_both_chrom)) %>%
      group_by(barcode, chrom) %>%
      summarize(both_sides = any(read_position <= position) & any(read_position >= position)) %>%
      filter(both_sides == TRUE) %>%
      pull(barcode) %>% unique()

    barcodes_trans_or_both_sides <- bc_tbl %>%
      filter(barcode %in% barcodes_both_chrom | barcode %in% barcodes_both_sides) %>%
      filter(barcode %in% consistent_barcodes) %>%
      pull(barcode) %>% unique()

    plot_df <- bc_tbl %>%
      filter(barcode %in% barcodes_trans_or_both_sides) %>%
      filter(barcode %in% consistent_barcodes) %>%
      left_join(min_max_read_positions,
                by = c("breakpoints", "chrom", "position")) %>%
      left_join(barcodes_dominant_haplotype,
                by = "barcode") %>%
      mutate(both_sides = barcode %in% barcodes_both_sides) %>%
      mutate(both_chrom = barcode %in% barcodes_both_chrom) %>%
      mutate(both_chrom_name = case_when(both_chrom ~ str_c("Translocation ", translocation_name),
                                         TRUE ~ "No Translocation")) %>%
      filter(read_position >= min & read_position <= max) %>%
      mutate(relative_to_breakpoint = read_position - position) %>%
      arrange(dominant_haplotype, read_position) #%>%
      #mutate(barcode = factor(barcode, barcode))

    p <- ggplot(plot_df, aes(x = read_position/1e6, y = barcode)) +
      geom_segment(aes(xend = (read_position + 1000)/1e6,
                       yend = barcode,
                       color = read_haplotype),
                 show.legend = FALSE) +
      scale_color_manual(values = c("#bdbdbd", "#ae8dc1", "#7fbf7b"))

    for (i in 1:length(breakpoints_list)) {
      p <- p + geom_vline(data = plot_df %>% filter(chrom == breakpoints_list[[i]][1]),
                          aes_string(xintercept = as.numeric(breakpoints_list[[i]][2])/1e6),
                          lty = 2)
    }

    p <- p +
      facet_grid(both_chrom_name ~ chrom, scales = "free") +
      labs(x = "Read Mapping Position (Mb)", y = "Barcode") +
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

    ggsave(str_c(main, filename, ".pdf"), p,
           width = my_width, height = my_height, useDingbats = FALSE)
  }

  large_barcode_plot(sv_barcodes_tbl, "27522_1", "chr14:105858088,chr4:1871962", "27522 (P)", "t(4;14)", "27522_1.t414.new", 3.25, 5)

  #large_barcode_plot(sv_barcodes_tbl, "27522_1", c("chr17:76201083,chr6:37194888", "chr4:152580898,chr6:37208806"), "27522 (P)", "t(4;6) + t(6;17)", "27522_1.t4617", 3.25, 7.25, gemtools_support_only = FALSE, min_max = FALSE )

  large_barcode_plot(sv_barcodes_tbl, "77570", c("chr11:69404963,chr14:105741942", "chr11:69404963,chr14:106269142"), "77570 (P)", "t(11;14)", "77570.t1114", 7, 10)


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
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "barwick_translocations.pdf"), width = 2, height = 2, useDingbats = FALSE)

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
          plot.margin = unit(c(0,0,0,0), "lines")) +
    ggsave(str_c(supp, "sample_translocations.pdf"), width = 5, height = 2, useDingbats = FALSE)

  rm(barcodes_tbl)

}
