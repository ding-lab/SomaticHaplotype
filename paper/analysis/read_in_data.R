################################################################################
# Read in all data for analysis
################################################################################

phase_proportion <- 0.91

# shh! (this is a library)
library(tidyverse)
library(viridis)
library(ggrepel)
library(fishplot)

# Load previous input data?
# Processed input data tables are saved in a date-stamped .Rdata object. This is
# a resusable input data object that is re-loaded at the beginning of each session.
# (We do not save the environment to .Rdata at the end of the session.)

# numbers used in manuscript
manuscript_numbers <- list()

last_updated <- "2020-05-06"
input_data_path_str <- str_c("data/collected_input_objects.", last_updated, ".RData")
if (file.exists(input_data_path_str)) {
  load(input_data_path_str)
  print(str_c("Data tables loaded from .RData file created ", last_updated, "."))
  rm(last_updated, input_data_path_str)

} else {

  # patient and sample names

  patient_sample_names_tbl <- read_tsv("data/patient_sample_names.tsv",
                                       col_types = "ccccll") %>%
    mutate(timepoint = factor(timepoint,
                              levels = c("SMM",
                                         "Primary",
                                         "Pre-transplant",
                                         "Post-transplant",
                                         "Remission",
                                         "Relapse-1",
                                         "Normal"),
                              labels = c("SMM",
                                         "Primary",
                                         "Pre-transplant",
                                         "Post-transplant",
                                         "Remission",
                                         "Relapse",
                                         "Normal"),
                              ordered = TRUE)) %>%
    mutate(display_name = case_when(timepoint == "SMM" ~ str_c(patient, " (S)"),
                                    timepoint == "Primary" ~ str_c(patient, " (P)"),
                                    timepoint == "Pre-transplant" ~ str_c(patient, " (PrT)"),
                                    timepoint == "Post-transplant" ~ str_c(patient, " (PoT)"),
                                    timepoint == "Remission" ~ str_c(patient, " (Rem)"),
                                    timepoint == "Relapse" ~ str_c(patient, " (Rel)"),
                                    timepoint == "Normal" ~ str_c(patient, " (N)"))) %>%
    mutate(display_name = factor(display_name, levels = display_name, ordered = TRUE))


  color_tbl <- patient_sample_names_tbl %>%
    filter(timepoint != "Normal") %>%
    select(sample) %>% arrange(sample) %>%
    mutate(sample_n = row_number(),
           my_color_100 = viridis(n(), alpha = 1)[row_number()],
           my_color_75 = viridis(n(), alpha = .75)[row_number()],
           my_color_50 = viridis(n(), alpha = .50)[row_number()],
           my_color_25 = viridis(n(), alpha = .25)[row_number()]) %>%
    bind_rows(patient_sample_names_tbl %>%
                filter(timepoint == "Normal") %>%
                select(sample) %>%
                mutate(sample_n = 0,
                       my_color_100 = "#000000", my_color_75 = "#636363",
                       my_color_50 = "#bdbdbd", my_color_25 = "#f0f0f0")) %>%
    mutate(my_shape = 16)

  patient_sample_names_tbl <- patient_sample_names_tbl %>%
    left_join(color_tbl, by = "sample")

  rm(color_tbl)

  # 1000G samples from 10x website

  normal_samples_tbl <- tibble(patient = c("NA12878", "NA19240"),
                               sample = c("NA12878", "NA19240"),
                               timepoint = "Normal",
                               vcf_column = NA,
                               sorted = FALSE,
                               cnv_maf_status = FALSE) %>%
    mutate(timepoint = factor(timepoint,
                              levels = c("SMM",
                                         "Primary",
                                         "Pre-transplant",
                                         "Post-transplant",
                                         "Remission",
                                         "Relapse-1",
                                         "Normal"),
                              labels = c("SMM",
                                         "Primary",
                                         "Pre-transplant",
                                         "Post-transplant",
                                         "Remission",
                                         "Relapse",
                                         "Normal"),
                              ordered = TRUE)) %>%
    mutate(sample_n = 0,
           my_color_100 = "#fc9272", my_color_75 = "#fc9272",
           my_color_50 = "#fc9272", my_color_25 = "#fc9272",
           my_shape = c(3,4))

  # chromosomes

  chromosome_tbl <- tibble(chromosome = factor(str_c("chr", seq(1:22)),
                                               levels = str_c("chr", seq(1:22)),
                                               ordered = TRUE))

  # chromosome lengths

  chromosome_length_tbl <- read_tsv("data/genome.fa.fai",
                                    col_names = c("contig", "size", "location",
                                                  "basesPerLine",
                                                  "bytesPerLine"),
                                    col_types = c("cidii")) %>%
    filter(contig %in% chromosome_tbl$chromosome) %>%
    mutate(contig = factor(contig,
                           levels = str_c("chr", seq(1:22)),
                           ordered = TRUE))
  # cannot use integer for location because .Machine$integer.max = 2147483647

  # centromeres

  centromere_column_names <- c("Region-Name", "Chromosome", "Chromosome-Start",
                               "Chromosome-Stop", "Scaffold-Role",
                               "Scaffold-GenBank-Accn", "Scaffold-RefSeq-Accn",
                               "Assembly-Unit")

  centromere_tbl <- read_tsv("data/centromeres.tsv", comment = "#",
                             col_names = centromere_column_names,
                             col_types = c("cciicccc")) %>%
    filter(`Scaffold-Role` == "CEN",
           Chromosome  %in% chromosome_tbl$chromosome) %>%
    mutate(Chromosome = factor(Chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  rm(centromere_column_names)

  # protein coding gene annotations
  protein_coding_genes_tbl <- read_tsv("data/gene_annotations.protein_coding.gtf",
                                       col_types = c("ccciicccc"))

  # 10x_longranger_summaries

  lr_summary_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    if ( is.null(lr_summary_tbl) ) {
      lr_summary_tbl <- read_csv(str_c("data/10x_longranger_summaries/",
                                       sample_id, "/summary.csv"),
                                 col_types = "ccidddidddiiddiidddd-ddddddddii")
    } else {
      new_row <- read_csv(str_c("data/10x_longranger_summaries/",
                                sample_id, "/summary.csv"),
                          col_types = "ccidddidddiiddiidddd-ddddddddii")
      lr_summary_tbl <- bind_rows(lr_summary_tbl, new_row)
      rm(new_row)
    }
  }

  lr_summary_tbl <- bind_cols(patient_sample_names_tbl,
                              lr_summary_tbl)

  # for 1000G sampleps from 10x

  lr_summary_1000G_tbl <- NULL

  for (sample_id in normal_samples_tbl$sample) {
    if ( is.null(lr_summary_1000G_tbl) ) {
      lr_summary_1000G_tbl <- read_csv(str_c("data/10x_longranger_summaries/",
                                       sample_id, "/summary.csv"),
                                 col_types = "ccidddidddiiddiidddd-ddddddddii")
    } else {
      new_row <- read_csv(str_c("data/10x_longranger_summaries/",
                                sample_id, "/summary.csv"),
                          col_types = "ccidddidddiiddiidddd-ddddddddii")
      lr_summary_1000G_tbl <- bind_rows(lr_summary_1000G_tbl, new_row)
      rm(new_row)
    }
  }

  lr_summary_1000G_tbl <- bind_cols(normal_samples_tbl,
                                    lr_summary_1000G_tbl)

  # phase set summary

  phase_set_summary_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    for (chromosome in chromosome_tbl$chromosome) {
      if ( is.null(phase_set_summary_tbl) ) {
        phase_set_summary_tbl <- read_tsv(str_c("data/summarize/",
                                                sample_id,
                                                "/",
                                                str_c(sample_id,
                                                      chromosome,
                                                      "phase_set_summary.tsv",
                                                      sep = ".")),
                                          col_types = "dddddddddd")  %>%
          mutate(sample = sample_id, chromosome = chromosome)
      } else {
        new_row <- read_tsv(str_c("data/summarize/",
                                  sample_id,
                                  "/",
                                  str_c(sample_id,
                                        chromosome,
                                        "phase_set_summary.tsv",
                                        sep = ".")),
                            col_types = "dddddddddd") %>%
          mutate(sample = sample_id, chromosome = chromosome)
        phase_set_summary_tbl <- bind_rows(phase_set_summary_tbl, new_row)
        rm(new_row)
      }
    }
  }

  phase_set_summary_tbl <- phase_set_summary_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(chromosome = factor(chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  # phase sets

  phase_sets_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    for (chromosome in chromosome_tbl$chromosome) {
      if ( is.null(phase_sets_tbl) ) {
        phase_sets_tbl <- read_tsv(str_c("data/summarize/",
                                         sample_id,
                                         "/",
                                         str_c(sample_id,
                                               chromosome,
                                               "phase_sets.tsv",
                                               sep = ".")),
                                   col_types = "cciiiiiiiii")  %>%
          mutate(sample = sample_id, chromosome = chromosome)
      } else {
        new_row <- read_tsv(str_c("data/summarize/",
                                  sample_id,
                                  "/",
                                  str_c(sample_id,
                                        chromosome,
                                        "phase_sets.tsv",
                                        sep = ".")),
                            col_types = "cciiiiiiiii") %>%
          mutate(sample = sample_id, chromosome = chromosome)
        phase_sets_tbl <- bind_rows(phase_sets_tbl, new_row)
        rm(new_row)
      }
    }
  }

  phase_sets_tbl <- phase_sets_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(chromosome = factor(chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  # extend stats

  extend_stats_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    for (extended_by in patient_sample_names_tbl$sample) {
      for (chromosome in chromosome_tbl$chromosome) {
        file_path1 <- str_c("data/extend/", sample_id, "/")
        file_path2 <- str_c(sample_id, chromosome, "extended_by_", sep = ".")
        file_path3 <- str_c(extended_by, "extend_stats.tsv", sep = ".")
        file_path <- str_c(file_path1, file_path2, file_path3)
        if (file.exists(file_path)) {
          if ( is.null(extend_stats_tbl) ) {
            extend_stats_tbl <- read_tsv(file_path,
                                         col_types = "ccddccdddddddddddci")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
          } else {
            new_row <- read_tsv(file_path, col_types = "ccddccdddddddddddci")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
            extend_stats_tbl <- bind_rows(extend_stats_tbl, new_row)
            rm(new_row)
          }
        }
      }
    }
  }

  extend_stats_tbl <- extend_stats_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(ps1_Chromosome = factor(ps1_Chromosome,
                                   levels = str_c("chr", seq(1:22)),
                                   ordered = TRUE)) %>%
    mutate(ps2_Chromosome = factor(ps2_Chromosome,
                                   levels = str_c("chr", seq(1:22)),
                                   ordered = TRUE))

  rm(file_path1, file_path2, file_path3, file_path, extended_by)

  # extend phase sets

  extend_phase_sets_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    for (extended_by in patient_sample_names_tbl$sample) {
      for (chromosome in chromosome_tbl$chromosome) {
        file_path1 <- str_c("data/extend/", sample_id, "/")
        file_path2 <- str_c(sample_id, chromosome, "extended_by_", sep = ".")
        file_path3 <- str_c(extended_by, "extended_phase_sets.tsv", sep = ".")
        file_path <- str_c(file_path1, file_path2, file_path3)
        if (file.exists(file_path)) {
          if ( is.null(extend_phase_sets_tbl) ) {
            extend_phase_sets_tbl <- read_tsv(file_path,
                                              col_types = "cciiiiiiiiiicli")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
          } else {
            new_row <- read_tsv(file_path, col_types = "cciiiiiiiiiicli")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
            extend_phase_sets_tbl <- bind_rows(extend_phase_sets_tbl, new_row)
            rm(new_row)
          }
        }
      }
    }
  }

  extend_phase_sets_tbl <- extend_phase_sets_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(chr = factor(chr,
                        levels = str_c("chr", seq(1:22)),
                        ordered = TRUE))

  rm(file_path1, file_path2, file_path3, file_path, extended_by)

  # coverage

  coverage_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    for (chromosome in chromosome_tbl$chromosome) {
      if ( is.null(coverage_tbl) ) {
        coverage_tbl <- read_tsv(str_c("data/coverage_n_variants/",
                                       sample_id,
                                       "/",
                                       str_c(sample_id,
                                             chromosome,
                                             "coverage_n_variants.tsv",
                                             sep = ".")),
                                 skip = 1,
                                 col_names = c("chromosome", "start", "stop",
                                               "coverage",
                                               "n_phased_heterozygotes",
                                               "n_heterozygotes",
                                               "n_homozygotes"),
                                 col_types = "ciidiii") %>%
          mutate(sample = sample_id, chromosome = chromosome)
      } else {
        new_row <- read_tsv(str_c("data/coverage_n_variants/",
                                  sample_id,
                                  "/",
                                  str_c(sample_id,
                                        chromosome,
                                        "coverage_n_variants.tsv",
                                        sep = ".")),
                            skip = 1,
                            col_names = c("chromosome", "start", "stop",
                                          "coverage",
                                          "n_phased_heterozygotes",
                                          "n_heterozygotes",
                                          "n_homozygotes"),
                            col_types = "ciidiii") %>%
          mutate(sample = sample_id, chromosome = chromosome)
        coverage_tbl <- bind_rows(coverage_tbl, new_row)
        rm(new_row)
      }
    }
  }

  coverage_tbl <- coverage_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(chromosome = factor(chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  # Sorted WGS CNV

  cnv_tbl <- NULL

  for (sample_id_lr in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    if (sample_id_lr == "27522_3") {
      sample_id <- "27522_4" # CNV time point (tp4) is 0.3 years later than lr WGS (tp3)
    } else {
      sample_id <- sample_id_lr
    }
    if ( is.null(cnv_tbl) ) {
      cnv_tbl <- read_tsv(str_c("data/sorted_wgs_cnv/", sample_id,
                                ".lambda90.noscale.CNV.bicseq2.sorted.WGS.tsv"),
                          col_types = "ciid") %>%
        mutate(sample = sample_id_lr)
    } else {
      new_cnv <- read_tsv(str_c("data/sorted_wgs_cnv/", sample_id,
                                ".lambda90.noscale.CNV.bicseq2.sorted.WGS.tsv"),
                          col_types = "ciid") %>%
        mutate(sample = sample_id_lr)
      cnv_tbl <- bind_rows(cnv_tbl, new_cnv)
      rm(new_cnv)
    }
  }

  cnv_tbl <- cnv_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(chrom = factor(chrom,
                          levels = str_c("chr", seq(1:22)),
                          ordered = TRUE))

  # Sorted WGS MAF

  maf_tbl <- NULL
  maf_col_types = "c--ccddcccccc--cc-----------------cccc------------c---------cccc-----------c----------------c------------------------------"

  for (sample_id_lr in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    if (sample_id_lr == "27522_3") {
      sample_id <- "27522_4" # MAF time point (tp4) is 0.3 years later than lr WGS (tp3)
    } else {
      sample_id <- sample_id_lr
    }
    if ( is.null(maf_tbl) ) {
      maf_tbl <- read_tsv(str_c("data/sorted_wgs_maf/", sample_id,
                                ".withmutect.maf"),
                          col_types = maf_col_types,
                          skip = 1) %>%
        mutate(sample = sample_id_lr)
    } else {
      new_maf <- read_tsv(str_c("data/sorted_wgs_maf/", sample_id,
                                ".withmutect.maf"),
                          col_types = maf_col_types,
                          skip = 1) %>%
        mutate(sample = sample_id_lr)
      maf_tbl <- bind_rows(maf_tbl, new_maf)
      rm(new_maf)
    }
  }

  rm(maf_col_types)

  maf_tbl <- maf_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample") %>%
    mutate(Chromosome = factor(Chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  # somatic events (WGS samples only)
  {
    # somatic barcodes variants

    barcodes_variants_mapq20_tbl <- list() # use a list to keep samples separate
                                           # because combined files are too big

    for (sample_id in patient_sample_names_tbl %>%
         filter(cnv_maf_status) %>% pull(sample)) {
      barcodes_variants_mapq20_tbl[[sample_id]] <- NULL
      for (chromosome in chromosome_tbl$chromosome) {
        if (is.null(barcodes_variants_mapq20_tbl[[sample_id]])) {
          barcodes_variants_mapq20_tbl[[sample_id]] <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.barcodes_variants.tsv", sep = ".")),
            col_types = "ccciccccccc") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.barcodes_variants.tsv", sep = ".")),
            col_types = "ccciccccccc") %>%
            mutate(sample = sample_id)
          barcodes_variants_mapq20_tbl[[sample_id]] <- bind_rows(barcodes_variants_mapq20_tbl[[sample_id]], q20)
          rm(q20)
        }
      }

      barcodes_variants_mapq20_tbl[[sample_id]] <- barcodes_variants_mapq20_tbl[[sample_id]] %>%
        left_join(patient_sample_names_tbl, by = "sample") %>%
        mutate(Chromosome = factor(Chromosome,
                                   levels = str_c("chr", seq(1:22)),
                                   ordered = TRUE))

    }

    # somatic phasing variants

    phasing_variants_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(cnv_maf_status) %>% pull(sample)) {
      for (chromosome in chromosome_tbl$chromosome) {
        if (is.null(phasing_variants_mapq20_tbl)) {
          phasing_variants_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.phasing_variants.tsv", sep = ".")),
            col_types = "cciccciccddddiiiiiddddiiiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.phasing_variants.tsv", sep = ".")),
            col_types = "cciccciccddddiiiiiddddiiiii") %>%
            mutate(sample = sample_id)
          phasing_variants_mapq20_tbl <- bind_rows(phasing_variants_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    phasing_variants_mapq20_tbl <- phasing_variants_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      filter(Chromosome %in% str_c("chr", seq(1:22))) %>%
      mutate(Chromosome = factor(Chromosome,
                                 levels = str_c("chr", seq(1:22)),
                                 ordered = TRUE)) %>%
      mutate(enough_coverage = (n_ALT_H1 + n_ALT_H2 >= 10 | barcode_ALT_H1 + barcode_ALT_H2 > 0),
             phased_by_linked_alleles = case_when(n_ALT_H1 + n_ALT_H2 < 10 ~ "NC",
                                                  n_ALT_H1 + n_ALT_H2 >= 10 & pct_ALT_on_H1 >= phase_proportion ~ "H1",
                                                  n_ALT_H1 + n_ALT_H2 >= 10 & pct_ALT_on_H2 >= phase_proportion ~ "H2",
                                                  TRUE ~ "NP"),
             phased_by_barcodes = case_when(barcode_ALT_H1 + barcode_ALT_H2 == 0 ~ "NC",
                                            barcode_ALT_H1 + barcode_ALT_H2 > 0 & pct_barcode_ALT_on_H1 == 1 ~ "H1",
                                            barcode_ALT_H1 + barcode_ALT_H2 > 0 & pct_barcode_ALT_on_H2 == 1 ~ "H2",
                                            TRUE ~ "NP"),
             phased_by = case_when(!enough_coverage ~ "NC",
                                   phased_by_linked_alleles == "NC" & phased_by_barcodes %in% c("H1", "H2") ~ "BC",
                                   phased_by_barcodes == "NC" & phased_by_linked_alleles %in% c("H1", "H2") ~ "LA",
                                   phased_by_linked_alleles == "NP" & phased_by_barcodes %in% c("H1", "H2") ~ "BC",
                                   phased_by_barcodes == "NP" & phased_by_linked_alleles %in% c("H1", "H2") ~ "LA",
                                   phased_by_linked_alleles %in% c("H1", "H2") & phased_by_barcodes %in% c("H1", "H2")  & phased_by_linked_alleles == phased_by_barcodes ~ "Both (agree)",
                                   phased_by_linked_alleles %in% c("H1", "H2") & phased_by_barcodes %in% c("H1", "H2")  & phased_by_linked_alleles != phased_by_barcodes ~ "Both (conflict)",
                                   TRUE ~ "NP"),
             phased = case_when(!enough_coverage ~ "Not enough coverage",
                                phased_by == "Both (conflict)" ~ "Conflict",
                                phased_by %in% c("BC", "LA", "Both (agree)") ~ "Phased",
                                TRUE ~ "Not phased"))

    # somatic somatic per phase set

    somatic_per_phase_set_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(cnv_maf_status) %>% pull(sample)) {
      for (chromosome in chromosome_tbl$chromosome) {
        if (is.null(somatic_per_phase_set_mapq20_tbl)) {
          somatic_per_phase_set_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.somatic_per_phase_set.tsv", sep = ".")),
            col_types = "cciiiiiiiiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.somatic_per_phase_set.tsv", sep = ".")),
            col_types = "cciiiiiiiiii") %>%
            mutate(sample = sample_id)
          somatic_per_phase_set_mapq20_tbl <- bind_rows(somatic_per_phase_set_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    somatic_per_phase_set_mapq20_tbl <- somatic_per_phase_set_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      mutate(chrom = factor(chrom,
                            levels = str_c("chr", seq(1:22)),
                            ordered = TRUE))

    # somatic somatic per phase set

    variant_pairs_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(cnv_maf_status) %>% pull(sample)) {
      for (chromosome in chromosome_tbl$chromosome) {
        if (is.null(variant_pairs_mapq20_tbl)) {
          variant_pairs_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.variant_pairs.with_cnv.tsv", sep = ".")),
            col_types = "ccccccciiiiddd") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.variant_pairs.with_cnv.tsv", sep = ".")),
            col_types = "ccccccciiiiddd") %>%
            mutate(sample = sample_id)
          variant_pairs_mapq20_tbl <- bind_rows(variant_pairs_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    variant_pairs_mapq20_tbl <- variant_pairs_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      separate(Variant1, into = c("chromosome1", "position1",
                                  "reference1", "alternate1"),
               sep = ":", remove = FALSE) %>%
      separate(Variant2, into = c("chromosome2", "position2",
                                  "reference2", "alternate2"),
               sep = ":", remove = FALSE) %>%
      mutate(chromosome1 = factor(chromosome1,
                                  levels = str_c("chr", seq(1:22)),
                                  ordered = TRUE),
             chromosome2 = factor(chromosome2,
                                  levels = str_c("chr", seq(1:22)),
                                  ordered = TRUE)) %>%
      mutate(position1 = as.numeric(position1),
             position2 = as.numeric(position2)) %>%
      rowwise() %>%
      mutate(distance_between_variants = position2 - position1,
             n_overlapping_barcodes = sum(n_bx_overlap_00, n_bx_overlap_01,
                                          n_bx_overlap_10, n_bx_overlap_11),
             n_barcodes_with_mutation = sum(n_bx_overlap_01, n_bx_overlap_10,
                                            n_bx_overlap_11)) %>%
      ungroup()

    # somatic somatic barcodes (sombx)

    sombx_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(cnv_maf_status) %>% pull(sample)) {
      if (is.null(sombx_mapq20_tbl)) {
        sombx_mapq20_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id,
                "mapq20.sombx.tsv", sep = ".")),
          col_types = "ccccccccccc") %>%
          mutate(sample = sample_id)
      } else {
        q20 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id,
                "mapq20.sombx.tsv", sep = ".")),
          col_types = "ccccccccccc") %>%
          mutate(sample = sample_id)
        sombx_mapq20_tbl <- bind_rows(sombx_mapq20_tbl, q20)
        rm(q20)
      }
    }

    sombx_mapq20_tbl <- sombx_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      mutate(chromosome = factor(chromosome,
                                 levels = str_c("chr", seq(1:22)),
                                 ordered = TRUE))
    }

  # somatic events (all tumor samples, important genes only)
  {
    # (all) somatic barcodes variants

    important_mutations_tbl <- read_delim("data/important_genes.mutations.with_annotation.txt",
                                          delim = ":",
                                          col_names = c("gene", "chr", "pos",
                                                        "ref", "alt", "protein"),
                                          col_types = "cciccc") %>%
      mutate(chr = factor(chr,
                          levels = str_c("chr", seq(1:22)),
                          ordered = TRUE))

    important_mutations_vaf_tbl <- read_tsv("data/important_mutations_vaf.txt",
                                          col_types = "cccccicccdc") %>%
      select(-c("timepoint")) %>%
      mutate(chr = factor(chr,
                          levels = str_c("chr", seq(1:22)),
                          ordered = TRUE)) %>%
      left_join(patient_sample_names_tbl,
                by = c("sample", "patient"))

    barcodes_variants_all_mapq20_tbl <- list() # use a list to keep samples separate
                                               # because combined files are too big

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {

      barcodes_variants_all_mapq20_tbl[[sample_id]] <- NULL
      for (chromosome in unique(important_mutations_tbl$chr)) {
        if (is.null(barcodes_variants_all_mapq20_tbl[[sample_id]])) {
          barcodes_variants_all_mapq20_tbl[[sample_id]] <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.barcodes_variants.tsv", sep = ".")),
            col_types = "ccciccccccc") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.barcodes_variants.tsv", sep = ".")),
            col_types = "ccciccccccc") %>%
            mutate(sample = sample_id)
          barcodes_variants_all_mapq20_tbl[[sample_id]] <- bind_rows(barcodes_variants_all_mapq20_tbl[[sample_id]], q20)
          rm(q20)
        }
      }

      barcodes_variants_all_mapq20_tbl[[sample_id]] <- barcodes_variants_all_mapq20_tbl[[sample_id]] %>%
        left_join(patient_sample_names_tbl, by = "sample") %>%
        mutate(Chromosome = factor(Chromosome,
                                   levels = str_c("chr", seq(1:22)),
                                   ordered = TRUE))

    }

    # (all) somatic phasing variants

    phasing_variants_all_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      for (chromosome in unique(important_mutations_tbl$chr)) {
        if (is.null(phasing_variants_all_mapq20_tbl)) {
          phasing_variants_all_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.phasing_variants.tsv", sep = ".")),
            col_types = "cciccciccddddiiiiiddddiiiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.phasing_variants.tsv", sep = ".")),
            col_types = "cciccciccddddiiiiiddddiiiii") %>%
            mutate(sample = sample_id)
          phasing_variants_all_mapq20_tbl <- bind_rows(phasing_variants_all_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    phasing_variants_all_mapq20_tbl <- phasing_variants_all_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      filter(Chromosome %in% str_c("chr", seq(1:22))) %>%
      mutate(Chromosome = factor(Chromosome,
                                 levels = str_c("chr", seq(1:22)),
                                 ordered = TRUE)) %>%
      mutate(enough_coverage = (n_ALT_H1 + n_ALT_H2 >= 10 | barcode_ALT_H1 + barcode_ALT_H2 > 0),
             phased_by_linked_alleles = case_when(n_ALT_H1 + n_ALT_H2 < 10 ~ "NC",
                                                  n_ALT_H1 + n_ALT_H2 >= 10 & pct_ALT_on_H1 >= phase_proportion ~ "H1",
                                                  n_ALT_H1 + n_ALT_H2 >= 10 & pct_ALT_on_H2 >= phase_proportion ~ "H2",
                                                  TRUE ~ "NP"),
             phased_by_barcodes = case_when(barcode_ALT_H1 + barcode_ALT_H2 == 0 ~ "NC",
                                            barcode_ALT_H1 + barcode_ALT_H2 > 0 & pct_barcode_ALT_on_H1 == 1 ~ "H1",
                                            barcode_ALT_H1 + barcode_ALT_H2 > 0 & pct_barcode_ALT_on_H2 == 1 ~ "H2",
                                            TRUE ~ "NP"),
             phased_by = case_when(!enough_coverage ~ "NC",
                                   phased_by_linked_alleles == "NC" & phased_by_barcodes %in% c("H1", "H2") ~ "BC",
                                   phased_by_barcodes == "NC" & phased_by_linked_alleles %in% c("H1", "H2") ~ "LA",
                                   phased_by_linked_alleles == "NP" & phased_by_barcodes %in% c("H1", "H2") ~ "BC",
                                   phased_by_barcodes == "NP" & phased_by_linked_alleles %in% c("H1", "H2") ~ "LA",
                                   phased_by_linked_alleles %in% c("H1", "H2") & phased_by_barcodes %in% c("H1", "H2")  & phased_by_linked_alleles == phased_by_barcodes ~ "Both (agree)",
                                   phased_by_linked_alleles %in% c("H1", "H2") & phased_by_barcodes %in% c("H1", "H2")  & phased_by_linked_alleles != phased_by_barcodes ~ "Both (conflict)",
                                   TRUE ~ "NP"),
               phased = case_when(!enough_coverage ~ "Not enough coverage",
                                  phased_by == "Both (conflict)" ~ "Conflict",
                                  phased_by %in% c("BC", "LA", "Both (agree)") ~ "Phased",
                                  TRUE ~ "Not phased"))

    # (all) somatic somatic per phase set

    somatic_per_phase_set_all_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      for (chromosome in unique(important_mutations_tbl$chr)) {
        if (is.null(somatic_per_phase_set_all_mapq20_tbl)) {
          somatic_per_phase_set_all_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.somatic_per_phase_set.tsv", sep = ".")),
            col_types = "cciiiiiiiiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.somatic_per_phase_set.tsv", sep = ".")),
            col_types = "cciiiiiiiiii") %>%
            mutate(sample = sample_id)
          somatic_per_phase_set_all_mapq20_tbl <- bind_rows(somatic_per_phase_set_all_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    somatic_per_phase_set_all_mapq20_tbl <- somatic_per_phase_set_all_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      mutate(chrom = factor(chrom,
                            levels = str_c("chr", seq(1:22)),
                            ordered = TRUE))


    # (all) somatic somatic per phase set

    variant_pairs_all_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      for (chromosome in unique(important_mutations_tbl$chr)) {
        if (is.null(variant_pairs_all_mapq20_tbl)) {
          variant_pairs_all_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.variant_pairs.tsv", sep = ".")),
            col_types = "ccccccciiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_all/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.variant_pairs.tsv", sep = ".")),
            col_types = "ccccccciiii") %>%
            mutate(sample = sample_id)
          variant_pairs_all_mapq20_tbl <- bind_rows(variant_pairs_all_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    variant_pairs_all_mapq20_tbl <- variant_pairs_all_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      separate(Variant1, into = c("chromosome1", "position1",
                                  "reference1", "alternate1"),
               sep = ":", remove = FALSE) %>%
      separate(Variant2, into = c("chromosome2", "position2",
                                  "reference2", "alternate2"),
               sep = ":", remove = FALSE) %>%
      mutate(chromosome1 = factor(chromosome1,
                                  levels = str_c("chr", seq(1:22)),
                                  ordered = TRUE),
             chromosome2 = factor(chromosome2,
                                  levels = str_c("chr", seq(1:22)),
                                  ordered = TRUE)) %>%
      mutate(position1 = as.numeric(position1),
             position2 = as.numeric(position2)) %>%
      rowwise() %>%
      mutate(distance_between_variants = position2 - position1,
             n_overlapping_barcodes = sum(n_bx_overlap_00, n_bx_overlap_01,
                                          n_bx_overlap_10, n_bx_overlap_11),
             n_barcodes_with_mutation = sum(n_bx_overlap_01, n_bx_overlap_10,
                                            n_bx_overlap_11)) %>%
      ungroup()

    # (all) somatic somatic barcodes (sombx)

    sombx_all_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      if (is.null(sombx_all_mapq20_tbl)) {
        sombx_all_mapq20_tbl <- read_tsv(str_c(
          str_c("data/somatic_all/", sample_id, "/"),
          str_c(sample_id,
                "mapq20.sombx.tsv", sep = ".")),
          col_types = "ccccccccccc") %>%
          mutate(sample = sample_id)
      } else {
        q20 <- read_tsv(str_c(
          str_c("data/somatic_all/", sample_id, "/"),
          str_c(sample_id,
                "mapq20.sombx.tsv", sep = ".")),
          col_types = "ccccccccccc") %>%
          mutate(sample = sample_id)
        sombx_all_mapq20_tbl <- bind_rows(sombx_all_mapq20_tbl, q20)
        rm(q20)
      }
    }

    sombx_all_mapq20_tbl <- sombx_all_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      mutate(chromosome = factor(chromosome,
                                 levels = str_c("chr", seq(1:22)),
                                 ordered = TRUE))
  }

  # somatic events (all tumor samples, driver_mm genes only)
  {
    # (drivers) somatic barcodes variants

    driver_mutations_tbl <- read_delim("data/driver_mm_genes.mutations.with_annotation.txt",
                                       delim = ":",
                                       col_names = c("gene", "chr", "pos",
                                                     "ref", "alt", "protein"),
                                       col_types = "cciccc") %>%
      filter(chr %in% str_c("chr", seq(1:22))) %>%
      mutate(chr = factor(chr,
                          levels = str_c("chr", seq(1:22)),
                          ordered = TRUE))

    driver_mutations_vaf_tbl <- read_tsv("data/driver_mm_mutations_vaf.txt",
                                         col_types = "cccccicccdc") %>%
      select(-c("timepoint")) %>%
      filter(chr %in% str_c("chr", seq(1:22))) %>%
      mutate(chr = factor(chr,
                          levels = str_c("chr", seq(1:22)),
                          ordered = TRUE)) %>%
      left_join(patient_sample_names_tbl,
                by = c("sample", "patient"))

    barcodes_variants_driver_mapq20_tbl <- list() # use a list to keep samples separate
    # because combined files are too big

    for (sample_id in patient_sample_names_tbl %>% filter(timepoint != "Normal") %>% pull(sample)) {
      barcodes_variants_driver_mapq20_tbl[[sample_id]] <- NULL
      for (chromosome in unique(driver_mutations_tbl$chr)) {
        if (is.null(barcodes_variants_driver_mapq20_tbl[[sample_id]])) {
          barcodes_variants_driver_mapq20_tbl[[sample_id]] <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.barcodes_variants.tsv", sep = ".")),
            col_types = "ccciccccccc") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.barcodes_variants.tsv", sep = ".")),
            col_types = "ccciccccccc") %>%
            mutate(sample = sample_id)
          barcodes_variants_driver_mapq20_tbl[[sample_id]] <- bind_rows(barcodes_variants_driver_mapq20_tbl[[sample_id]], q20)
          rm(q20)
        }
      }

      barcodes_variants_driver_mapq20_tbl[[sample_id]] <- barcodes_variants_driver_mapq20_tbl[[sample_id]] %>%
        left_join(patient_sample_names_tbl, by = "sample") %>%
        filter(Chromosome %in% str_c("chr", seq(1:22))) %>%
        mutate(Chromosome = factor(Chromosome,
                                   levels = str_c("chr", seq(1:22)),
                                   ordered = TRUE))

    }

    # (driver) somatic phasing variants

    phasing_variants_driver_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      for (chromosome in unique(driver_mutations_tbl$chr)) {
        if (is.null(phasing_variants_driver_mapq20_tbl)) {
          phasing_variants_driver_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.phasing_variants.tsv", sep = ".")),
            col_types = "cciccciccddddiiiiiddddiiiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.phasing_variants.tsv", sep = ".")),
            col_types = "cciccciccddddiiiiiddddiiiii") %>%
            mutate(sample = sample_id)
          phasing_variants_driver_mapq20_tbl <- bind_rows(phasing_variants_driver_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    phasing_variants_driver_mapq20_tbl <- phasing_variants_driver_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      filter(Chromosome %in% str_c("chr", seq(1:22))) %>%
      mutate(Chromosome = factor(Chromosome,
                                 levels = str_c("chr", seq(1:22)),
                                 ordered = TRUE)) %>%
      mutate(enough_coverage = (n_ALT_H1 + n_ALT_H2 >= 10 | barcode_ALT_H1 + barcode_ALT_H2 > 0),
             phased_by_linked_alleles = case_when(n_ALT_H1 + n_ALT_H2 < 10 ~ "NC",
                                                  n_ALT_H1 + n_ALT_H2 >= 10 & pct_ALT_on_H1 >= phase_proportion ~ "H1",
                                                  n_ALT_H1 + n_ALT_H2 >= 10 & pct_ALT_on_H2 >= phase_proportion ~ "H2",
                                                  TRUE ~ "NP"),
             phased_by_barcodes = case_when(barcode_ALT_H1 + barcode_ALT_H2 == 0 ~ "NC",
                                            barcode_ALT_H1 + barcode_ALT_H2 > 0 & pct_barcode_ALT_on_H1 == 1 ~ "H1",
                                            barcode_ALT_H1 + barcode_ALT_H2 > 0 & pct_barcode_ALT_on_H2 == 1 ~ "H2",
                                            TRUE ~ "NP"),
             phased_by = case_when(!enough_coverage ~ "NC",
                                   phased_by_linked_alleles == "NC" & phased_by_barcodes %in% c("H1", "H2") ~ "BC",
                                   phased_by_barcodes == "NC" & phased_by_linked_alleles %in% c("H1", "H2") ~ "LA",
                                   phased_by_linked_alleles == "NP" & phased_by_barcodes %in% c("H1", "H2") ~ "BC",
                                   phased_by_barcodes == "NP" & phased_by_linked_alleles %in% c("H1", "H2") ~ "LA",
                                   phased_by_linked_alleles %in% c("H1", "H2") & phased_by_barcodes %in% c("H1", "H2")  & phased_by_linked_alleles == phased_by_barcodes ~ "Both (agree)",
                                   phased_by_linked_alleles %in% c("H1", "H2") & phased_by_barcodes %in% c("H1", "H2")  & phased_by_linked_alleles != phased_by_barcodes ~ "Both (conflict)",
                                   TRUE ~ "NP"),
               phased = case_when(!enough_coverage ~ "Not enough coverage",
                                  phased_by == "Both (conflict)" ~ "Conflict",
                                  phased_by %in% c("BC", "LA", "Both (agree)") ~ "Phased",
                                  TRUE ~ "Not phased"))

    # (driver) somatic somatic per phase set

    somatic_per_phase_set_driver_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      for (chromosome in unique(driver_mutations_tbl$chr)) {
        if (is.null(somatic_per_phase_set_driver_mapq20_tbl)) {
          somatic_per_phase_set_driver_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.somatic_per_phase_set.tsv", sep = ".")),
            col_types = "cciiiiiiiiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.somatic_per_phase_set.tsv", sep = ".")),
            col_types = "cciiiiiiiiii") %>%
            mutate(sample = sample_id)
          somatic_per_phase_set_driver_mapq20_tbl <- bind_rows(somatic_per_phase_set_driver_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    somatic_per_phase_set_driver_mapq20_tbl <- somatic_per_phase_set_driver_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      filter(chrom %in% str_c("chr", seq(1:22))) %>%
      mutate(chrom = factor(chrom,
                            levels = str_c("chr", seq(1:22)),
                            ordered = TRUE))

    # (driver) somatic somatic per phase set

    variant_pairs_driver_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      for (chromosome in unique(driver_mutations_tbl$chr)) {
        if (is.null(variant_pairs_driver_mapq20_tbl)) {
          variant_pairs_driver_mapq20_tbl <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.variant_pairs.tsv", sep = ".")),
            col_types = "ccccccciiii") %>%
            mutate(sample = sample_id)
        } else {
          q20 <- read_tsv(str_c(
            str_c("data/somatic_driver/", sample_id, "/"),
            str_c(sample_id, chromosome,
                  "mapq20.variant_pairs.tsv", sep = ".")),
            col_types = "ccccccciiii") %>%
            mutate(sample = sample_id)
          variant_pairs_driver_mapq20_tbl <- bind_rows(variant_pairs_driver_mapq20_tbl, q20)
          rm(q20)
        }
      }
    }

    variant_pairs_driver_mapq20_tbl <- variant_pairs_driver_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      separate(Variant1, into = c("chromosome1", "position1",
                                  "reference1", "alternate1"),
               sep = ":", remove = FALSE) %>%
      separate(Variant2, into = c("chromosome2", "position2",
                                  "reference2", "alternate2"),
               sep = ":", remove = FALSE) %>%
      filter(chromosome1 %in% str_c("chr", seq(1:22))) %>%
      filter(chromosome2 %in% str_c("chr", seq(1:22))) %>%
      mutate(chromosome1 = factor(chromosome1,
                                  levels = str_c("chr", seq(1:22)),
                                  ordered = TRUE),
             chromosome2 = factor(chromosome2,
                                  levels = str_c("chr", seq(1:22)),
                                  ordered = TRUE)) %>%
      mutate(position1 = as.numeric(position1),
             position2 = as.numeric(position2)) %>%
      rowwise() %>%
      mutate(distance_between_variants = position2 - position1,
             n_overlapping_barcodes = sum(n_bx_overlap_00, n_bx_overlap_01,
                                          n_bx_overlap_10, n_bx_overlap_11),
             n_barcodes_with_mutation = sum(n_bx_overlap_01, n_bx_overlap_10,
                                            n_bx_overlap_11)) %>%
      ungroup()

    # (driver) somatic somatic barcodes (sombx)

    sombx_driver_mapq20_tbl <- NULL

    for (sample_id in patient_sample_names_tbl %>%
         filter(timepoint != "Normal") %>% pull(sample)) {
      if (is.null(sombx_driver_mapq20_tbl)) {
        sombx_driver_mapq20_tbl <- read_tsv(str_c(
          str_c("data/somatic_driver/", sample_id, "/"),
          str_c(sample_id,
                "mapq20.sombx.tsv", sep = ".")),
          col_types = "ccccccccccc") %>%
          mutate(sample = sample_id)
      } else {
        q20 <- read_tsv(str_c(
          str_c("data/somatic_driver/", sample_id, "/"),
          str_c(sample_id,
                "mapq20.sombx.tsv", sep = ".")),
          col_types = "ccccccccccc") %>%
          mutate(sample = sample_id)
        sombx_driver_mapq20_tbl <- bind_rows(sombx_driver_mapq20_tbl, q20)
        rm(q20)
      }
    }

    sombx_driver_mapq20_tbl <- sombx_driver_mapq20_tbl %>%
      left_join(patient_sample_names_tbl, by = "sample") %>%
      filter(chromosome %in% str_c("chr", seq(1:22))) %>%
      mutate(chromosome = factor(chromosome,
                                 levels = str_c("chr", seq(1:22)),
                                 ordered = TRUE))
  }

  # IBD segment ancestry

  segment_ancestry_tbl <- NULL

  for (sample_id in c("NA12878")) {
    for (chromosome in chromosome_tbl$chromosome) {
      if ( is.null(segment_ancestry_tbl) ) {
        segment_ancestry_tbl <- read_tsv(str_c("data/ancestry/",
                                               sample_id,
                                               "/",
                                               str_c(sample_id,
                                                     chromosome,
                                                     "IBD_segment_ancestry.tsv",
                                                     sep = ".")),
                                         col_types = "iciiilccc")  %>%
          mutate(sample = sample_id)
      } else {
        new_row <- read_tsv(str_c("data/ancestry/",
                                  sample_id,
                                  "/",
                                  str_c(sample_id,
                                        chromosome,
                                        "IBD_segment_ancestry.tsv",
                                        sep = ".")),
                            col_types = "iciiilccc") %>%
          mutate(sample = sample_id)
        segment_ancestry_tbl <- bind_rows(segment_ancestry_tbl, new_row)
        rm(new_row)
      }
    }
  }

  segment_ancestry_tbl <- segment_ancestry_tbl %>%
    left_join(normal_samples_tbl, by = "sample") %>%
    mutate(Chromosome = factor(Chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  # IBD segment overlap

  segment_overlap_tbl <- NULL

  for (sample_id in c("NA12878")) {
    for (chromosome in chromosome_tbl$chromosome) {
      if ( is.null(segment_overlap_tbl) ) {
        segment_overlap_tbl <- read_tsv(str_c("data/ancestry/",
                                              sample_id,
                                              "/",
                                              str_c(sample_id,
                                                    chromosome,
                                                    "IBD_segment_overlap.tsv",
                                                    sep = ".")),
                                        col_types = "icciciciiddciiidddccccc") %>%
          mutate(sample = sample_id)
      } else {
        new_row <- read_tsv(str_c("data/ancestry/",
                                  sample_id,
                                  "/",
                                  str_c(sample_id,
                                        chromosome,
                                        "IBD_segment_overlap.tsv",
                                        sep = ".")),
                            col_types = "icciciciiddciiidddccccc") %>%
          mutate(sample = sample_id)
        segment_overlap_tbl <- bind_rows(segment_overlap_tbl, new_row)
        rm(new_row)
      }
    }
  }

  segment_overlap_tbl <- segment_overlap_tbl %>%
    left_join(normal_samples_tbl, by = "sample") %>%
    mutate(Chromosome = factor(Chromosome,
                               levels = str_c("chr", seq(1:22)),
                               ordered = TRUE))

  sv_barcodes_tbl <- list()
  #for (sample_id in c("27522_1", "27522_3", "77570")) {
  for (sample_id in c("27522_1", "77570")) {
    sv_barcodes_tbl[[sample_id]] <- read_tsv(str_c("data/SVs/", sample_id, ".sv_reads.tsv"),
                                             col_types = "ccicicdd") %>%
      mutate(chrom = factor(chrom,
                            levels = str_c("chr", seq(1:22)),
                            ordered = TRUE)) %>%
      mutate(read_haplotype = factor(read_haplotype,
                                     levels = c(0,1,2),
                                     ordered = TRUE)) %>%
      mutate(supports_trans = factor(supports_trans,
                                     levels = c(0,1),
                                     ordered = TRUE))
  }

  sv_haplotypes_tbl <- NULL

  sv_haplotypes_tbl <- read_tsv("data/SVs/Manta_translocation_11samples_highconfidence_v033120/haplotype_support.tsv",
                                col_names = c("sample",
                                              "first_chromosome", "first_position",
                                              "second_chromosome", "second_position",
                                              "phase_set_pairs",
                                              "barcode_support_00_11_12_21_22",
                                              "total_barcodes",
                                              "haplotypes_consistent"),
                                col_types = "cciciccil") %>%
    filter(!(sample %in% c("27522_1", "77570"))) %>% # these samples have Manta calls based on sorted/unsorted lrWGS sample, not sorted WGS
    mutate(sample = case_when(sample == "27522_4" ~ "27522_3", # WGS timepoint 27522_4 matches closest to lrWGS timepoint 27522_3
                              TRUE ~ sample)) %>%
    mutate(first_chromosome = factor(first_chromosome,
                                     levels = str_c("chr", seq(1:22)),
                                     ordered = TRUE),
           second_chromosome = factor(second_chromosome,
                                      levels = str_c("chr", seq(1:22)),
                                      ordered = TRUE)) %>%
    left_join(patient_sample_names_tbl, by = "sample")


  # clean up

  rm(sample_id)
  rm(sample_id_lr)
  rm(chromosome)

  # save current session info
  dir.create("session_info", recursive = TRUE, showWarnings = FALSE)
  sink(str_c("session_info/session_info.", last_updated, ".txt"))
  print(devtools::session_info())
  sink()

  # save current objects to load next time
  save(list = ls(pattern = "_tbl"), file = str_c(input_data_path_str), envir = .GlobalEnv)
  rm(last_updated, input_data_path_str)

}

# Steven Foltz (envest) January/February 2020
