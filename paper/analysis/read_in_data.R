################################################################################
# Read in all data for analysis
################################################################################

# shh! (this is a library)
library(tidyverse)

# Load previous input data?
# We do not save .Rdata at end of session. This is a resusable input
# data object that is re-loaded at the beginning of each session.

last_updated <- "2020-01-10"
input_data_path_str <- str_c("data/collected_input_objects.", last_updated, ".RData")
if (TRUE) {
  load(input_data_path_str)
  rm(last_updated, input_data_path_str)
  print("Input data tables loaded from last updated .RData file.")
} else {

  # patient and sample names

  patient_sample_names_tbl <- read_tsv("data/patient_sample_names.tsv",
                                       col_types = "ccccll")

  # chromosomes

  chromosome_tbl <- tibble(chromosome = str_c("chr", seq(1,22)))

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
    left_join(patient_sample_names_tbl, by = "sample")

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
    left_join(patient_sample_names_tbl, by = "sample")

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
                                         col_types = "ccddccddddddddddcd")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
          } else {
            new_row <- read_tsv(file_path, col_types = "ccddccddddddddddcd")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
            extend_stats_tbl <- bind_rows(extend_stats_tbl, new_row)
            rm(new_row)
          }
        }
      }
    }
  }

  extend_stats_tbl <- extend_stats_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

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
                                              col_types = "cciiiiiiiiiicl")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
          } else {
            new_row <- read_tsv(file_path, col_types = "cciiiiiiiiiicl")  %>%
              mutate(sample = sample_id, extended_by = extended_by)
            extend_phase_sets_tbl <- bind_rows(extend_phase_sets_tbl, new_row)
            rm(new_row)
          }
        }
      }
    }
  }

  extend_phase_sets_tbl <- extend_phase_sets_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

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
                                 col_types = "ciidiii")  %>%
          mutate(sample = sample_id, chromosome = chromosome)
      } else {
        new_row <- read_tsv(str_c("data/coverage_n_variants/",
                                  sample_id,
                                  "/",
                                  str_c(sample_id,
                                        chromosome,
                                        "coverage_n_variants.tsv",
                                        sep = ".")),
                            col_types = "ciidiii") %>%
          mutate(sample = sample_id, chromosome = chromosome)
        coverage_tbl <- bind_rows(coverage_tbl, new_row)
        rm(new_row)
      }
    }
  }

  coverage_tbl <- coverage_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

  # somatic barcodes variants

  barcodes_variants_mapq0_tbl <- list() # use a list to keep samples separate
  barcodes_variants_mapq20_tbl <- list() # because combined files are too big

  for (sample_id in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    barcodes_variants_mapq0_tbl[[sample_id]] <- NULL
    barcodes_variants_mapq20_tbl[[sample_id]] <- NULL
    for (chromosome in chromosome_tbl$chromosome) {
      if (is.null(barcodes_variants_mapq0_tbl[[sample_id]])) {
        barcodes_variants_mapq0_tbl[[sample_id]] <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.barcodes_variants.tsv", sep = ".")),
          col_types = "ccciccccccc") %>%
          mutate(sample = sample_id)
        barcodes_variants_mapq20_tbl[[sample_id]] <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.barcodes_variants.tsv", sep = ".")),
          col_types = "ccciccccccc") %>%
          mutate(sample = sample_id)
      } else {
        q0 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.barcodes_variants.tsv", sep = ".")),
          col_types = "ccciccccccc") %>%
          mutate(sample = sample_id)
        q20 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.barcodes_variants.tsv", sep = ".")),
          col_types = "ccciccccccc") %>%
          mutate(sample = sample_id)
        barcodes_variants_mapq0_tbl[[sample_id]] <- bind_rows(barcodes_variants_mapq0_tbl[[sample_id]], q0)
        barcodes_variants_mapq20_tbl[[sample_id]] <- bind_rows(barcodes_variants_mapq20_tbl[[sample_id]], q20)
        rm(q0)
        rm(q20)
      }
    }

    barcodes_variants_mapq0_tbl[[sample_id]] <- barcodes_variants_mapq0_tbl[[sample_id]] %>%
      left_join(patient_sample_names_tbl, by = "sample")
    barcodes_variants_mapq20_tbl[[sample_id]] <- barcodes_variants_mapq20_tbl[[sample_id]] %>%
      left_join(patient_sample_names_tbl, by = "sample")

  }

  # somatic phasing variants

  phasing_variants_mapq0_tbl <- NULL
  phasing_variants_mapq20_tbl <- NULL

  for (sample_id in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    for (chromosome in chromosome_tbl$chromosome) {
      if (is.null(phasing_variants_mapq0_tbl)) {
        phasing_variants_mapq0_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.phasing_variants.tsv", sep = ".")),
          col_types = "cciccciccddddiiiiiddddiiiii") %>%
          mutate(sample = sample_id)
        phasing_variants_mapq20_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.phasing_variants.tsv", sep = ".")),
          col_types = "cciccciccddddiiiiiddddiiiii") %>%
          mutate(sample = sample_id)
      } else {
        q0 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.phasing_variants.tsv", sep = ".")),
          col_types = "cciccciccddddiiiiiddddiiiii") %>%
          mutate(sample = sample_id)
        q20 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.phasing_variants.tsv", sep = ".")),
          col_types = "cciccciccddddiiiiiddddiiiii") %>%
          mutate(sample = sample_id)
        phasing_variants_mapq0_tbl <- bind_rows(phasing_variants_mapq0_tbl, q0)
        phasing_variants_mapq20_tbl <- bind_rows(phasing_variants_mapq20_tbl, q20)
        rm(q0)
        rm(q20)
      }
    }
  }

  phasing_variants_mapq0_tbl <- phasing_variants_mapq0_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")
  phasing_variants_mapq20_tbl <- phasing_variants_mapq20_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

  # somatic somatic per phase set

  somatic_per_phase_set_mapq0_tbl <- NULL
  somatic_per_phase_set_mapq20_tbl <- NULL

  for (sample_id in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    for (chromosome in chromosome_tbl$chromosome) {
      if (is.null(somatic_per_phase_set_mapq0_tbl)) {
        somatic_per_phase_set_mapq0_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.somatic_per_phase_set.tsv", sep = ".")),
          col_types = "cciiiiiiiiii") %>%
          mutate(sample = sample_id)
        somatic_per_phase_set_mapq20_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.somatic_per_phase_set.tsv", sep = ".")),
          col_types = "cciiiiiiiiii") %>%
          mutate(sample = sample_id)
      } else {
        q0 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.somatic_per_phase_set.tsv", sep = ".")),
          col_types = "cciiiiiiiiii") %>%
          mutate(sample = sample_id)
        q20 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.somatic_per_phase_set.tsv", sep = ".")),
          col_types = "cciiiiiiiiii") %>%
          mutate(sample = sample_id)
        somatic_per_phase_set_mapq0_tbl <- bind_rows(somatic_per_phase_set_mapq0_tbl, q0)
        somatic_per_phase_set_mapq20_tbl <- bind_rows(somatic_per_phase_set_mapq20_tbl, q20)
        rm(q0)
        rm(q20)
      }
    }
  }

  somatic_per_phase_set_mapq0_tbl <- somatic_per_phase_set_mapq0_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")
  somatic_per_phase_set_mapq20_tbl <- somatic_per_phase_set_mapq20_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

  # somatic somatic per phase set

  variant_pairs_mapq0_tbl <- NULL
  variant_pairs_mapq20_tbl <- NULL

  for (sample_id in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    for (chromosome in chromosome_tbl$chromosome) {
      if (is.null(variant_pairs_mapq0_tbl)) {
        variant_pairs_mapq0_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.variant_pairs.tsv", sep = ".")),
          col_types = "ccccccciiii") %>%
          mutate(sample = sample_id)
        variant_pairs_mapq20_tbl <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.variant_pairs.tsv", sep = ".")),
          col_types = "ccccccciiii") %>%
          mutate(sample = sample_id)
      } else {
        q0 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq0.variant_pairs.tsv", sep = ".")),
          col_types = "ccccccciiii") %>%
          mutate(sample = sample_id)
        q20 <- read_tsv(str_c(
          str_c("data/somatic/", sample_id, "/"),
          str_c(sample_id, chromosome,
                "mapq20.variant_pairs.tsv", sep = ".")),
          col_types = "ccccccciiii") %>%
          mutate(sample = sample_id)
        variant_pairs_mapq0_tbl <- bind_rows(variant_pairs_mapq0_tbl, q0)
        variant_pairs_mapq20_tbl <- bind_rows(variant_pairs_mapq20_tbl, q20)
        rm(q0)
        rm(q20)
      }
    }
  }

  variant_pairs_mapq0_tbl <- variant_pairs_mapq0_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")
  variant_pairs_mapq20_tbl <- variant_pairs_mapq20_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

  # somatic somatic barcodes (sombx)

  sombx_mapq0_tbl <- NULL
  sombx_mapq20_tbl <- NULL

  for (sample_id in patient_sample_names_tbl %>%
       filter(cnv_maf_status) %>% pull(sample)) {
    if (is.null(sombx_mapq0_tbl)) {
      sombx_mapq0_tbl <- read_tsv(str_c(
        str_c("data/somatic/", sample_id, "/"),
        str_c(sample_id,
              "mapq0.sombx.tsv", sep = ".")),
        col_types = "ccccccccccc") %>%
        mutate(sample = sample_id)
      sombx_mapq20_tbl <- read_tsv(str_c(
        str_c("data/somatic/", sample_id, "/"),
        str_c(sample_id,
              "mapq20.sombx.tsv", sep = ".")),
        col_types = "ccccccccccc") %>%
        mutate(sample = sample_id)
    } else {
      q0 <- read_tsv(str_c(
        str_c("data/somatic/", sample_id, "/"),
        str_c(sample_id,
              "mapq0.sombx.tsv", sep = ".")),
        col_types = "ccccccccccc") %>%
        mutate(sample = sample_id)
      q20 <- read_tsv(str_c(
        str_c("data/somatic/", sample_id, "/"),
        str_c(sample_id,
              "mapq20.sombx.tsv", sep = ".")),
        col_types = "ccccccccccc") %>%
        mutate(sample = sample_id)
      sombx_mapq0_tbl <- bind_rows(sombx_mapq0_tbl, q0)
      sombx_mapq20_tbl <- bind_rows(sombx_mapq20_tbl, q20)
      rm(q0)
      rm(q20)
    }
  }

  sombx_mapq0_tbl <- sombx_mapq0_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")
  sombx_mapq20_tbl <- sombx_mapq20_tbl %>%
    left_join(patient_sample_names_tbl, by = "sample")

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
    left_join(patient_sample_names_tbl, by = "sample")

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
    left_join(patient_sample_names_tbl, by = "sample")

  # IBD segment ancestry

  segment_ancestry_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
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
    left_join(patient_sample_names_tbl, by = "sample")

  # IBD segment overlap

  segment_overlap_tbl <- NULL

  for (sample_id in patient_sample_names_tbl$sample) {
    for (chromosome in chromosome_tbl$chromosome) {
      if ( is.null(segment_overlap_tbl) ) {
        segment_overlap_tbl <- read_tsv(str_c("data/ancestry/",
                                              sample_id,
                                              "/",
                                              str_c(sample_id,
                                                    chromosome,
                                                    "IBD_segment_overlap.tsv",
                                                    sep = ".")),
                                        col_types = "icciciciiddciiidddccccc")  %>%
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
    left_join(patient_sample_names_tbl, by = "sample")

  # clean up

  rm(sample_id)
  rm(sample_id_lr)
  rm(chromosome)

  # save current objects to load next time

  save(list = ls(pattern = "_tbl"), file = str_c(input_data_path_str), envir = .GlobalEnv)
  rm(last_updated, input_data_path_str)

}

sessionInfo()
# Steven Foltz (envest) January 2020
