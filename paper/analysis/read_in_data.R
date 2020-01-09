################################################################################
# Read in all data for analysis
################################################################################

# shh! (this is a library)
library(tidyverse)

# patient and sample names

patient_sample_names <- read_tsv("data/patient_sample_names.tsv",
                                 col_types = "ccccll")

# chromosomes

chromosome_vector <- str_c("chr", seq(1,22))

# 10x_longranger_summaries

lr_summary_tbl <- NULL

for (sample_id in patient_sample_names$sample) {
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

lr_summary_tbl <- bind_cols(patient_sample_names,
                            lr_summary_tbl)

# phase set summary

phase_set_summary_tbl <- NULL

for (sample_id in patient_sample_names$sample) {
  for (chromosome in chromosome_vector) {
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
  left_join(patient_sample_names, by = "sample")

# phase sets

phase_sets_tbl <- NULL

for (sample_id in patient_sample_names$sample) {
  for (chromosome in chromosome_vector) {
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
  left_join(patient_sample_names, by = "sample")

# extend stats

extend_stats_tbl <- NULL

for (sample_id in patient_sample_names$sample) {
  for (extended_by in patient_sample_names$sample) {
    for (chromosome in chromosome_vector) {
      file_path1 <- str_c("data/extend/", sample_id, "/")
      file_path2 <- str_c(sample_id, chromosome, "extended_by_", sep = ".")
      file_path3 <- str_c(extended_by, "extend_stats.tsv", sep = ".")
      file_path <- str_c(file_path1, file_path2, file_path3)
      if (file.exists(file_path)) {
        if ( is.null(extend_stats_tbl) ) {
          extend_stats_tbl <- read_tsv(file_path, col_types = "ccddccddddddddddcd")  %>%
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

rm(filepath1)
rm(filepath2)
rm(filepath3)
rm(filepath)

extend_stats_tbl <- extend_stats_tbl %>%
  left_join(patient_sample_names, by = "sample")

# extend phase sets

extend_phase_sets_tbl <- NULL

for (sample_id in patient_sample_names$sample) {
  for (extended_by in patient_sample_names$sample) {
    for (chromosome in chromosome_vector) {
      file_path1 <- str_c("data/extend/", sample_id, "/")
      file_path2 <- str_c(sample_id, chromosome, "extended_by_", sep = ".")
      file_path3 <- str_c(extended_by, "extended_phase_sets.tsv", sep = ".")
      file_path <- str_c(file_path1, file_path2, file_path3)
      if (file.exists(file_path)) {
        if ( is.null(extend_phase_sets_tbl) ) {
          extend_phase_sets_tbl <- read_tsv(file_path, col_types = "cciiiiiiiiiicl")  %>%
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

rm(extended_by)

extend_phase_sets_tbl <- extend_phase_sets_tbl %>%
  left_join(patient_sample_names, by = "sample")

# coverage

coverage_tbl <- NULL

for (sample_id in patient_sample_names$sample) {
  for (chromosome in chromosome_vector) {
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
  left_join(patient_sample_names, by = "sample")

# somatic barcodes variants

barcodes_variants_mapq0_tbl <- list() # use a list to keep samples separate
barcodes_variants_mapq20_tbl <- list() # because combined files are too big

for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  barcodes_variants_mapq0_tbl[patient_id] <- NULL
  barcodes_variants_mapq20_tbl[patient_id] <- NULL
  for (chromosome in chromosome_vector) {
    if (is.null(barcodes_variants_mapq0_tbl[patient_id])) {
      barcodes_variants_mapq0_tbl[patient_id] <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.barcodes_variants.tsv", sep = ".")),
        col_types = "ccciccccccc") %>%
        mutate(sample = sample_id)
      barcodes_variants_mapq20_tbl[patient_id] <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.barcodes_variants.tsv", sep = ".")),
        col_types = "ccciccccccc") %>%
        mutate(sample = sample_id)
    } else {
      q0 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.barcodes_variants.tsv", sep = ".")),
        col_types = "ccciccccccc") %>%
        mutate(sample = sample_id)
      q20 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.barcodes_variants.tsv", sep = ".")),
        col_types = "ccciccccccc") %>%
        mutate(sample = sample_id)
      barcodes_variants_mapq0_tbl[patient_id] <- bind_rows(barcodes_variants_mapq0_tbl[patient_id], q0)
      barcodes_variants_mapq20_tbl[patient_id] <- bind_rows(barcodes_variants_mapq0_tbl[patient_id], q20)
      rm(q0)
      rm(q20)
    }
  }

  barcodes_variants_mapq0_tbl[patient_id] %>%
    left_join(patient_sample_names, by = "sample")
  barcodes_variants_mapq20_tbl[patient_id] %>%
    left_join(patient_sample_names, by = "sample")

}

# somatic phasing variants

phasing_variants_mapq0_tbl <- NULL
phasing_variants_mapq20_tbl <- NULL

for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  for (chromosome in chromosome_vector) {
    if (is.null(phasing_variants_mapq0_tbl)) {
      phasing_variants_mapq0_tbl <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.phasing_variants.tsv", sep = ".")),
        col_types = "cciccciccddddiiiiiddddiiiii") %>%
        mutate(sample = sample_id)
      phasing_variants_mapq20_tbl <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.phasing_variants.tsv", sep = ".")),
        col_types = "cciccciccddddiiiiiddddiiiii") %>%
        mutate(sample = sample_id)
    } else {
      q0 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.phasing_variants.tsv", sep = ".")),
        col_types = "cciccciccddddiiiiiddddiiiii") %>%
        mutate(sample = sample_id)
      q20 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.phasing_variants.tsv", sep = ".")),
        col_types = "cciccciccddddiiiiiddddiiiii") %>%
        mutate(sample = sample_id)
      phasing_variants_mapq0_tbl <- bind_rows(phasing_variants_mapq0_tbl, q0)
      phasing_variants_mapq20_tbl <- bind_rows(phasing_variants_mapq0_tbl, q20)
      rm(q0)
      rm(q20)
    }
  }
}

phasing_variants_mapq0_tbl %>%
  left_join(patient_sample_names, by = "sample")
phasing_variants_mapq20_tbl %>%
  left_join(patient_sample_names, by = "sample")

# somatic somatic per phase set

somatic_per_phase_set_mapq0_tbl <- NULL
somatic_per_phase_set_mapq20_tbl <- NULL

for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  for (chromosome in chromosome_vector) {
    if (is.null(somatic_per_phase_set_mapq0_tbl)) {
      somatic_per_phase_set_mapq0_tbl <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.somatic_per_phase_set.tsv", sep = ".")),
        col_types = "cciiiiiiiiii") %>%
        mutate(sample = sample_id)
      somatic_per_phase_set_mapq20_tbl <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.somatic_per_phase_set.tsv", sep = ".")),
        col_types = "cciiiiiiiiii") %>%
        mutate(sample = sample_id)
    } else {
      q0 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.somatic_per_phase_set.tsv", sep = ".")),
        col_types = "cciiiiiiiiii") %>%
        mutate(sample = sample_id)
      q20 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.somatic_per_phase_set.tsv", sep = ".")),
        col_types = "cciiiiiiiiii") %>%
        mutate(sample = sample_id)
      somatic_per_phase_set_mapq0_tbl <- bind_rows(somatic_per_phase_set_mapq0_tbl, q0)
      somatic_per_phase_set_mapq20_tbl <- bind_rows(somatic_per_phase_set_mapq0_tbl, q20)
      rm(q0)
      rm(q20)
    }
  }
}

somatic_per_phase_set_mapq0_tbl %>%
  left_join(patient_sample_names, by = "sample")
somatic_per_phase_set_mapq20_tbl %>%
  left_join(patient_sample_names, by = "sample")

# somatic somatic per phase set

variant_pairs_mapq0_tbl <- NULL
variant_pairs_mapq20_tbl <- NULL

for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  for (chromosome in chromosome_vector) {
    if (is.null(variant_pairs_mapq0_tbl)) {
      variant_pairs_mapq0_tbl <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.variant_pairs.tsv", sep = ".")),
        col_types = "ccccccciiii") %>%
        mutate(sample = sample_id)
      variant_pairs_mapq20_tbl <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.variant_pairs.tsv", sep = ".")),
        col_types = "ccccccciiii") %>%
        mutate(sample = sample_id)
    } else {
      q0 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq0.variant_pairs.tsv", sep = ".")),
        col_types = "ccccccciiii") %>%
        mutate(sample = sample_id)
      q20 <- read_tsv(str_c(
        str_c("data/somatic/", patient_id, "/"),
        str_c(patient_id, chromosome,
              "mapq20.variant_pairs.tsv", sep = ".")),
        col_types = "ccccccciiii") %>%
        mutate(sample = sample_id)
      variant_pairs_mapq0_tbl <- bind_rows(variant_pairs_mapq0_tbl, q0)
      variant_pairs_mapq20_tbl <- bind_rows(variant_pairs_mapq0_tbl, q20)
      rm(q0)
      rm(q20)
    }
  }
}

variant_pairs_mapq0_tbl %>%
  left_join(patient_sample_names, by = "sample")
variant_pairs_mapq20_tbl %>%
  left_join(patient_sample_names, by = "sample")

# somatic somatic barcodes (sombx)

sombx_mapq0_tbl <- NULL
sombx_mapq20_tbl <- NULL

for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  if (is.null(sombx_mapq0_tbl)) {
    sombx_mapq0_tbl <- read_tsv(str_c(
      str_c("data/somatic/", patient_id, "/"),
      str_c(patient_id,
            "mapq0.sombx.tsv", sep = ".")),
      col_types = "ccccccccccc") %>%
      mutate(sample = sample_id)
    sombx_mapq20_tbl <- read_tsv(str_c(
      str_c("data/somatic/", patient_id, "/"),
      str_c(patient_id,
            "mapq20.sombx.tsv", sep = ".")),
      col_types = "ccccccccccc") %>%
      mutate(sample = sample_id)
  } else {
    q0 <- read_tsv(str_c(
      str_c("data/somatic/", patient_id, "/"),
      str_c(patient_id,
            "mapq0.sombx.tsv", sep = ".")),
      col_types = "ccccccccccc") %>%
      mutate(sample = sample_id)
    q20 <- read_tsv(str_c(
      str_c("data/somatic/", patient_id, "/"),
      str_c(patient_id,
            "mapq20.sombx.tsv", sep = ".")),
      col_types = "ccccccccccc") %>%
      mutate(sample = sample_id)
    sombx_mapq0_tbl <- bind_rows(sombx_mapq0_tbl, q0)
    sombx_mapq20_tbl <- bind_rows(sombx_mapq0_tbl, q20)
    rm(q0)
    rm(q20)
  }
}

sombx_mapq0_tbl %>%
  left_join(patient_sample_names, by = "sample")
sombx_mapq20_tbl %>%
  left_join(patient_sample_names, by = "sample")

# Sorted WGS CNV

cnv_tbl <- NULL

for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  if ( is.null(cnv_tbl) ) {
    cnv_tbl <- read_csv(str_c("data/sorted_wgs_cnv/", sample_id,
                              ".lambda90.noscale.CNV.bicseq2.sorted.WGS.tsv"),
                        col_types = "ciid") %>%
      mutate(sample = sample_id)
  } else {
    new_cnv <- read_csv(str_c("data/sorted_wgs_cnv/", sample_id,
                              ".lambda90.noscale.CNV.bicseq2.sorted.WGS.tsv"),
                        col_types = "ciid") %>%
      mutate(sample = sample_id)
    cnv_tbl <- bind_rows(cnv_tbl, new_cnv)
    rm(new_cnv)
  }
}

cnv_tbl <- cnv_tbl %>%
  left_join(patient_sample_names, by = "sample")

# Sorted WGS MAF

maf_tbl <- NULL
maf_col_types = "ccccccccccc-iiiccccccccccc--ccc----------------cici--i---------c-cc-----------"
for (patient_id in patient_sample_names %>%
     filter(cnv_maf_status) %>% pull(sample)) {
  if ( is.null(maf_tbl) ) {
    maf_tbl <- read_csv(str_c("data/sorted_wgs_maf/", sample_id,
                              ".withmutect.maf"),
                        col_types = maf_col_types) %>%
      mutate(sample = sample_id)
  } else {
    new_maf <- read_csv(str_c("data/sorted_wgs_maf/", sample_id,
                              ".withmutect.maf"),
                        col_types = maf_col_types) %>%
      mutate(sample = sample_id)
    maf_tbl <- bind_rows(maf_tbl, new_maf)
    rm(new_maf)
  }
}

rm(maf_col_types)

maf_tbl <- maf_tbl %>%
  left_join(patient_sample_names, by = "sample")

# ancestry 1

# ancestry 2

rm(patient_id)
rm(sample_id)



# Steven Foltz (envest) January 2020
