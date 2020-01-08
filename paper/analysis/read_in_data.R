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
          print(chromosome)
        }
      }
    }
  }
}
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

# Steven Foltz (envest) January 2020
