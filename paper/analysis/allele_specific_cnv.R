overlaps_TF <- function(start1, end1, start2, end2){
  if (start2 <= start1 & start1 <= end2 & end2 <= end1 |
      start1 <= start2 & end2 <= end1 |
      start1 <= start2 & start2 <= end1 & end2 >= end1 |
      start2 <= start1 & end2 >= end1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

get_ps_most_overlap <- function(my_phase_sets_tbl, my_sample, my_chrom, my_start, my_end){

  x <- my_phase_sets_tbl %>%
    filter(sample == my_sample) %>%
    filter(chr == my_chrom) %>%
    filter(length_variants >= 1000,
           n_variants_total >= 100) %>%
    rowwise() %>%
    filter(overlaps_TF(first_variant_pos, last_variant_pos, my_start, my_end))

  if (nrow(x) == 0) {
    return(NA)
  } else {
    x %>%
      rowwise() %>%
      mutate(overlap = min(last_variant_pos, my_end) - max(first_variant_pos, my_start)) %>%
      arrange(desc(overlap)) %>%
      head(1) %>%
      mutate(ps_most_overlap = str_c(c(ps_id, first_variant_pos, last_variant_pos), collapse = ",")) %>%
      pull(ps_most_overlap) %>% return()
  }

}

cnvs_overlapping_phase_sets <- cnv_tbl %>%
  filter(!is.na(chrom),
         end - start >= 1000,
         log2.copyRatio <= -0.25 | log2.copyRatio >= 0.2) %>%
  rowwise() %>%
  mutate(ps_most_overlap = get_ps_most_overlap(phase_sets_tbl, sample, chrom, start, end)) %>%
  filter(!is.na(ps_most_overlap)) %>%
  separate(ps_most_overlap, into = c("ps", "ps_start", "ps_end"), sep = ",") %>%
  mutate(ps_start = as.numeric(ps_start),
         ps_end = as.numeric(ps_end)) %>%
  rowwise() %>%
  mutate(overlap_start = max(start, ps_start),
         overlap_end = min(end, ps_end)) %>%
  filter(overlap_end - overlap_start >= 1000) %>%
  select(chrom, start, end, sample, log2.copyRatio,
         ps, ps_start, ps_end, overlap_start, overlap_end)

write_tsv(cnvs_overlapping_phase_sets,
          "data/cnvs_overlapping_phase_sets.tsv")
