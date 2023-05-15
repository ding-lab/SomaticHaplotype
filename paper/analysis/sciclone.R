cnv_sciclone_list <- cnv_tbl %>%
  mutate(segment_mean = 2*(2^log2.copyRatio)) %>%
  select(chr = chrom,
         start,
         stop = end,
         segment_mean) %>%
  filter(!is.na(chr)) %>%
  split(f = cnv_tbl$sample[!is.na(cnv_tbl$chrom)])

vaf_sciclone_list <- sombx_mapq20_tbl %>%
  rowwise() %>%
  mutate(pos = as.numeric(position),
         ref_reads = length(str_split(ref_barcodes_H1,
                                      pattern = ";",
                                      simplify = TRUE)) +
           length(str_split(ref_barcodes_H2,
                            pattern = ";",
                            simplify = TRUE)) +
           length(str_split(ref_barcodes_None,
                            pattern = ";",
                            simplify = TRUE)),
         var_reads = length(str_split(alt_barcodes_H1,
                                      pattern = ";",
                                      simplify = TRUE)) +
           length(str_split(alt_barcodes_H2,
                            pattern = ";",
                            simplify = TRUE)) +
           length(str_split(alt_barcodes_None,
                            pattern = ";",
                            simplify = TRUE)),
         vaf = var_reads/(ref_reads + var_reads)) %>%
  select(chr = chromosome,
         pos,
         ref_reads,
         var_reads,
         vaf) %>%
  ungroup() %>%
  split(f = sombx_mapq20_tbl$sample)

sciclone_results_list <- purrr::map2(vaf_sciclone_list,
                                     cnv_sciclone_list,
                                     function(x,y) sciClone::sciClone(vafs = as.data.frame(x),
                                                                      copyNumberCalls = as.data.frame(y),
                                                                      sampleNames = "c1",
                                                                      useSexChrs = FALSE))

cluster_vaf_per_sample <- purrr::map2(vaf_sciclone_list,
                                      sciclone_results_list,
                                      function(x,y) data.frame(x %>% filter(!is.na(chr)),
                                                               cluster = y@vafs.merged$cluster) %>%
                                        filter(!is.na(cluster)) %>%
                                        group_by(cluster) %>%
                                        summarize(n_var = n(),
                                                  mean_vaf = mean(vaf)))

tumor_purity_estimates_list <-  cluster_vaf_per_sample %>%
  purrr::map(function(x) 2*max(x$mean_vaf))
