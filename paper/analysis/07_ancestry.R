################################################################################
# Use ancestry information
################################################################################

main = "figures/07_ancestry/main/"
supp = "figures/07_ancestry/supplementary/"

# Create directories
dir.create(main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)


segment_overlap_tbl %>%
  filter(IBD_or_HBD == "IBD") %>%
  group_by(Sample1_ID, Sample1_haplotype, Sample2_ID, Sample2_haplotype, IBD_start_position, IBD_end_positions) %>% summarize(count = n(), ps = str_c(Phase_set_key, collapse = ", ")) %>% filter(count > 1) %>% View()
