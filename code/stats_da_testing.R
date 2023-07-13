{
  source('code/load_data.R')
  library(phyloseq)
  library(microbiomeMarker)
  library(ANCOMBC)
  library(beepr)
}

tax_table <- taxonomy %>% 
  mutate(taxon = str_replace_all(taxon, '\\.', ' ')) %>% 
  pivot_wider(names_from = 'level', values_from = 'taxon') %>% 
  filter(Kingdom == 'Bacteria') %>% 
  mutate(Family = paste(Kingdom, Phylum, Class, Order, Family, sep = '.')) %>%
  mutate(Phylum = paste(Kingdom, Phylum, sep = '.')) %>%
  column_to_rownames(var = 'featureid') %>% 
  as.matrix()

asv_table <- table %>% 
  column_to_rownames(var = 'featureid') %>% 
  as.matrix()

sample_metadata <- metadata %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         gm = ifelse(gm == 1, 'GM-Low', 'GM-High')) %>% 
  column_to_rownames(var='sampleid') %>% 
  as.data.frame()

ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(tax_table)
METADATA <- sample_data(sample_metadata) 

physeq <- phyloseq(ASV, TAX, METADATA)
physeq_room_temp <- subset_samples(physeq, study == 'room_temp')
gm_1_physeq_room_temp <- subset_samples(physeq_room_temp, gm == 'GM-Low')

gm_4_physeq_room_temp <- subset_samples(physeq_room_temp, gm == 'GM-High')
gm_4_physeq_room_temp <- subset_samples(gm_4_physeq_room_temp, cage != '143')

sample_data(gm_4_physeq_room_temp)

gm1_family_out = ancombc2(
  data = gm_1_physeq_room_temp, 
  fix_formula = 'time_point',
  rand_formula = '(time_point | cage)',
  tax_level = 'Family',
  p_adj_method = "BH", 
  group = "time_point", 
  alpha = 0.05, 
  lib_cut = 0, 
  struc_zero = TRUE, 
  assay_name = 'counts', 
  pseudo_sens = TRUE,
  verbose = TRUE,
  global = T,
  pairwise = T
)

gm1_family_out$res %>% 
  write_tsv('stats/room_temp/ancombc2_gm1_family_res.tsv')
gm1_family_out$zero_ind %>% 
  write_tsv('stats/room_temp/ancombc2_gm1_family_sz.tsv')
gm1_family_out$res_pair %>% 
  write_tsv('stats/room_temp/ancombc2_gm1_family_res_pair.tsv')

 
gm1_phylum_out = ancombc2(
  data = gm_1_physeq_room_temp, 
  fix_formula = 'time_point',
  rand_formula = '(time_point | cage)',
  tax_level = 'Phylum',
  p_adj_method = "BH", 
  group = "time_point", 
  alpha = 0.05, 
  lib_cut = 0, 
  struc_zero = TRUE, 
  assay_name = 'counts', 
  pseudo_sens = TRUE,
  verbose = TRUE,
  global = T,
  pairwise = T
)

gm1_phylum_out$res %>% 
  write_tsv('stats/room_temp/ancombc2_gm1_phylum_res.tsv')
gm1_phylum_out$zero_ind %>% 
  write_tsv('stats/room_temp/ancombc2_gm1_phylum_sz.tsv')
gm1_phylum_out$res_pair %>% 
  write_tsv('stats/room_temp/ancombc2_gm1_phylum_res_pair.tsv')


####

gm4_family_out = ancombc2(
  data = gm_4_physeq_room_temp, 
  fix_formula = 'time_point',
  rand_formula = '(time_point | cage)',
  tax_level = 'Phylum',
  p_adj_method = "BH", 
  group = "time_point", 
  alpha = 0.05, 
  lib_cut = 0, 
  struc_zero = TRUE, 
  assay_name = 'counts', 
  pseudo_sens = TRUE,
  verbose = TRUE,
  global = T,
  pairwise = T
)

gm4_family_out$res %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_family_res.tsv')
gm4_family_out$zero_ind %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_family_sz.tsv')
gm4_family_out$res_pair %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_family_res_pair.tsv')


gm4_phylum_out = ancombc2(
  data = gm_4_physeq_room_temp, 
  fix_formula = 'time_point',
  tax_level = 'Phylum',
  p_adj_method = "BH", 
  group = "time_point", 
  alpha = 0.05, 
  lib_cut = 0, 
  struc_zero = TRUE, 
  assay_name = 'counts', 
  pseudo_sens = TRUE,
  verbose = TRUE,
  global = T,
  pairwise = T
)

gm4_phylum_out$res_global %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_phylum_res_global.tsv')
gm4_phylum_out$res %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_phylum_res.tsv')
gm4_phylum_out$zero_ind %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_phylum_sz.tsv')
gm4_phylum_out$res_pair %>% 
  write_tsv('stats/room_temp/ancombc2_gm4_phylum_res_pair.tsv')


taxa_ids <- taxonomy %>% 
  pivot_wider(names_from = 'level',
              values_from = 'taxon') %>% 
  mutate(Family = paste(Kingdom, Phylum, Class, Order, Family, sep = '.'),
         Phylum = paste(Kingdom, Phylum, sep = '.')) %>% 
  select(featureid, Phylum, Family)

gm1_room_temp_ids <- metadata %>% 
  filter(study == 'room_temp') %>% 
  filter(gm == 1)

gm1_taxa_table <- table %>% 
  select(featureid, any_of(gm1_room_temp_ids$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid') %>% 
  group_by(featureid) %>% 
  left_join(., taxa_ids)

gm1_taxa_table %>% 
  select(featureid, sampleid, value, Family) %>% 
  group_by(sampleid, Family) %>% 
  summarise(count = sum(value), .groups = 'drop') %>% 
  right_join(., gm1_room_temp_ids) %>% 
  group_by(Family) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point) %>% 
  as_tibble() %>% 
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  arrange(p.adj.bh) %>% 
  write_tsv('stats/room_temp/anova_res_gm_low_family.tsv')
gm1_taxa_table %>% 
  select(featureid, sampleid, value, Phylum) %>% 
  group_by(sampleid, Phylum) %>% 
  summarise(count = sum(value), .groups = 'drop') %>% 
  right_join(., gm1_room_temp_ids) %>% 
  group_by(Phylum) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point) %>% 
  as_tibble() %>% 
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  arrange(p.adj.bh) %>% 
  write_tsv('stats/room_temp/anova_res_gm_low_phylum.tsv')
gm4_room_temp_ids <- metadata %>% 
  filter(study == 'room_temp') %>% 
  filter(gm == 4)

gm4_taxa_table <- table %>% 
  select(featureid, any_of(gm4_room_temp_ids$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid') %>% 
  group_by(featureid) %>% 
  left_join(., taxa_ids)

gm4_taxa_table %>% 
  select(featureid, sampleid, value, Family) %>% 
  group_by(sampleid, Family) %>% 
  summarise(count = sum(value), .groups = 'drop') %>% 
  right_join(., gm4_room_temp_ids) %>% 
  group_by(Family) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point) %>% 
  as_tibble() %>% 
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  arrange(p.adj.bh) %>% 
  write_tsv('stats/room_temp/anova_res_gm_high_family.tsv')
gm4_taxa_table %>% 
  select(featureid, sampleid, value, Phylum) %>% 
  group_by(sampleid, Phylum) %>% 
  summarise(count = sum(value), .groups = 'drop') %>% 
  right_join(., gm4_room_temp_ids) %>% 
  group_by(Phylum) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point) %>% 
  as_tibble() %>% 
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  arrange(p.adj.bh) %>% 
  write_tsv('stats/room_temp/anova_res_gm_high_phylum.tsv')
#### ---


f_names <- taxonomy %>% 
  pivot_wider(names_from = 'level',
              values_from = 'taxon') %>% 
  mutate(Family = paste(Kingdom, Phylum, Class, Order, Family, sep = '.')) %>% 
  pivot_longer(-featureid, names_to = 'level', values_to = 'taxon') %>% 
  filter(level == 'Family')
p_names <- taxonomy %>% 
  pivot_wider(names_from = 'level',
              values_from = 'taxon') %>% 
  mutate(Phylum = paste(Kingdom, Phylum, sep = '.')) %>% 
  pivot_longer(-featureid, names_to = 'level', values_to = 'taxon') %>% 
  filter(level == 'Phylum')

f_table <- table %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  left_join(., f_names) %>% 
  group_by(sampleid, taxon) %>% 
  summarize(count = sum(count)) %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  left_join(., metadata) 
  
p_table <- table %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  left_join(., p_names) %>% 
  group_by(sampleid, taxon) %>% 
  summarize(count = sum(count)) %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  left_join(., metadata) 

f_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 1) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point * sample_type) %>% 
  as_tibble() %>%
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  write_tsv('stats/colon_position/colon_position_gm1_family_anova_res.tsv')

f_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 1) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::tukey_hsd(count ~ time_point * sample_type) %>% 
  as_tibble() %>% 
  write_tsv('stats/colon_position/colon_position_gm1_family_tukey_res.tsv')

f_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 4) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point * sample_type) %>% 
  as_tibble() %>%
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  write_tsv('stats/colon_position/colon_position_gm4_family_anova_res.tsv')

f_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 4) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::tukey_hsd(count ~ time_point * sample_type) %>% 
  as_tibble() %>% 
  write_tsv('stats/colon_position/colon_position_gm4_family_tukey_res.tsv')

p_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 1) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point * sample_type) %>% 
  as_tibble() %>%
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  write_tsv('stats/colon_position/colon_position_gm1_phylum_anova_res.tsv')

p_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 1) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::tukey_hsd(count ~ time_point * sample_type) %>% 
  as_tibble() %>% 
  write_tsv('stats/colon_position/colon_position_gm1_phylum_tukey_res.tsv')

p_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 4) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::anova_test(count ~ time_point * sample_type) %>% 
  as_tibble() %>%
  mutate(p.adj.bh = p.adjust(p, method = 'BH')) %>% 
  write_tsv('stats/colon_position/colon_position_gm4_phylum_anova_res.tsv')

p_table %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 4) %>%
  group_by(taxon) %>% 
  filter(sum(count) > 0) %>% 
  rstatix::tukey_hsd(count ~ time_point * sample_type) %>% 
  as_tibble() %>% 
  write_tsv('stats/colon_position/colon_position_gm4_phylum_tukey_res.tsv')

  

