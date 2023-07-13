source('code/load_data.R')
library(RColorBrewer)
## Need to make
# - Tree
# - Phylum level classification
# - 1 Cecum
# - 2 Proximal
# - 3 Mid
# - 4 Distal
# - 5 Differ across sample type

mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))

### GM - High

gm4_ids <- metadata %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 4)

gm4_featureids <- table %>% 
  select(featureid, any_of(gm4_ids$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid') %>% 
  group_by(featureid) %>% 
  filter(sum(value) > 0) %>% 
  pull(featureid) %>% 
  unique()

gm4_collapsed_taxonomy <- taxonomy %>% 
  filter(featureid %in% gm4_featureids) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'level', values_from = 'taxon') %>% 
  mutate(Family = paste(Kingdom, Phylum, Class, Order, Family, sep = '.'),
         Phylum = paste(Kingdom, Phylum, sep = '.')) %>% 
  select(featureid, Phylum, Family)

gm4_enriched <- table %>% 
  select(featureid, any_of(gm4_ids$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid') %>% 
  group_by(featureid) %>% 
  filter(sum(value) > 0) %>% 
  ungroup() %>% 
  left_join(., gm4_collapsed_taxonomy) %>% 
  group_by(sampleid, Phylum) %>% 
  summarize(count = sum(value), .groups = 'drop') %>% 
  left_join(., metadata) %>% 
  group_by(Phylum, sample_type, time_point) %>% 
  summarise(avg_count = mean(count)) %>% 
  pivot_wider(names_from = 'time_point', values_from = 'avg_count') %>%
  rename(sample_1 = 'sample_type') %>% 
  mutate(greater = ifelse(`4p` > `7a`, 'PM', 'AM')) 

gm4_taxa_names <- gm4_collapsed_taxonomy %>% 
  select(Phylum, Family) %>% 
  distinct() %>% 
  filter(Phylum != 'ssigned.ssigned')

gm4_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/tree.txt', col_names = F)

gm4_phylum_length <- gm4_taxa_names %>% 
  pull(Phylum) %>% 
  unique() %>% 
  length()

pal <- as_tibble(mycolors(gm4_phylum_length))

gm4_pal <- gm4_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  mutate(clade_attribute = 'clade_marker_color') %>% 
  cbind(., pal) 

gm4_pal %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring_inner_clade_color.txt',
            col_names = F)

gm4_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  mutate(clade_attribute = 'clade_marker_size',
         clade_marker_size = 50) %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring_inner_clade_marker_size.txt',
            col_names = F)
gm4_phylum_anova_res <- read_tsv('stats/colon_position/colon_position_gm4_phylum_anova_res.tsv')
gm4_phylum_tukey_res <- read_tsv('stats/colon_position/colon_position_gm4_phylum_tukey_res.tsv')

gm4_phylum_anova_res %>% 
  filter(Effect == 'sample_type') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 5,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring5_differ_across_sample_location.txt', 
            col_names = F)
gm4_phylum_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 6,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring6_differ_across_time_point_shape.txt', 
            col_names = F)
gm4_phylum_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = 6,
         ring_color = '#BEBDB8') %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring6_differ_across_time_point_color.txt', 
            col_names = F)
gm4_phylum_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_height',
         ring_level = 6,
         ring_height = .6) %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring6_differ_across_time_point_height.txt', 
            col_names = F)

gm4_phylum_to_plot <- gm4_phylum_tukey_res %>% 
  filter(term == 'time_point:sample_type') %>% 
  separate(group1, into = c('time_1', 'sample_1'),remove = F) %>% 
  separate(group2, into = c('time_2', 'sample_2'), remove = F) %>% 
  filter(sample_1 == sample_2) %>% 
  filter(p.adj < 0.05) %>% 
  select(taxon, sample_1) %>% 
  mutate(enriched_key = paste(taxon, sample_1, sep = '_'))

am_pm_key <- gm4_enriched %>% 
  mutate(enriched_key = paste(Phylum, sample_1, sep = '_')) %>% 
  filter(enriched_key %in% gm4_phylum_to_plot$enriched_key) %>% 
  ungroup() %>% 
  select(enriched_key, greater)

gm4_phylum_to_plot %>% 
  left_join(., am_pm_key) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = case_when(sample_1 == 'cecum' ~ 1,
                                sample_1 == 'proximal' ~ 2,
                                sample_1 == 'mid' ~ 3,
                                sample_1 == 'distal' ~ 4),
         ring_color = case_when(greater == 'AM' ~ '#0000FF',
                                greater == 'PM' ~ '#9696FF')) %>% 
  select(taxon, ring_atribute, ring_level, ring_color) %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/phylum/inputs/ring1_4_tukey_res_gm_low_family.txt',
            col_names = F)



## GM-Low

gm1_ids <- metadata %>% 
  filter(study == 'colon_position') %>% 
  filter(gm == 1)

gm1_featureids <- table %>% 
  select(featureid, any_of(gm1_ids$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid') %>% 
  group_by(featureid) %>% 
  filter(sum(value) > 0) %>% 
  pull(featureid) %>% 
  unique()

gm1_collapsed_taxonomy <- taxonomy %>% 
  filter(featureid %in% gm1_featureids) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'level', values_from = 'taxon') %>% 
  mutate(Family = paste(Kingdom, Phylum, Class, Order, Family, sep = '.'),
         Phylum = paste(Kingdom, Phylum, sep = '.')) %>% 
  select(featureid, Phylum, Family)

gm1_enriched <- table %>% 
  select(featureid, any_of(gm1_ids$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid') %>% 
  group_by(featureid) %>% 
  filter(sum(value) > 0) %>% 
  ungroup() %>% 
  left_join(., gm1_collapsed_taxonomy) %>% 
  group_by(sampleid, Phylum) %>% 
  summarize(count = sum(value), .groups = 'drop') %>% 
  left_join(., metadata) %>% 
  group_by(Phylum, sample_type, time_point) %>% 
  summarise(avg_count = mean(count)) %>% 
  pivot_wider(names_from = 'time_point', values_from = 'avg_count') %>%
  rename(sample_1 = 'sample_type') %>% 
  mutate(greater = ifelse(`4p` > `7a`, 'PM', 'AM')) 

gm1_taxa_names <- gm1_collapsed_taxonomy %>% 
  select(Phylum, Family) %>% 
  distinct()

gm1_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/tree.txt', col_names = F)

gm1_phylum_names <- gm1_taxa_names %>%
  pull(Phylum) %>%
  unique() 

gm4_pal %>% 
  filter(Phylum %in% gm1_phylum_names ) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring_inner_clade_color.txt',
            col_names = F)

gm1_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  mutate(clade_attribute = 'clade_marker_size',
         clade_marker_size = 50) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring_inner_clade_marker_size.txt',
            col_names = F)
gm1_phylum_anova_res <- read_tsv('stats/colon_position/colon_position_gm1_phylum_anova_res.tsv')
gm1_phylum_tukey_res <- read_tsv('stats/colon_position/colon_position_gm1_phylum_tukey_res.tsv')

gm1_phylum_anova_res %>% 
  filter(Effect == 'sample_type') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 5,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring5_differ_across_sample_location.txt', 
            col_names = F)
gm1_phylum_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 6,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring6_differ_across_time_point_shape.txt', 
            col_names = F)
gm1_phylum_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = 6,
         ring_color = '#BEBDB8') %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring6_differ_across_time_point_color.txt', 
            col_names = F)
gm1_phylum_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_height',
         ring_level = 6,
         ring_height = .6) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring6_differ_across_time_point_height.txt', 
            col_names = F)

gm1_phylum_to_plot <- gm1_phylum_tukey_res %>% 
  filter(term == 'time_point:sample_type') %>% 
  separate(group1, into = c('time_1', 'sample_1'),remove = F) %>% 
  separate(group2, into = c('time_2', 'sample_2'), remove = F) %>% 
  filter(sample_1 == sample_2) %>% 
  # filter(p.adj < 0.05) %>% 
  filter(taxon == 'Bacteria.Actinobacteriota') # need to pull all sample locations so graphlan can build tree. only ring4 is sig differ. will replace colors of other tissue with white
  # select(taxon, sample_1) %>% 
  # mutate(enriched_key = paste(taxon, sample_1, sep = '_'))

# am_pm_key <- gm1_enriched %>% 
#   mutate(enriched_key = paste(Phylum, sample_1, sep = '_')) %>% 
#   filter(enriched_key %in% gm1_phylum_to_plot$enriched_key) %>% 
#   ungroup() %>% 
#   select(enriched_key, greater) 
#   Will manually code in.

gm1_phylum_to_plot %>% 
  # left_join(., am_pm_key) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = case_when(sample_1 == 'cecum' ~ 1,
                                sample_1 == 'proximal' ~ 2,
                                sample_1 == 'mid' ~ 3,
                                sample_1 == 'distal' ~ 4),
         ring_color = case_when(sample_1 == 'distal' ~ '#FF0000',
                                TRUE ~ '#FFFFFF')) %>% 
  select(taxon, ring_atribute, ring_level, ring_color) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/phylum/inputs/ring1_4_tukey_res_gm_low_phylum.txt',
            col_names = F)

#USE LEGENDS FROM FAMILY PLOTS
