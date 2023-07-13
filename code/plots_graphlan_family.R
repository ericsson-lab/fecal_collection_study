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
  group_by(sampleid, Family) %>% 
  summarize(count = sum(value), .groups = 'drop') %>% 
  left_join(., metadata) %>% 
  group_by(Family, sample_type, time_point) %>% 
  summarise(avg_count = mean(count)) %>% 
  pivot_wider(names_from = 'time_point', values_from = 'avg_count') %>%
  rename(sample_1 = 'sample_type') %>% 
  mutate(greater = ifelse(`4p` > `7a`, 'PM', 'AM')) 

gm4_taxa_names <- gm4_collapsed_taxonomy %>% 
  select(Phylum, Family) %>% 
  distinct() %>% 
  filter(Phylum != 'ssigned.ssigned')

gm4_taxa_names %>% 
  select(Family) %>% 
  distinct() %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/tree.txt', col_names = F)

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
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring_inner_clade_color.txt',
            col_names = F)

gm4_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  mutate(clade_attribute = 'clade_marker_size',
         clade_marker_size = 30) %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring_inner_clade_marker_size.txt',
            col_names = F)
gm4_family_anova_res <- read_tsv('stats/colon_position/colon_position_gm4_family_anova_res.tsv')
gm4_family_tukey_res <- read_tsv('stats/colon_position/colon_position_gm4_family_tukey_res.tsv')

gm4_family_anova_res %>% 
  filter(Effect == 'sample_type') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 5,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring5_differ_across_sample_location.txt', 
            col_names = F)
gm4_family_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 6,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring6_differ_across_time_point_shape.txt', 
            col_names = F)
gm4_family_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = 6,
         ring_color = '#BEBDB8') %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring6_differ_across_time_point_color.txt', 
            col_names = F)
gm4_family_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_height',
         ring_level = 6,
         ring_height = .6) %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring6_differ_across_time_point_height.txt', 
            col_names = F)

gm4_family_to_plot <- gm4_family_tukey_res %>% 
  filter(term == 'time_point:sample_type') %>% 
  separate(group1, into = c('time_1', 'sample_1'),remove = F) %>% 
  separate(group2, into = c('time_2', 'sample_2'), remove = F) %>% 
  filter(sample_1 == sample_2) %>% 
  filter(p.adj < 0.05) %>% 
  select(taxon, sample_1) %>% 
  mutate(enriched_key = paste(taxon, sample_1, sep = '_'))

am_pm_key <- gm4_enriched %>% 
  mutate(enriched_key = paste(Family, sample_1, sep = '_')) %>% 
  filter(enriched_key %in% gm4_family_to_plot$enriched_key) %>% 
  ungroup() %>% 
  select(enriched_key, greater)

gm4_family_to_plot %>% 
  left_join(., am_pm_key) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = case_when(sample_1 == 'cecum' ~ 1,
                                sample_1 == 'proximal' ~ 2,
                                sample_1 == 'mid' ~ 3,
                                sample_1 == 'distal' ~ 4),
         ring_color = case_when(greater == 'AM' ~ '#0000FF',
                                greater == 'PM' ~ '#9696FF')) %>% 
  select(taxon, ring_atribute, ring_level, ring_color) %>% 
  write_tsv('plots/colon_position/graphlan/gm_high/family/inputs/ring1_4_tukey_res_gm_low_family.txt',
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
  group_by(sampleid, Family) %>% 
  summarize(count = sum(value), .groups = 'drop') %>% 
  left_join(., metadata) %>% 
  group_by(Family, sample_type, time_point) %>% 
  summarise(avg_count = mean(count)) %>% 
  pivot_wider(names_from = 'time_point', values_from = 'avg_count') %>%
  rename(sample_1 = 'sample_type') %>% 
  mutate(greater = ifelse(`4p` > `7a`, 'PM', 'AM')) 

gm1_taxa_names <- gm1_collapsed_taxonomy %>% 
  select(Phylum, Family) %>% 
  distinct()

gm1_taxa_names %>% 
  select(Family) %>% 
  distinct() %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/tree.txt', col_names = F)

gm1_phylum_names <- gm1_taxa_names %>%
  pull(Phylum) %>%
  unique() 
# 
# pal <- as_tibble(mycolors(gm1_phylum_length))
gm4_pal %>% 
  filter(Phylum %in% gm1_phylum_names ) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring_inner_clade_color.txt',
            col_names = F)
  
gm1_taxa_names %>% 
  select(Phylum) %>% 
  distinct() %>% 
  mutate(clade_attribute = 'clade_marker_size',
         clade_marker_size = 30) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring_inner_clade_marker_size.txt',
            col_names = F)
gm1_family_anova_res <- read_tsv('stats/colon_position/colon_position_gm1_family_anova_res.tsv')
gm1_family_tukey_res <- read_tsv('stats/colon_position/colon_position_gm1_family_tukey_res.tsv')

gm1_family_anova_res %>% 
  filter(Effect == 'sample_type') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 5,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring5_differ_across_sample_location.txt', 
            col_names = F)

gm1_family_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_shape',
         ring_level = 6,
         ring_shape = 'v') %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring6_differ_across_time_point_shape.txt', 
            col_names = F)
gm1_family_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = 6,
         ring_color = '#BEBDB8') %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring6_differ_across_time_point_color.txt', 
            col_names = F)
gm1_family_anova_res %>% 
  filter(Effect == 'time_point') %>% 
  filter(p.adj.bh < 0.05) %>% 
  select(taxon) %>% 
  mutate(ring_atribute = 'ring_height',
         ring_level = 6,
         ring_height = .6) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring6_differ_across_time_point_height.txt', 
            col_names = F)

gm1_family_to_plot <- gm1_family_tukey_res %>% 
  filter(term == 'time_point:sample_type') %>% 
  separate(group1, into = c('time_1', 'sample_1'),remove = F) %>% 
  separate(group2, into = c('time_2', 'sample_2'), remove = F) %>% 
  filter(sample_1 == sample_2) %>% 
  filter(p.adj < 0.05) %>% 
  select(taxon, sample_1) %>% 
  mutate(enriched_key = paste(taxon, sample_1, sep = '_'))

am_pm_key <- gm1_enriched %>% 
  mutate(enriched_key = paste(Family, sample_1, sep = '_')) %>% 
  filter(enriched_key %in% gm1_family_to_plot$enriched_key) %>% 
  ungroup() %>% 
  select(enriched_key, greater)

gm1_family_to_plot %>% 
  left_join(., am_pm_key) %>% 
  mutate(ring_atribute = 'ring_color',
         ring_level = case_when(sample_1 == 'cecum' ~ 1,
                                sample_1 == 'proximal' ~ 2,
                                sample_1 == 'mid' ~ 3,
                                sample_1 == 'distal' ~ 4),
         ring_color = case_when(greater == 'AM' ~ '#FF0000',
                                greater == 'PM' ~ '#FFa5a5')) %>% 
  select(taxon, ring_atribute, ring_level, ring_color) %>% 
  write_tsv('plots/colon_position/graphlan/gm_low/family/inputs/ring1_4_tukey_res_gm_low_family.txt',
            col_names = F)

order <- gm4_pal %>% 
  separate(Phylum, into = c('Kingdom', 'Phylum')) %>% 
  select(Phylum, value)
gm4_pal %>% 
  separate(Phylum, into = c('Kingdom', 'Phylum')) %>% 
  select(Phylum, value) %>% 
  mutate(Phylum = case_when(Phylum == 'Bacteria'~ 'Unresolved Bacteria',
                            TRUE ~ Phylum)) %>% 
  ggplot(aes(x = 0, y = factor(Phylum, levels = c(rev(`Phylum`))))) +
  geom_tile(aes(fill = value), color = 'black', size = 1) +
  scale_fill_identity() +
  scale_y_discrete(position = "right") +
  theme_void() +
  theme(
    axis.text.y = element_text(face = 'bold', hjust = 0, size = 12,
                               margin = margin(0,0,0,5)),
    aspect.ratio = 9,
  )
ggsave('plots/colon_position/graphlan/phylum_key.png', width = 3, height =3)


enriched <- tibble(time = c('AM', 'PM'),
                   low = c('#FF0000', '#FFa5a5'),
                   high = c('#0000FF', '#9696FF' ))

enriched %>% 
  pivot_longer(-time) %>% 
  ggplot(aes(x = factor(name, levels = c('low', 'high')), y = factor(time, levels = c('PM', 'AM')))) +
  geom_tile(aes(fill = value), color = 'black', size = 1.5) +
  scale_fill_identity() +
  theme_void() +
  scale_x_discrete(position = 'top', 
                   labels = c(
                     expression(bold('GM'['Low'])),
                     expression(bold('GM'['High']))
                   )) +
  scale_y_discrete(position = 'right') +
  theme(
    axis.text.x = element_text(face = 'bold', size = 21, hjust = 0.5, margin = margin(0,0,1,0)),
    axis.text.y = element_text(face = 'bold', size = 21, hjust = 0, margin = margin(0,0,0,1))
  )
ggsave('plots/colon_position/graphlan/time_key.png', width = 2.5, height =2.5)



