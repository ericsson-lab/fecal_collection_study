{
  source('code/load_data.R')
  
  library(RColorBrewer)
  library(ggVennDiagram)
}

# First identify sex-dependent differces in taxonomy
t0_samples <- metadata %>% 
  filter(time_point == '7a' & 
           sample_type %in% c('distal', 'feces'))

genus_taxonomy <- taxonomy %>% 
  filter(level == 'Genus')

# Started filtering by prevalance, set > 0%
prevalance_list <- table %>% 
  select(featureid, any_of(t0_samples$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  group_by(featureid) %>% 
  filter(sum(count) > 0) %>% 
  left_join(., genus_taxonomy) %>% 
  group_by(sampleid, taxon) %>% 
  summarize(count = sum(count), .groups = 'drop') %>% 
  group_by(sampleid, taxon) %>% 
  mutate(present = ifelse(count > 0, 1, 0)) %>% 
  group_by(taxon) %>% 
  left_join(., metadata) %>% 
  group_by(gm, taxon) %>% 
  summarize(prevalance = sum(present)/n()) %>% 
  filter(prevalance > 0)

# Taxa in each GM
gm1_taxa <- prevalance_list %>% filter(gm == 1) %>% pull(taxon)
gm4_taxa <- prevalance_list %>% filter(gm == 4) %>% pull(taxon)

# Plot Venn diagram for supplemental figure
plot_sex_venn <- function(list, palette_name) {
  taxa_list <- table %>% 
    select(featureid, any_of(t0_samples$sampleid)) %>% 
    pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
    group_by(featureid) %>% 
    filter(sum(count) > 0) %>% 
    left_join(., genus_taxonomy) %>% 
    filter(taxon %in% {{list}}) %>% 
    left_join(., metadata) %>% 
    group_by(sex, taxon) %>% 
    summarize(present = sum(count) > 0) 
  list_m <- taxa_list %>% 
    filter(sex == 'M' & present == TRUE) %>% 
    pull(taxon)
  list_f <- taxa_list %>% 
    filter(sex == 'F' & present == TRUE) %>% 
    pull(taxon)

  taxa_venn <- list(
    Male = list_m,
    Female = list_f
  )
  
  ggVennDiagram(taxa_venn, label = 'both', label_percent_digit = 2, set_color = 'black', label_color = 'black',
                label_size = 6, set_size = 8, edge_size = 3) +
    scale_fill_distiller(palette = palette_name, direction = 1) +
    scale_color_manual(values = c('black', 'black')) +
    theme(
      legend.position = 'none'
    )
}

# Make supplemental figures 
figS1d_gm1_sex_venn <- plot_sex_venn(list = gm1_taxa, palette_name = 'Reds') 
figS1e_gm4_sex_venn <- plot_sex_venn(list = gm4_taxa, palette_name = 'Blues') 

# calculate abundance
abundance <- table %>% 
  select(featureid, any_of(t0_samples$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  group_by(featureid) %>% 
  filter(sum(count) > 0) %>% 
  left_join(., genus_taxonomy) %>% 
  group_by(sampleid, taxon) %>% 
  summarize(count = sum(count), .groups = 'drop') %>%  
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  left_join(., metadata) %>% 
  group_by(gm, taxon) %>% 
  summarize(abundance = mean(rel_abund)) 

# GM low taxa
taxa_list_1 <- table %>% 
  select(featureid, any_of(t0_samples$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  group_by(featureid) %>% 
  filter(sum(count) > 0) %>% 
  left_join(., genus_taxonomy) %>% 
  filter(taxon %in% gm1_taxa) %>% 
  left_join(., metadata) %>% 
  group_by(sex, taxon) %>% 
  summarize(present = sum(count) > 0) 
list_m_1 <- taxa_list_1 %>% 
  filter(sex == 'M' & present == TRUE) %>% 
  pull(taxon)
list_f_1 <- taxa_list_1 %>% 
  filter(sex == 'F' & present == TRUE) %>% 
  pull(taxon)

gm1_shared <- intersect(list_m_1, list_f_1)
gm1_m <- setdiff(list_m_1, list_f_1)
gm1_f <- setdiff(list_f_1, list_m_1)

# GM high taxa
taxa_list_4 <- table %>% 
  select(featureid, any_of(t0_samples$sampleid)) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  group_by(featureid) %>% 
  filter(sum(count) > 0) %>% 
  left_join(., genus_taxonomy) %>% 
  filter(taxon %in% gm4_taxa) %>% 
  left_join(., metadata) %>% 
  group_by(sex, taxon) %>% 
  summarize(present = sum(count) > 0) 
list_m_4 <- taxa_list_4 %>% 
  filter(sex == 'M' & present == TRUE) %>% 
  pull(taxon)
list_f_4 <- taxa_list_4 %>% 
  filter(sex == 'F' & present == TRUE) %>% 
  pull(taxon)

gm4_shared <- intersect(list_m_4, list_f_4)
gm4_m <- setdiff(list_m_4, list_f_4)
gm4_f <- setdiff(list_f_4, list_m_4)

# Make prevalance/abundance figure
figS1f_prevelance_abundance <- prevalance_list %>% 
  left_join(., abundance) %>% 
  mutate(sex = case_when(gm == 1 & taxon %in% gm1_shared ~ 'Core',
                         gm == 1 & taxon %in% gm1_m ~ 'Male',
                         gm == 1 & taxon %in% gm1_f ~ 'Female',
                         gm == 4 & taxon %in% gm4_shared ~ 'Core',
                         gm == 4 & taxon %in% gm4_m ~ 'Male',
                         gm == 4 & taxon %in% gm4_f ~ 'Female')) %>%
  ggplot(aes(x = abundance, y = prevalance)) +
  geom_point(aes(color = factor(gm), 
                 shape = factor(sex, levels = c('Core',
                                                'Female', 
                                                'Male'))), 
             size = 4, alpha = 0.5) +
  scale_color_manual(values = c('red', 'blue'),
                     labels = c(expression(bold('GM'['Low'])),
                                expression(bold('GM'['High']))),
                     name = 'GM') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05)),
                     labels = scales::percent,
                     breaks = seq(0,1,0.25)) +
  scale_x_continuous(expand = expansion(mult=c(0.05, 0.1)),
                     labels = scales::percent) +
  scale_shape_manual(values = c(0, 17, 16), name = 'Sex') +
  ggprism::theme_prism() +
  labs(x = 'Mean Abundance', y = 'Prevalence') +
  theme(
    aspect.ratio = 2/2,
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.text = element_text(hjust = 0, size = 14, face = 'bold'),
    legend.title = element_text(hjust = 0.5, size = 14)
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1)),
    shape = guide_legend(override.aes = list(alpha = 1))
  )

prevalance_list %>% 
  left_join(., abundance) %>% 
  mutate(sex = case_when(gm == 1 & taxon %in% gm1_shared ~ 'Core',
                         gm == 1 & taxon %in% gm1_m ~ 'Male',
                         gm == 1 & taxon %in% gm1_f ~ 'Female',
                         gm == 4 & taxon %in% gm4_shared ~ 'Core',
                         gm == 4 & taxon %in% gm4_m ~ 'Male',
                         gm == 4 & taxon %in% gm4_f ~ 'Female')) %>% 
  filter(sex == 'Core') %>% 
  summarise(
    mean_prev = round(mean(prevalance), 4),
    median_prev = round(median(prevalance), 4),
    mean_abund = round(mean(abundance), 4),
    median_abudn = round(median(abundance), 4) ) %>% 
  clipr::write_clip()


#####
## Room temp study
#####
time_metadata <- metadata %>% 
  filter(study == 'room_temp')
# pull family level names
f_names <- taxonomy %>% 
  filter(level == 'Family') 

mean_rel_abund_table <- table %>% 
  select(featureid, time_metadata$sampleid) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  left_join(., f_names) %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  group_by(sampleid, taxon) %>% 
  summarize(rel_abund = sum(rel_abund), .groups = 'drop') %>% 
  group_by(sampleid) %>% 
  left_join(., time_metadata) %>% 
  group_by(taxon, time_point, gm) %>% 
  summarise(mean_rel_abund = mean(rel_abund)) %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         gm = ifelse(gm == 1, 'GM-Low', 'GM-High'))

# plotting only abudance > 1%
# < 1% grouped into 'Other'
taxon_pool <- mean_rel_abund_table %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean_rel_abund) < 0.01, .groups = 'drop') 

data_to_plot <- inner_join(mean_rel_abund_table, taxon_pool, by = "taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(taxon, time_point, gm) %>% 
  summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
  arrange(desc(mean_rel_abund)) %>% 
  mutate(taxon = case_when(taxon == 'uncultured' ~ 'Uncultured', 
                           TRUE ~ taxon)) %>% 
  mutate(gm = case_when(gm == 'GM-Low' ~ 'bold(GM[Low])',
                        gm == 'GM-High' ~ 'bold(GM[High])')) 

taxa_list <- unique(data_to_plot$taxon)
taxa_list <- taxa_list[!taxa_list %in% c("Other")]

data_to_plot$taxon <- factor(data_to_plot$taxon,
                             levels = c(taxa_list, "Other"))

spectral_palette <- colorRampPalette(brewer.pal(11, "Spectral"))

num_of_colors <- taxa_list %>% 
  length()

fig3_family <-data_to_plot %>% 
  group_by(gm, time_point) %>% 
  ggplot(aes(x = time_point, y = mean_rel_abund, fill = taxon)) +
  geom_area() +
  ggprism::theme_prism() +
  scale_x_continuous(expand = c(0,0), breaks = c(0,1,2,3,6,9)) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c(spectral_palette(num_of_colors),"#606060"), name = 'Family') +
  labs(y = 'Mean Relative\nAbundance (%)',
       x = 'Hours Post Collection') +
  guides(fill=guide_legend(ncol=2)) +
  facet_wrap(factor(gm, levels = c('bold(GM[Low])', 'bold(GM[High])'))~., nrow = 1, scales = 'free', labeller = label_parsed) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 12),
    legend.title = element_text(face = 'bold', color = 'black', size = 12, hjust = 0.5),
    aspect.ratio = 3/2,
    strip.text = element_text(size = 14)
  ) 

# Pull phylum names
p_names <- taxonomy %>% 
  filter(level == 'Phylum') 

mean_rel_abund_table <- table %>% 
  select(featureid, time_metadata$sampleid) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  left_join(., p_names) %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  group_by(sampleid, taxon) %>% 
  summarize(rel_abund = sum(rel_abund), .groups = 'drop') %>% 
  group_by(sampleid) %>% 
  left_join(., time_metadata) %>% 
  group_by(taxon, time_point, gm) %>% 
  summarise(mean_rel_abund = mean(rel_abund)) %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         gm = ifelse(gm == 1, 'GM-Low', 'GM-High'))

# Relative abundance > 0.1%
taxon_pool <- mean_rel_abund_table %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean_rel_abund) < 0.001, .groups = 'drop') 

data_to_plot <- inner_join(mean_rel_abund_table, taxon_pool, by = "taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(taxon, time_point, gm) %>% 
  summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
  arrange(desc(mean_rel_abund)) %>% 
  mutate(taxon = case_when(taxon == 'Firmicutes' ~ 'Bacillota',
                           taxon == 'Proteobacteria' ~ 'Pseudomonadota',
                           TRUE ~ taxon),
         gm = case_when(gm == 'GM-Low' ~ 'bold(GM[Low])',
                        gm == 'GM-High' ~ 'bold(GM[High])')) 

taxa_list <- unique(data_to_plot$taxon)
taxa_list <- taxa_list[!taxa_list %in% c("Other")]

data_to_plot$taxon <- factor(data_to_plot$taxon,
                             levels = c(taxa_list, "Other"))

spectral_palette <- colorRampPalette(brewer.pal(11, "Spectral"))

num_of_colors <- taxa_list %>% 
  length()

fig3_phylum <- data_to_plot %>% 
  group_by(gm, time_point) %>% 
  ggplot(aes(x = time_point, y = mean_rel_abund, fill = taxon)) +
  geom_area() +
  ggprism::theme_prism() +
  scale_x_continuous(expand = c(0,0), breaks = c(0,1,2,3,6,9)) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c(spectral_palette(num_of_colors),"#606060"), name = 'Phylum') +
  labs(y = 'Mean Relative\nAbundance (%)',
       x = 'Hours Post Collection') +
  facet_wrap(factor(gm, levels = c('bold(GM[Low])', 'bold(GM[High])'))~., nrow = 1, scales = 'free', labeller = label_parsed) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 12),
    legend.title = element_text(face = 'bold', color = 'black', size = 12, hjust = 0.5),
    aspect.ratio = 3/2,
    strip.text = element_text(size = 14)
  ) 


#######
## Study position
#######
# position_metadata <- metadata %>% 
#   filter(study == 'colon_position')
# 
# # family level plot
# f_names <- taxonomy %>% 
#   filter(level == 'Family') 
# 
# mean_rel_abund_table <- table %>% 
#   select(featureid, any_of(position_metadata$sampleid)) %>%
#   pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
#   left_join(., f_names) %>% 
#   group_by(sampleid) %>% 
#   mutate(rel_abund = count/sum(count)) %>% 
#   group_by(sampleid, taxon) %>% 
#   summarize(rel_abund = sum(rel_abund), .groups = 'drop') %>% 
#   group_by(sampleid) %>% 
#   left_join(., position_metadata) %>% 
#   group_by(taxon, gm, sample_type, time_point) %>% 
#   summarise(mean_rel_abund = mean(rel_abund)) %>% 
#   mutate(time_point = case_when(time_point == '7a' ~ 'AM',
#                                 time_point == '4p' ~ 'PM'),
#          gm = ifelse(gm == 1, 'GM-Low', 'GM-High'),
#          sample_type = str_to_sentence(sample_type),
#          gm = case_when(gm == 'GM-Low' ~ 'bold(GM[Low])',
#                         gm == 'GM-High' ~ 'bold(GM[High])'))
# 
# taxon_pool <- mean_rel_abund_table %>% 
#   group_by(taxon) %>% 
#   summarise(pool = max(mean_rel_abund) < 0.01, .groups = 'drop') 
# 
# data_to_plot <- inner_join(mean_rel_abund_table, taxon_pool, by = "taxon") %>% 
#   mutate(taxon = if_else(pool, "Other", taxon)) %>% 
#   group_by(taxon, gm, sample_type, time_point) %>% 
#   summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
#   arrange(desc(mean_rel_abund)) %>% 
#   mutate(taxon = case_when(taxon == 'uncultured' ~ 'Uncultured', 
#                            TRUE ~ taxon)) 
# 
# taxa_list <- unique(data_to_plot$taxon)
# taxa_list <- taxa_list[!taxa_list %in% c("Other")]
# 
# data_to_plot$taxon <- factor(data_to_plot$taxon,
#                              levels = c(taxa_list, "Other"))
# 
# spectral_palette <- colorRampPalette(brewer.pal(11, "Spectral"))
# 
# num_of_colors <- taxa_list %>% 
#   length()
# 
# data_to_plot %>% 
#   mutate(time_point = paste0('bold(',time_point, ')')) %>% 
#   ggplot(aes(x = factor(sample_type, 
#                         levels = c('Cecum', 'Proximal', 'Mid', 'Distal')), 
#              y = mean_rel_abund, fill = taxon, group = taxon)) +
#   geom_area() +
#   ggprism::theme_prism() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_fill_manual(values = c(spectral_palette(num_of_colors),"#606060"), name = 'Taxon') +
#   labs(y = 'Mean Relative Abundance (%)',
#        x = 'Hours Post Collection') +
#   theme(
#     axis.text.y = element_text(size = 10),
#     axis.text.x = element_text(size = 10, ),
#     axis.title = element_text(size = 14),
#     axis.title.x = element_blank(),
#     strip.text = element_text(size = 14),
#     legend.text = element_text(face = 'bold', color = 'black', size = 10),
#     legend.title = element_text(face = 'bold', color = 'black', size = 10, hjust = 0),
#     panel.spacing.x = unit(5, "mm"),
#     aspect.ratio = 2.75/2,
#   ) +
#   ggh4x::facet_grid2(time_point ~ factor(gm, levels = c('bold(GM[Low])', 'bold(GM[High])'))  ,labeller = label_parsed, scales = "free", independent = "all")
# 
# ggsave('plots/colon_position/f_level_taxonomy.png', 
#        width = 8, height = 6, bg = 'white')  
# 
# 
# 
# 
# 
# #
# #
# 
# source('code/load_data.R')
# library(RColorBrewer)
# 
# position_metadata <- metadata %>% 
#   filter(study == 'colon_position')
# 
# f_names <- taxonomy %>% 
#   filter(level == 'Family') 
# 
# mean_rel_abund_table <- table %>% 
#   select(featureid, any_of(position_metadata$sampleid)) %>%
#   pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
#   left_join(., f_names) %>% 
#   group_by(sampleid) %>% 
#   mutate(rel_abund = count/sum(count)) %>% 
#   group_by(sampleid, taxon) %>% 
#   summarize(rel_abund = sum(rel_abund), .groups = 'drop') %>% 
#   group_by(sampleid) %>% 
#   left_join(., position_metadata) %>% 
#   filter(cage != 172143) %>% 
#   group_by(taxon, gm, sample_type, time_point) %>% 
#   summarise(mean_rel_abund = mean(rel_abund)) %>% 
#   mutate(time_point = case_when(time_point == '7a' ~ 'AM',
#                                 time_point == '4p' ~ 'PM'),
#          gm = ifelse(gm == 1, 'GM-Low', 'GM-High'),
#          sample_type = str_to_sentence(sample_type))
# 
# 
# taxon_pool <- mean_rel_abund_table %>% 
#   group_by(taxon) %>% 
#   summarise(pool = max(mean_rel_abund) < 0.01, .groups = 'drop') 
# 
# data_to_plot <- inner_join(mean_rel_abund_table, taxon_pool, by = "taxon") %>% 
#   mutate(taxon = if_else(pool, "Other", taxon)) %>% 
#   group_by(taxon, gm, sample_type, time_point) %>% 
#   summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
#   arrange(desc(mean_rel_abund)) %>% 
#   mutate(taxon = case_when(taxon == 'uncultured' ~ 'Uncultured', 
#                            TRUE ~ taxon)) 
# 
# taxa_list <- unique(data_to_plot$taxon)
# taxa_list <- taxa_list[!taxa_list %in% c("Other")]
# 
# data_to_plot$taxon <- factor(data_to_plot$taxon,
#                              levels = c(taxa_list, "Other"))
# 
# spectral_palette <- colorRampPalette(brewer.pal(11, "Spectral"))
# 
# num_of_colors <- taxa_list %>% 
#   length()
# 
# data_to_plot %>% 
#   ggplot(aes(x = factor(sample_type, 
#                         levels = c('Cecum', 'Proximal', 'Mid', 'Distal')), 
#              y = mean_rel_abund, fill = taxon, group = taxon)) +
#   geom_area() +
#   ggprism::theme_prism() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_fill_manual(values = c(spectral_palette(num_of_colors),"#606060"), name = 'Taxon') +
#   labs(y = 'Mean Relative Abundance (%)',
#        x = 'Hours Post Collection') +
#   theme(
#     axis.text.y = element_text(size = 10),
#     axis.text.x = element_text(size = 10, ),
#     axis.title = element_text(size = 10),
#     axis.title.x = element_blank(),
#     strip.text = element_text(size = 10),
#     legend.text = element_text(face = 'bold', color = 'black', size = 10),
#     legend.title = element_text(face = 'bold', color = 'black', size = 12, hjust = 0),
#     panel.spacing.x = unit(5, "mm")
#     # aspect.ratio = 3/2
#   ) +
#   ggh4x::facet_grid2(gm ~ time_point , scales = "free", independent = "all")
# 
# 
# ggsave('plots/colon_position/f_level_taxonomy.png', 
#        width = 8, height = 5, bg = 'white')  
# 
#  
# 
# p_names <- taxonomy %>% 
#   filter(level == 'Phylum') 
# 
# mean_rel_abund_table <- table %>% 
#   select(featureid, any_of(position_metadata$sampleid)) %>%
#   pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
#   left_join(., p_names) %>% 
#   group_by(sampleid) %>% 
#   mutate(rel_abund = count/sum(count)) %>% 
#   group_by(sampleid, taxon) %>% 
#   summarize(rel_abund = sum(rel_abund), .groups = 'drop') %>% 
#   group_by(sampleid) %>% 
#   left_join(., position_metadata) %>% 
#   filter(cage != 172143) %>% 
#   group_by(taxon, gm, sample_type, time_point) %>% 
#   summarise(mean_rel_abund = mean(rel_abund)) %>% 
#   mutate(time_point = case_when(time_point == '7a' ~ 'AM',
#                                 time_point == '4p' ~ 'PM'),
#          gm = ifelse(gm == 1, 'GM-Low', 'GM-High'),
#          sample_type = str_to_sentence(sample_type),
#          gm = case_when(gm == 'GM-Low' ~ 'bold(GM[Low])',
#                         gm == 'GM-High' ~ 'bold(GM[High])'))
# 
# 
# taxon_pool <- mean_rel_abund_table %>% 
#   group_by(taxon) %>% 
#   summarise(pool = max(mean_rel_abund) < 0.01, .groups = 'drop') 
# 
# data_to_plot <- inner_join(mean_rel_abund_table, taxon_pool, by = "taxon") %>% 
#   mutate(taxon = if_else(pool, "Other", taxon)) %>% 
#   group_by(taxon, gm, sample_type, time_point) %>% 
#   summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
#   arrange(desc(mean_rel_abund)) %>% 
#   mutate(taxon = case_when(taxon == 'uncultured' ~ 'Uncultured', 
#                            taxon == 'Firmicutes' ~ 'Bacillota',
#                            taxon == 'Proteobacteria' ~ 'Pseudomonadota',
#                            TRUE ~ taxon)) 
# 
# taxa_list <- unique(data_to_plot$taxon)
# taxa_list <- taxa_list[!taxa_list %in% c("Other")]
# 
# data_to_plot$taxon <- factor(data_to_plot$taxon,
#                              levels = c(taxa_list, "Other"))
# 
# spectral_palette <- colorRampPalette(brewer.pal(11, "Spectral"))
# 
# num_of_colors <- taxa_list %>% 
#   length()
# 
# data_to_plot %>% 
#   mutate(time_point = paste0('bold(',time_point, ')')) %>% 
#   ggplot(aes(x = factor(sample_type, 
#                         levels = c('Cecum', 'Proximal', 'Mid', 'Distal')), 
#              y = mean_rel_abund, fill = taxon, group = taxon)) +
#   geom_area() +
#   ggprism::theme_prism() +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_fill_manual(values = c(spectral_palette(num_of_colors),"#606060"), name = 'Taxon') +
#   labs(y = 'Mean Relative Abundance (%)',
#        x = 'Hours Post Collection') +
#   theme(
#     axis.text.y = element_text(size = 10),
#     axis.text.x = element_text(size = 10, ),
#     axis.title = element_text(size = 14),
#     axis.title.x = element_blank(),
#     strip.text = element_text(size = 14),
#     legend.text = element_text(face = 'bold', color = 'black', size = 14),
#     legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
#     panel.spacing.x = unit(5, "mm"),
#     aspect.ratio = 2.75/2,
#   ) +
#   ggh4x::facet_grid2(time_point ~ factor(gm, levels = c('bold(GM[Low])', 'bold(GM[High])'))  ,labeller = label_parsed, scales = "free", independent = "all")
# 
# ggsave('plots/colon_position/p_level_taxonomy.png', 
#        width = 8, height = 6, bg = 'white')  

