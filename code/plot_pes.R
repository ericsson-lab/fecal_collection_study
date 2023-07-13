source('code/data_alpha_stats.R')
library(rstatix)

## running stats for indicated values 
## Pull partial eta squared value to plot at end
chao1_pes <- alpha_stats %>% 
  filter(study == 'room_temp') %>% 
  anova_test(chao1 ~ gm * time_point, effect.size = 'pes') %>% 
  as_tibble() %>% 
  rename(effect = "Effect",
         chao1 = 'pes') %>% 
  select(effect, chao1)

shannon_pes <- alpha_stats %>% 
  filter(study == 'room_temp') %>% 
  anova_test(shannon ~ gm * time_point, effect.size = 'pes') %>% 
  as_tibble() %>% 
  rename(effect = "Effect",
         shannon = 'pes') %>% 
  select(effect, shannon)

source('code/data_beta_stats.R')

bc_dist <- generate_dist("bray")

room_temp_samples <- metadata %>% 
  filter(study == 'room_temp')

room_temp_bc_dist <- usedist::dist_subset(bc_dist, room_temp_samples$sampleid )

bc_res <- vegan::adonis2(room_temp_bc_dist ~ gm * time_point, 
                         data = room_temp_samples, permutations = 9999)  

# calculating partial eta squared manually
bc_res_table <- bc_res %>% 
  as_tibble(rownames = 'effect')
ss_error <- bc_res_table$SumOfSqs[4]
bc_permanova_ges <- bc_res_table[1:3,] %>% 
  mutate(pes = SumOfSqs/(sum(SumOfSqs) + ss_error)) %>% 
  select(effect, pes) %>% 
  rename(bc_permanova = 'pes')

bc_permanova_ges$effect <- factor(bc_permanova_ges$effect, levels = rev(c('gm',  'time_point',
                                                                          'gm:time_point')))
fig6a_time <- bc_permanova_ges %>% 
  left_join(., chao1_pes) %>% 
  left_join(., shannon_pes) %>% 
  filter(effect %in% c('gm', 'time_point')) %>% 
  mutate(effect = case_when(effect == 'gm' ~ 'GM',
                            effect == 'time_point' ~ 'Time Point')) %>% 
  pivot_longer(-effect) %>% 
  
  ggplot(aes(x = effect, y = value, shape = name, color = name)) +
  geom_point(size = 5) +
  stat_summary(aes(group = effect), geom = 'errorbar', fun.min = mean, fun.max = mean,
               color = 'red',
               show.legend = F, size = 1.5, width = 0.3) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0,1), 
                     expand = expansion(mult = c(0, 0.1))) +
  scale_shape_manual(values = c(15,16,17),
                     labels = c('PERMANOVA', 'Chao1', 'Shannon')) +
  scale_color_manual(values = c('gray', '#808080', 'black'), labels = c('PERMANOVA', 'Chao1', 'Shannon')) +
  coord_flip() +
  labs(y = expression(bold('Partial Eta Squared (\u03B7'[p]^2*')'))) +
  ggprism::theme_prism() +
  theme(
    axis.title.y = element_blank(), aspect.ratio = 3/2,
    legend.text = element_text(face = 'bold', color = 'black', size = 14)
  )


bc_permanova_ges %>% 
  left_join(., chao1_pes) %>% 
  left_join(., shannon_pes) %>% 
  mutate(avg = (bc_permanova + chao1 +shannon)/3) %>% 
  mutate(bc_permanova = round(bc_permanova, 3),
         chao1 = round(chao1, 3),
         shannon = round(shannon, 3)) %>% 
  clipr::write_clip()


## For position study
chao1_pos_pes <- alpha_stats %>% 
  filter(study == 'colon_position') %>% 
  anova_test(chao1 ~ gm * time_point * sample_type, effect.size = 'pes') %>% 
  as_tibble() %>% 
  rename(effect = "Effect",
         chao1 = 'pes') %>% 
  select(effect, chao1)

shannon_pos_pes <- alpha_stats %>% 
  filter(study == 'colon_position') %>% 
  anova_test(shannon ~ gm * time_point * sample_type, effect.size = 'pes') %>% 
  as_tibble() %>% 
  rename(effect = "Effect",
         shannon = 'pes') %>% 
  select(effect, shannon)

bc_dist <- generate_dist("bray")
bc_dist_sample <- attributes(bc_dist)$Labels

pos_samples <- metadata %>% 
  filter(study == 'colon_position') %>% 
  filter(sampleid %in% bc_dist_sample)

pos_bc_dist <- usedist::dist_subset(bc_dist, pos_samples$sampleid )

pos_bc_res <- vegan::adonis2(pos_bc_dist ~ gm * time_point * sample_type, 
                             data = pos_samples, permutations = 9999, method = 'bray')  

pos_bc_res_table <- pos_bc_res %>% 
  as_tibble(rownames = 'effect')

pos_ss_error <- pos_bc_res_table$SumOfSqs[8]

pos_bc_permanova_ges <- pos_bc_res_table[1:7,] %>% 
  mutate(pes = SumOfSqs/(sum(SumOfSqs) + ss_error)) %>% 
  select(effect, pes) %>% 
  rename(bc_permanova = 'pes')

pos_bc_permanova_ges$effect <- factor(pos_bc_permanova_ges$effect, levels = rev(c('gm',  'time_point', 'sample_type',
                                                                                  'gm:time_point', 'gm:sample_type',
                                                                                  'time_point:sample_type',
                                                                                  'gm:time_point:sample_type')))
fig6b_pos <- pos_bc_permanova_ges %>% 
  left_join(., chao1_pos_pes) %>% 
  left_join(., shannon_pos_pes) %>% 
  filter(effect %in% c('gm', 'time_point', 'sample_type')) %>% 
  mutate(effect = case_when(effect == 'gm' ~ 'GM',
                            effect == 'time_point' ~ 'Time Point',
                            effect == 'sample_type' ~ 'Sample Type')) %>% 
  pivot_longer(-effect) %>% 
  
  ggplot(aes(x = effect, y = value, shape = name, color = name)) +
  geom_point(size = 5) +
  stat_summary(aes(group = effect), geom = 'errorbar', fun.min = mean, fun.max = mean,
               color = 'red',
               show.legend = F, size = 1.5, width = 0.3) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0,1), 
                     expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = c('gray', '#808080', 'black'), labels = c('PERMANOVA', 'Chao1', 'Shannon')) +
  scale_shape_manual(values = c(15,16,17),
                     labels = c('PERMANOVA', 'Chao1', 'Shannon')) +
  coord_flip() +
  labs(y = expression(bold('Partial Eta Squared (\u03B7'[p]^2*')'))) +
  ggprism::theme_prism()  +
  theme(
    
    axis.title.y = element_blank(), aspect.ratio = 3/2,
    legend.text = element_text(face = 'bold', color = 'black', size = 14),
  )

pos_bc_permanova_ges %>% 
  left_join(., chao1_pos_pes) %>% 
  left_join(., shannon_pos_pes) %>% 
  mutate(avg = round((bc_permanova + chao1 +shannon)/3,3)) %>% 
  mutate(bc_permanova = round(bc_permanova, 3),
         chao1 = round(chao1, 3),
         shannon = round(shannon, 3)) %>% 
  clipr::write_clip()

