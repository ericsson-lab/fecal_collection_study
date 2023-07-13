{
  source('code/data_beta_stats.R')
  
  library(usedist)
}
## Figure 1C PCoa
bc_dist <- generate_dist("bray")
# Pull T0 IDs
t0_samples <- metadata %>% 
  filter(time_point == '7a' & 
           sample_type %in% c('distal', 'feces')) 

t0_bc_dist <- dist_subset(bc_dist, t0_samples$sampleid )
t0_room_temp_pcoa <- generate_pcoa(t0_bc_dist)

t0_room_temp_pcoa_data <- t0_room_temp_pcoa[[1]]
t0_room_temp_pcoa_p_var <- t0_room_temp_pcoa[[2]]

# plot 1C PCoA
fig1c_t0_pcoa <- t0_room_temp_pcoa_data %>%
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         order = case_when(time_point == '7a' ~ 0,
                           time_point == '8a' ~ 1,
                           time_point == '9a' ~ 2,
                           time_point == '10a' ~ 3,
                           time_point == '1p' ~ 6,
                           time_point == '4p' ~ 9),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High'),
         sex = case_when(sex == 'M' ~ 'Male',
                         sex == 'F' ~ 'Female')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = factor(gm, levels = c('GM-Low', 'GM-High'))), size = 4) +
  scale_color_manual(values = c("red", 'blue'), name = 'GM', labels = c(expression(bold('GM'['Low']),
                                                                                   bold('GM'['High'])))) +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(t0_room_temp_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(t0_room_temp_pcoa_p_var[2],2)}%')) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    aspect.ratio = 1
  ) +
  guides(
    shape = guide_legend(order = 2),
    color = guide_legend(order = 1)
  )

# Plot Supplementary Figure 1C
figS1C_t0_pcoa_supplemental <- t0_room_temp_pcoa_data %>%
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         order = case_when(time_point == '7a' ~ 0,
                           time_point == '8a' ~ 1,
                           time_point == '9a' ~ 2,
                           time_point == '10a' ~ 3,
                           time_point == '1p' ~ 6,
                           time_point == '4p' ~ 9),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High'),
         sex = case_when(sex == 'M' ~ 'Male',
                         sex == 'F' ~ 'Female')) %>% 
  ggplot(aes(x = PCo1, y = PCo2, shape = sex)) +
  geom_point(aes(color = factor(gm, levels = c('GM-Low', 'GM-High'))), size = 4) +
  scale_color_manual(values = c("red", 'blue'), name = 'GM', labels = c(expression(bold('GM'['Low']),
                                                                                   bold('GM'['High'])))) +
  scale_shape_manual(values = c(16,17), name = 'Sex') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(t0_room_temp_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(t0_room_temp_pcoa_p_var[2],2)}%')) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    aspect.ratio = 1
  ) +
  guides(
    shape = guide_legend(order = 2),
    color = guide_legend(order = 1)
  )

# Run adonis PERMANOVA for main figure
adonis2(t0_bc_dist ~ gm , data = t0_samples, permutations = 9999) %>% 
  as_tibble(rownames = 'effect') %>% 
  select(effect, F, `Pr(>F)`) %>% 
  mutate(F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

# Run adonis PERMANOVA for supplementary figure
adonis2(t0_bc_dist ~ gm * sex, data = t0_samples, permutations = 9999) %>% 
  as_tibble(rownames = 'effect') %>% 
  select(effect, F, `Pr(>F)`) %>% 
  mutate(F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

#######
### Room Temp Data
#######

time_metadata <- metadata %>% 
  filter(study == 'room_temp') %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         order = case_when(time_point == 0 ~ 0,
                           time_point == 1 ~ 1,
                           time_point == 2 ~ 2,
                           time_point == 3 ~ 3,
                           time_point == 6 ~ 4,
                           time_point == 9 ~ 5),
         gm = case_when(gm == 1 ~ 'GM-Low',
                          gm == 4 ~ 'GM-High'))


# Pivot BC_dist long
bc_dist_long <- dist_subset(bc_dist, time_metadata$sampleid ) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid") %>% 
  pivot_longer(-sampleid,
               values_to = "dist",
               names_to = "b") %>% 
  filter(sampleid %in% time_metadata$sampleid  &
           b %in% time_metadata$sampleid )
# Widen dist matrix
room_temp_bc <- bc_dist_long %>% 
  pivot_wider(names_from = 'b', values_from = 'dist') %>% 
  column_to_rownames(var = 'sampleid') %>% 
  as.matrix()
# T0 IDs for room temp study
t0_ids <- time_metadata %>%
  filter(time_point == 0)
# B column needed for dist comparisons
b_metadata <- time_metadata %>% 
  rename(b = 'sampleid') %>% 
  rename(cage_b = 'cage') %>% 
  select(b, cage_b)

## Distance from D0 sample for intracage
fig2d_t0_distance <- bc_dist_long %>% 
  filter(b %in% t0_ids$sampleid) %>% 
  right_join(., time_metadata, by = "sampleid") %>%
  left_join(., b_metadata, by = 'b') %>% 
  filter(cage == cage_b) %>% 
  ggplot(aes(x = time_point, y = dist, color = factor(gm, levels = c('GM-Low', 'GM-High')))) +
  stat_summary(fun = mean,  geom = "line", linewidth = 2.5) +
  stat_summary(fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x),
               geom = "errorbar", linewidth = 1, aes(group = (gm)),
               color = 'black', width = 0.2) +
  scale_color_manual(values = c("red", "blue"), name = 'GM',
                     labels = c(expression(bold('GM'['Low']),
                                           bold('GM'['High'])))) +
  scale_x_continuous(expand = c(0.00,0.05), breaks = c(0,1,2,3,6,9)) +
  scale_y_continuous(expand = expansion( mult = c(0, 0.01)), limits = c(0,1)) +
  scale_linetype_manual(values = c(1,2), name = 'Sex') +
  
  ggprism::theme_prism() +
  labs(x = "Hours Post Collection", y = "Distance From\nInitial Sample") +
  theme(
    axis.text = element_text(face = 'bold'),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14), 
    strip.text = element_text(size = 14),
    legend.position = 'right'
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(order = 2, override.aes = list(linewidth = 2))
    
  )

## Dist from 0 statistics for Fig2D
bc_dist_long %>% 
  filter(b %in% t0_ids$sampleid) %>% 
  right_join(., time_metadata, by = "sampleid") %>%
  left_join(., b_metadata, by = 'b') %>% 
  filter(cage == cage_b) %>% 
  filter(time_point > 0) %>%
  rstatix::anova_test(dist ~ gm * time_point, effect.size = 'pes') %>% 
  as_tibble() %>% 
  clipr::write_clip()

bc_dist_long %>% 
  filter(b %in% t0_ids$sampleid) %>% 
  right_join(., time_metadata, by = "sampleid") %>%
  left_join(., b_metadata, by = 'b') %>% 
  filter(cage == cage_b) %>% 
  filter(time_point > 0) %>%
  rstatix::tukey_hsd(dist ~ gm * as.factor(time_point)) %>% 
  as_tibble() %>% 
  write_tsv('stats/room_temp/dist_from0_tukey_hsd_results.tsv')

## Distance from previous time point
time_a <- time_metadata %>% 
  select(sampleid, cage, order) 
time_b <- time_metadata %>% 
  select(sampleid, time_point, cage, order, sex, gm) %>% 
  rename(b = 'sampleid',
         time_point_b = 'time_point',
         order_b = 'order',
         cage_b = 'cage')

fig2e_prev_t_distance <- bc_dist_long %>% 
  left_join(., time_a, by = 'sampleid') %>% 
  left_join(., time_b, by = 'b') %>% 
  mutate(t_diff = order_b - order) %>% 
  filter(t_diff == 1) %>% 
  filter(cage == cage_b &
           cage != 143) %>% 
  mutate(sex = ifelse(sex == 'M', 'Male', 'Female')) %>% 
  ggplot(aes(x = time_point_b, y = dist, 
             color = factor(gm, levels = c('GM-Low', 'GM-High')), 
             group = factor(gm, levels = c('GM-Low', 'GM-High')))) +
  stat_summary(fun = mean,  geom = "line", linewidth = 2.5) +
  stat_summary(fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x),
               geom = "errorbar", linewidth = 1,
               color = 'black', width = 0.2) +
  scale_color_manual(values = c("red", "blue"), name = 'GM',
                     labels = c(expression(bold('GM'['Low']),
                                           bold('GM'['High'])))) +
  scale_linetype_manual(values = c(1,2), name = 'Sex') +
  scale_x_continuous(expand = c(0.05,0.05), breaks = c(0,1,2,3,6,9)) +
  scale_y_continuous(expand = expansion( mult = c(0, 0.01)), limits = c(0,1)) +
  ggprism::theme_prism() +
  labs(x = "Hours Post Collection", y = "Distance From\nPrevious Sample") +
  theme(
    axis.text = element_text(face = 'bold'),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14), 
    strip.text = element_text(size = 14),
    legend.position = 'right'
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(order = 2, override.aes = list(linewidth = 2))
    
  )

# Test diff from previous time point
bc_dist_long %>% 
  left_join(., time_a, by = 'sampleid') %>% 
  left_join(., time_b, by = 'b') %>% 
  mutate(t_diff = order_b - order) %>% 
  filter(t_diff == 1) %>% 
  filter(cage == cage_b) %>% 
  rstatix::anova_test(dist ~ gm * time_point_b, effect.size = 'pes') %>% 
  as_tibble() %>% 
  select(Effect, F, p, pes) %>% 
  clipr::write_clip()

bc_dist_long %>% 
  left_join(., time_a, by = 'sampleid') %>% 
  left_join(., time_b, by = 'b') %>% 
  mutate(t_diff = order_b - order) %>% 
  filter(t_diff == 1) %>% 
  filter(cage == cage_b)%>%
  rstatix::tukey_hsd(dist ~gm * as.factor(time_point_b)) %>% 
  as_tibble() %>% 
  write_tsv('stats/room_temp/dist_from_previous_measurement_tukey_hsd_results.tsv')

## Format data for Room Temp Study PCoA
time_samples <- time_metadata
time_bc_dist <- dist_subset(room_temp_bc, time_samples$sampleid )
time_room_temp_pcoa <- generate_pcoa(time_bc_dist)

time_room_temp_pcoa_data <- time_room_temp_pcoa[[1]]
time_room_temp_pcoa_p_var <- time_room_temp_pcoa[[2]]

fig2c_time_pcoa <- time_room_temp_pcoa_data %>%
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         order = case_when(time_point == '7a' ~ 0,
                           time_point == '8a' ~ 1,
                           time_point == '9a' ~ 2,
                           time_point == '10a' ~ 3,
                           time_point == '1p' ~ 6,
                           time_point == '4p' ~ 9),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = factor(gm, levels = c('GM-Low', 'GM-High')), alpha = time_point), size = 4) +
  scale_color_manual(values = c("red", "blue"), name = 'GM',
                     labels = c(expression(bold('GM'['Low']),
                                           bold('GM'['High'])))) +  scale_alpha_continuous(breaks = c(0, 1,2,3,6,9), name = 'Hours Post\nCollection') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(time_room_temp_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(time_room_temp_pcoa_p_var[2],2)}%')) +
  guides(
  ) +
  theme(
    legend.text = element_text(hjust = 0.5,face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0)
  ) +
  guides(
    color = guide_legend(order = 1),
    alpha = guide_legend(order = 2)
  )

# Room temp PCoA for main figure
time_all_res <- adonis2(time_bc_dist ~ gm * time_point, data = time_room_temp_pcoa_data,
        permutations = 9999,method = 'bray' )
time_all_res %>% 
  as_tibble(rownames = 'effect') %>% 
  mutate(SumOfSqs = round(SumOfSqs, 4),
         R2 = round(R2, 4),
         F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

## Begin parsing by GM
# Start with GM-Low
gm1_samples <- time_metadata %>% 
  filter(gm == 'GM-Low')
gm1_bc_dist <- dist_subset(room_temp_bc, gm1_samples$sampleid )
gm1_room_temp_pcoa <- generate_pcoa(gm1_bc_dist)

gm1_room_temp_pcoa_data <- gm1_room_temp_pcoa[[1]]
gm1_room_temp_pcoa_p_var <- gm1_room_temp_pcoa[[2]]

figS2a_gm1_pcoa <- gm1_room_temp_pcoa_data %>%
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         order = case_when(time_point == '7a' ~ 0,
                           time_point == '8a' ~ 1,
                           time_point == '9a' ~ 2,
                           time_point == '10a' ~ 3,
                           time_point == '1p' ~ 6,
                           time_point == '4p' ~ 9),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = gm, alpha = time_point), size = 4) +
  scale_color_manual(values = c("red"), 'GM') +
  scale_alpha_continuous(breaks = c(0,1,2,3,6,9), name = 'Hours Post\nCollection') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(gm1_room_temp_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(gm1_room_temp_pcoa_p_var[2],2)}%')) +
  guides(
    color = 'none'
  ) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    aspect.ratio = 2/2
  ) +
  guides(
    color = 'none',
    shape = guide_legend(order = 2),
    alpha = guide_legend(order = 3)
  )

# GM-Low PERMANOVA for supplemental
gm_low_adonis_res <- vegan::adonis2(gm1_bc_dist ~ time_point, data =  gm1_samples,
                                    permutations = 9999)
gm_low_adonis_res %>% 
  as_tibble(rownames = 'effect') %>% 
  select(effect, R2, F, `Pr(>F)`) %>% 
  mutate(R2 = round(R2, 4),
         F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

# Next with GM-High
gm4_samples <- time_metadata %>% 
  filter(gm == 'GM-High')
gm4_bc_dist <- dist_subset(room_temp_bc, gm4_samples$sampleid )
gm4_room_temp_pcoa <- generate_pcoa(gm4_bc_dist)

gm4_room_temp_pcoa_data <- gm4_room_temp_pcoa[[1]]
gm4_room_temp_pcoa_p_var <- gm4_room_temp_pcoa[[2]]

figS2b_gm4_pcoa <- gm4_room_temp_pcoa_data %>%
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                                time_point == '8a' ~ 1,
                                time_point == '9a' ~ 2,
                                time_point == '10a' ~ 3,
                                time_point == '1p' ~ 6,
                                time_point == '4p' ~ 9),
         order = case_when(time_point == '7a' ~ 0,
                           time_point == '8a' ~ 1,
                           time_point == '9a' ~ 2,
                           time_point == '10a' ~ 3,
                           time_point == '1p' ~ 6,
                           time_point == '4p' ~ 9),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = gm, alpha = time_point), size = 4) +
  scale_color_manual(values = c("blue"), 'GM') +
  scale_alpha_continuous(breaks = c(0, 1,2,3,6,9), name = 'Hours Post\nCollection') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(gm4_room_temp_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(gm4_room_temp_pcoa_p_var[2],2)}%')) +
  guides(
    color = 'none'
  ) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    aspect.ratio = 2/2,
  ) +
  guides(
    color = 'none',
    shape = guide_legend(order = 2),
    alpha = guide_legend(order = 3)
  )

# PERMANOVA for GM-high
gm_high_adonis_res <- vegan::adonis2(gm4_bc_dist ~time_point, data =  gm4_samples,
                                    permutations = 9999)
gm_high_adonis_res %>% 
  as_tibble(rownames = 'effect') %>% 
  select(effect, R2, F, `Pr(>F)`) %>% 
  mutate(R2 = round(R2, 4),
         F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()


## Generate pairwise distances between time points
time_a <- time_metadata %>% 
  select(sampleid, cage, order, time_point) 
time_b <- time_metadata %>% 
  select(sampleid, time_point, cage, order, gm) %>% 
  rename(b = 'sampleid',
         time_point_b = 'time_point',
         order_b = 'order',
         cage_b = 'cage')

figS2c_gm1_bc_mat <- bc_dist_long %>% 
  left_join(., time_a, by = 'sampleid') %>% 
  left_join(., time_b, by = 'b') %>% 
  filter(cage == cage_b) %>% 
  group_by(gm, time_point, time_point_b) %>% 
  summarize(mean_dist = mean(dist)) %>% 
  filter(gm == 'GM-Low') %>%
  ggplot(aes(x = factor(time_point), y = factor(time_point_b), fill = mean_dist)) +
  geom_tile() +
  geom_text(aes(label = round(mean_dist, 3)), fontface = 'bold') +
  coord_flip() +
  scale_fill_gradient2(low = 'white', high = 'red', limits = c(0,1), na.value = 'white', name = 'Avg. Intracage\nBray-Curtis\nDissimilarity') +
  ggprism::theme_prism() +
  labs(x = 'Hours Post Collection', y = 'Hours Post Collection') +
  theme(
    legend.text = element_text(face = "bold", color = 'black', size = 10),
    legend.title = element_text(face = 'bold', color = 'black', size = 10, hjust = 0),
    strip.text = element_text(size = 10),
    aspect.ratio = 2/2
  )

figS2d_gm4_bc_mat <- bc_dist_long %>% 
  left_join(., time_a, by = 'sampleid') %>% 
  left_join(., time_b, by = 'b') %>% 
  filter(cage == cage_b) %>% 
  group_by(gm, time_point, time_point_b) %>% 
  summarize(mean_dist = mean(dist)) %>%
  filter(gm == 'GM-High') %>%
  ggplot(aes(x = factor(time_point), y = factor(time_point_b), fill = mean_dist)) +
  geom_tile() +
  geom_text(aes(label = round(mean_dist, 3)), fontface = 'bold') +
  coord_flip() +
  scale_fill_gradient2(low = 'white', high = 'blue', limits = c(0,1), na.value = 'white', name = 'Avg. Intracage\nBray-Curtis\nDissimilarity') +
  ggprism::theme_prism() +
  labs(x = 'Hours Post Collection', y = 'Hours Post Collection') +
  theme(
    legend.text = element_text(face = "bold", color = 'black', size = 10),
    legend.title = element_text(face = 'bold', color = 'black', size = 10, hjust = 0),
    strip.text = element_text(size = 10),
    aspect.ratio = 2/2
  )

bc_dist_long %>% 
  left_join(., time_a, by = 'sampleid') %>% 
  left_join(., time_b, by = 'b') %>% 
  filter(cage == cage_b) %>% 
  group_by(gm, time_point, time_point_b) %>% 
  summarize(mean_dist = mean(dist)) %>% 
  group_by(gm) %>% 
  summarise(mean = round(mean(mean_dist),4),
            sd = round(sd(mean_dist),4)) %>% 
  clipr::write_clip()

####
# Position Data
####

# Gather metadata
position_metadata <- metadata %>% 
  filter(study == 'colon_position') %>% 
  mutate(order = case_when(sample_type == 'cecum' ~ 1,
                           time_point == 'proximal' ~ 2,
                           time_point == 'mid' ~ 3,
                           time_point == 'distal' ~ 4),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High'))

bc_dist_samples <- attributes(bc_dist)$Labels

# Pull GM-High samples
gm4_position_samples <- position_metadata %>% 
  filter(gm == 'GM-High')

gm4_dist_samples <- gm4_position_samples %>% 
  filter(sampleid %in% bc_dist_samples)
gm4_bc_dist <- dist_subset(bc_dist, gm4_dist_samples$sampleid)
gm4_position_pcoa <- generate_pcoa(gm4_bc_dist)

gm4_position_pcoa_data <- gm4_position_pcoa[[1]] %>% 
  mutate(sample_type = str_to_sentence(sample_type))
gm4_position_pcoa_p_var <- gm4_position_pcoa[[2]]

figS4b_gm4_pos_pcoa <- gm4_position_pcoa_data %>%
  mutate(gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = interaction(time_point), 
                 shape = factor(sample_type, levels = c('Cecum', 'Proximal', 'Mid', 'Distal'))), 
             size = 4) +
  scale_color_manual(values = c('blue',  '#9696FF'), 
                     name = 'Timepoint', 
                     labels = c('AM',
                                'PM')) +
  scale_shape_manual(values = c(15,16,17,18),
                     name = 'Sample Type') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(gm4_position_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(gm4_position_pcoa_p_var[2],2)}%')) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14)
    
  )

# next pull samples for GM-Low
gm1_position_samples <- position_metadata %>% 
  filter(gm == 'GM-Low')

gm1_dist_samples <- gm1_position_samples %>% 
  filter(sampleid %in% bc_dist_samples)
gm1_bc_dist <- dist_subset(bc_dist, gm1_dist_samples$sampleid)
gm1_position_pcoa <- generate_pcoa(gm1_bc_dist)

gm1_position_pcoa_data <- gm1_position_pcoa[[1]] %>% 
  mutate(sample_type = str_to_sentence(sample_type))
gm1_position_pcoa_p_var <- gm1_position_pcoa[[2]]

figS4a_gm1_pos_pcoa <- gm1_position_pcoa_data %>%
  mutate(gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = interaction(time_point), 
                 shape = factor(sample_type, levels = c('Cecum', 'Proximal', 'Mid', 'Distal'))), size = 4) +
  scale_color_manual(values = c('red','#FFA5A5'), 
                     name = 'Timepoint', 
                     labels = c('AM',
                                'PM')) +
  scale_shape_manual(values = c(15,16,17,18), name ='Sample Type') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(gm1_position_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(gm1_position_pcoa_p_var[2],2)}%')) +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14),
  )

adonis2(gm1_bc_dist ~ sample_type * time_point, 
        data =gm1_position_pcoa_data,
        permutations = 9999) %>% 
  as_tibble(rownames = 'effect') %>% 
  select(effect, F, `Pr(>F)`) %>% 
  mutate(F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

adonis2(gm4_bc_dist ~ sample_type * time_point, 
        data =gm4_position_pcoa_data,
        permutations = 9999) %>% 
  as_tibble(rownames = 'effect') %>% 
  select(effect, F, `Pr(>F)`) %>% 
  mutate(F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

## Create pairwise comparisons between sample locations 
# First in GM_High
bc_dist_long <- bc_dist %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid") %>% 
  pivot_longer(-sampleid,
               values_to = "dist",
               names_to = "b") %>% 
  filter(sampleid %in% position_metadata$sampleid  &
           b %in% position_metadata$sampleid)

room_temp_bc <- bc_dist_long %>% 
  pivot_wider(names_from = 'b', values_from = 'dist') %>% 
  column_to_rownames(var = 'sampleid') %>% 
  as.matrix()

sample_a <- position_metadata %>% 
  select(sampleid, gm, animal,sex, sample_type, time_point, cage) %>% 
  rename(gm_a = 'gm',
         animal_a = 'animal',
         sample_type_a = 'sample_type',
         time_point_a = 'time_point')
sample_b <- position_metadata %>% 
  select(sampleid, gm, animal,sample_type, time_point) %>% 
  rename(b = 'sampleid',
         gm_b = 'gm',
         animal_b = 'animal',
         sample_type_b = 'sample_type',
         time_point_b = 'time_point')

figS4d_gm4_am_pm <- bc_dist_long %>% 
  left_join(., sample_a, by = 'sampleid') %>% 
  left_join(., sample_b, by = 'b') %>% 
  mutate(time_point_a = case_when(time_point_a == '7a' ~ 'AM',
                                  time_point_a == '4p' ~ 'PM'),
         sample_type_a = str_to_sentence(sample_type_a),
         sample_type_b = str_to_sentence(sample_type_b),) %>% 
  filter(animal_a == animal_b) %>% 
  group_by(gm_a, sample_type_a, sample_type_b, time_point_a, time_point_b) %>% 
  summarize(mean_dist = mean(dist)) %>% 
  filter(gm_a == 'GM-High') %>% 
  filter(sample_type_a != sample_type_b) %>% 
  ggplot(aes(x = factor(sample_type_a), y = factor(sample_type_b), fill = mean_dist)) +
  geom_tile() +
  geom_text(aes(label = round(mean_dist, 3)), fontface = 'bold')  +
  facet_wrap(~time_point_a, nrow = 2, scales = 'free') +
  coord_flip() +
  scale_fill_gradient2(low = 'white', high = 'blue', limits = c(0,1), name = 'Average\nIntrasubject\nDissimilarity') +
  ggprism::theme_prism() +
  theme(
    legend.text = element_text(face = "bold", color = 'black', size = 12),
    legend.title = element_text(face = 'bold', color = 'black', size = 12, hjust = 0),
    strip.text = element_blank(),
    aspect.ratio = 2/2,
    axis.title = element_blank()
  )

figS4c_gm1_am_pm  <- bc_dist_long %>% 
  left_join(., sample_a, by = 'sampleid') %>% 
  left_join(., sample_b, by = 'b') %>% 
  mutate(time_point_a = case_when(time_point_a == '7a' ~ 'AM',
                                  time_point_a == '4p' ~ 'PM'),
         sample_type_a = str_to_sentence(sample_type_a),
         sample_type_b = str_to_sentence(sample_type_b),) %>% 
  filter(animal_a == animal_b) %>% 
  filter(sample_type_a != sample_type_b) %>% 
  group_by(gm_a, sample_type_a, sample_type_b, time_point_a, time_point_b) %>% 
  summarize(mean_dist = mean(dist)) %>% 
  filter(gm_a == 'GM-Low') %>% 
  ggplot(aes(x = factor(sample_type_a), y = factor(sample_type_b), fill = mean_dist)) +
  geom_tile() +
  geom_text(aes(label = round(mean_dist, 3)), fontface = 'bold')  +
  facet_wrap(~time_point_a, nrow = 2, scales = 'free') +
  coord_flip() +
  scale_fill_gradient2(low = 'white', high = 'red', limits = c(0,1), name = 'Average\nIntrasubject\nDissimilarity') +
  ggprism::theme_prism() +
  theme(
    legend.text = element_text(face = "bold", color = 'black', size = 12),
    legend.title = element_text(face = 'bold', color = 'black', size = 12, hjust = 0),
    strip.text = element_blank(),
    aspect.ratio = 2/2,
    axis.title = element_blank()
  )

bc_dist_long %>% 
  left_join(., sample_a, by = 'sampleid') %>% 
  left_join(., sample_b, by = 'b') %>% 
  mutate(time_point_a = case_when(time_point_a == '7a' ~ 'AM',
                                  time_point_a == '4p' ~ 'PM'),
         sample_type_a = str_to_sentence(sample_type_a),
         sample_type_b = str_to_sentence(sample_type_b),) %>% 
  filter(animal_a == animal_b) %>% 
  group_by(gm_a, sample_type_a, sample_type_b, time_point_a, time_point_b) %>% 
  summarize(mean_dist = mean(dist)) %>% 
  group_by(gm_a, time_point_a) %>% 
  summarize(avg_dist = round(mean(mean_dist), 4), 
            sd_dist = round(std_error(mean_dist), 4)) %>% 
  clipr::write_clip()


# Plot sample distance from cecum
cecum_ids <- position_metadata %>%
  filter(sample_type == 'cecum')

b_metadata <- position_metadata %>% 
  rename(b = 'sampleid',
         animal_b = 'animal') %>% 
  select(b, animal_b)

fig4e_dist_from_cecum  <- bc_dist_long %>% 
  filter(b %in% cecum_ids$sampleid) %>% 
  right_join(., position_metadata, by = "sampleid") %>% 
  left_join(., b_metadata, by = 'b') %>% 
  mutate(sample_type = str_to_sentence(sample_type)) %>% 
  filter(animal == animal_b) %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 'AM',
                                time_point == '4p' ~ 'PM')) %>% 
  ggplot(aes(x = factor(sample_type, levels = c('Cecum', 'Proximal',
                                                'Mid', 'Distal')), 
             y = dist, 
             color = interaction(time_point,
                                 factor(gm, levels = c('GM-Low', 'GM-High'))
                                 ))) +
  stat_summary(fun = mean,  geom = "line", linewidth = 3,
               aes(group = interaction(gm, time_point),
                   linetype = sex))  +
  stat_summary(aes(group = interaction(gm, sample_type, time_point)),
               fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x), geom = 'errorbar',
               color = 'black', width = 0.2, linewidth = 1) +
  scale_color_manual(values = c('red','#FFA5A5', 'blue',  '#9696FF'), 
                     name = 'Timepoint/GM', 
                     labels = c(expression(bold('AM: GM'['Low'])),
                                expression(bold('PM: GM'['Low'])),
                                expression(bold('AM: GM'['High'])),
                                expression(bold('PM: GM'['High'])))) +
  scale_y_continuous(expand = expansion( mult = c(0, 0.01)), limits = c(0,1)) +
  scale_linetype_manual(values = c(1,2), name = 'Sex')+
  ggprism::theme_prism() +
  labs(x = "Sample Site", y = "Distance\nFrom Cecum") +
  theme(
    axis.text = element_text(face = 'bold'),
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14), 
    strip.text = element_text(size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14), 
    aspect.ratio = 2/2
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(order = 2, override.aes = list(linewidth = 2))
  )

# Distance from cecum stats for main figure
bc_dist_long %>% 
  filter(b %in% cecum_ids$sampleid) %>% 
  right_join(., position_metadata, by = "sampleid") %>% 
  left_join(., b_metadata, by = 'b') %>% 
  mutate(sample_type = str_to_sentence(sample_type)) %>% 
  filter(animal == animal_b) %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 'AM',
                                time_point == '4p' ~ 'PM')) %>% 
  filter(sample_type != 'Cecum') %>% 
  rstatix::anova_test(dist ~ gm * time_point * sample_type) %>% 
  clipr::write_clip()

bc_dist_long %>% 
  filter(b %in% cecum_ids$sampleid) %>% 
  right_join(., position_metadata, by = "sampleid") %>% 
  left_join(., b_metadata, by = 'b') %>% 
  mutate(sample_type = str_to_sentence(sample_type)) %>% 
  filter(animal == animal_b) %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 'AM',
                                time_point == '4p' ~ 'PM')) %>% 
  filter(sample_type != 'Cecum') %>% 
  rstatix::tukey_hsd(dist ~ gm * time_point * sample_type) %>% 
  write_tsv('stats/colon_position/dist_from_cecum_tukey_post_hoc.tsv')


# Figure 4D PCoa
bc_dist_samples <- attributes(bc_dist)$Labels
pos_dist_samples <- position_metadata %>% 
  filter(sampleid %in% bc_dist_samples)
pos_bc_dist <- dist_subset(bc_dist, pos_dist_samples$sampleid)
pos_position_pcoa <- generate_pcoa(pos_bc_dist)

pos_position_pcoa_data <- pos_position_pcoa[[1]] %>% 
  mutate(sample_type = str_to_sentence(sample_type))
pos_position_pcoa_p_var <- pos_position_pcoa[[2]]

fig4d_pos_pcoa <- pos_position_pcoa_data %>%
  mutate(gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High'),
         time_point = case_when(time_point == '7a' ~ 'AM',
                                time_point == '4p' ~ 'PM')) %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(color = interaction(time_point, factor(gm, levels = c('GM-Low',
                                                                       'GM-High')) ), shape = factor(sample_type, levels = c('Cecum', 'Proximal', 'Mid', 'Distal'))), size = 4) +
  scale_color_manual(values = c('red','#FFA5A5', 'blue',  '#9696FF'), 
                     name = 'Timepoint/GM', 
                     labels = c(expression(bold('AM: GM'['Low'])),
                                expression(bold('PM: GM'['Low'])),
                                expression(bold('AM: GM'['High'])),
                                expression(bold('PM: GM'['High'])))) +
  scale_shape_manual(values = c(15,16,17,18),
                     name = 'Sample Type') +
  ggprism::theme_prism() +
  labs(x = glue::glue('PCo1 - {round(pos_position_pcoa_p_var[1],2)}%'),
       y = glue::glue('PCo2 - {round(pos_position_pcoa_p_var[2],2)}%')) +
  guides(
    color = guide_legend(order = 1),
  ) +
  theme(
    axis.text = element_text(face = 'bold'),
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14), 
    strip.text = element_text(size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14), 
    aspect.ratio = 2/2
  )

pos_all_adonis_res <- adonis2(pos_bc_dist ~ gm * sample_type * time_point, data =pos_position_pcoa_data,
        permutations = 9999, method = 'bray')

pos_all_adonis_res %>% 
  as_tibble(rownames = 'effect') %>% 
  mutate(SumOfSqs = round(SumOfSqs, 4),
         R2 = round(R2, 4),
         F = round(F, 4),
         `Pr(>F)` = round(`Pr(>F)`, 4)) %>% 
  clipr::write_clip()

# Plot PCoA per sample location
# Start with GM-High
gm4_tissue_pcoa <- function(tissue, shape){
  gm4_dist_samples <- gm4_position_samples %>% 
    filter(sampleid %in% bc_dist_samples &
             sample_type == tissue )
  gm4_bc_dist <- dist_subset(bc_dist, gm4_dist_samples$sampleid)
  gm4_position_pcoa <- generate_pcoa(gm4_bc_dist)
  
  gm4_position_pcoa_data <- gm4_position_pcoa[[1]] %>% 
    mutate(sample_type = str_to_sentence(sample_type))
  
  am_ids <- gm4_position_pcoa_data %>% 
    filter(time_point == '7a') %>% 
    pull(sampleid)
  pm_ids <- gm4_position_pcoa_data %>% 
    filter(time_point == '4p') %>% 
    pull(sampleid)
  cd_value <- dist_between_centroids(gm4_bc_dist,am_ids, pm_ids) 
  cd_label <- paste('CD = ', round(cd_value, 3), sep = "")
  
  centroid <- gm4_position_pcoa_data %>% 
    group_by(time_point) %>% 
    summarize(x = mean(PCo1),
              y = mean(PCo2))
  pm_centroid <- centroid[1,] %>% 
    rename(x.pm = 'x', y.pm = 'y')
  am_centroid <- centroid[2,] %>% 
    rename(x.am = 'x', y.am = 'y')
  line <- cbind(am_centroid, pm_centroid) %>% 
    select(x.am, x.pm, y.am, y.pm)
  
  gm4_position_pcoa_p_var <- gm4_position_pcoa[[2]]
  adonis_res <- adonis2(gm4_bc_dist ~ time_point, data = gm4_position_pcoa_data, 
                        permutations = 9999,method = 'bray')
  f_value <- round(adonis_res$F[1],2)
  p_value <- round(adonis_res$`Pr(>F)`[1],4)

  gm4_pos_pcoa <- gm4_position_pcoa_data %>%
    mutate(gm = case_when(gm == 1 ~ 'GM-Low',
                          gm == 4 ~ 'GM-High')) %>% 
    ggplot(aes(x = PCo1, y = PCo2)) +
    geom_segment(aes(x = x.am, y = y.am, xend = x.pm, yend = y.pm), 
                 data = line, linewidth = 1, linetype = 1) +
    geom_point(aes(x = x.pm, y = y.pm), size = 4,
               color = 'blue',  shape = 7, stroke = 2,
               data = pm_centroid) +
    geom_point(aes(x = x.am, y = y.am), size = 4,
               color = '#9696FF',  shape = 7, stroke = 2,
               data = am_centroid) +
    geom_point(aes(color = time_point), size = 4, shape = shape) +
    scale_color_manual(values = c('blue',  '#9696FF'), 
                       name = expression(bold('GM'['High'])), 
                       labels = c('AM',
                                  'PM')) +
    annotate(geom = "text", 
             y= max(gm4_position_pcoa_data$PCo2) * 1.3, 
             x = max(gm4_position_pcoa_data$PCo1) * 1.1,
             label=cd_label, hjust=1, size = 4, 
             color = 'black', fontface = 'bold') +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    ggprism::theme_prism() +
    labs(x = glue::glue('PCo1 - {round(gm4_position_pcoa_p_var[1],2)}%'),
         y = glue::glue('PCo2 - {round(gm4_position_pcoa_p_var[2],2)}%'),
         caption = glue::glue_col('F = {f_value}\np = {p_value}'),
         ) +
         
    guides(
      color = guide_legend(override.aes = list(shape = 15))
    ) +
    theme(
      legend.text = element_text(face = 'bold', color = 'black', size = 14),
      legend.title = element_text(face = 'bold', color = 'black', size = 14),
      plot.caption = element_text(face = 'bold', color = 'black', size = 12),
      plot.title = element_text(face = 'bold', color = 'black', size = 14),
      aspect.ratio = 1
    )

  return(gm4_pos_pcoa)
}

gm1_tissue_pcoa <- function(tissue, shape){
  gm1_dist_samples <- gm1_position_samples %>% 
    filter(sampleid %in% bc_dist_samples &
             sample_type == tissue)
  gm1_bc_dist <- dist_subset(bc_dist, gm1_dist_samples$sampleid)
  gm1_position_pcoa <- generate_pcoa(gm1_bc_dist)
  
  gm1_position_pcoa_data <- gm1_position_pcoa[[1]] %>% 
    mutate(sample_type = str_to_sentence(sample_type))
  
  am_ids <- gm1_position_pcoa_data %>% 
    filter(time_point == '7a') %>% 
    pull(sampleid)
  pm_ids <- gm1_position_pcoa_data %>% 
    filter(time_point == '4p') %>% 
    pull(sampleid)
  cd_value <- dist_between_centroids(gm1_bc_dist,am_ids, pm_ids) 
  cd_label <- paste('CD = ', round(cd_value, 3), sep = "")
  
  centroid <- gm1_position_pcoa_data %>% 
    group_by(time_point) %>% 
    summarize(x = mean(PCo1),
              y = mean(PCo2))
  pm_centroid <- centroid[1,] %>% 
    rename(x.pm = 'x', y.pm = 'y')
  am_centroid <- centroid[2,] %>% 
    rename(x.am = 'x', y.am = 'y')
  line <- cbind(am_centroid, pm_centroid) %>% 
    select(x.am, x.pm, y.am, y.pm)

  gm1_position_pcoa_p_var <- gm1_position_pcoa[[2]]
  adonis_res <- adonis2(gm1_bc_dist ~ time_point, data = gm1_position_pcoa_data, 
                        permutations = 9999,method = 'bray')  
  f_value <- round(adonis_res$F[1],2)
  
  p_value <- round(adonis_res$`Pr(>F)`[1], 4)

  gm1_pos_pcoa <- gm1_position_pcoa_data %>%
    mutate(gm = case_when(gm == 1 ~ 'GM-Low',
                          gm == 4 ~ 'GM-High')) %>% 
    ggplot(aes(x = PCo1, y = PCo2)) +
    geom_point(aes(color = time_point), size = 4, shape = shape) +
    geom_segment(aes(x = x.am, y = y.am, xend = x.pm, yend = y.pm), 
                 data = line, linewidth = 1, linetype = 1) +
    geom_point(aes(x = x.pm, y = y.pm), size = 4,
               color = 'red',  shape = 7, stroke = 2,
               data = pm_centroid) +
    geom_point(aes(x = x.am, y = y.am), size = 4,
               color = '#FFA5A5',  shape = 7, stroke = 2,
               data = am_centroid) +
    annotate(geom = "text", 
             y= max(gm1_position_pcoa_data$PCo2) * 1.1, 
             x = max(gm1_position_pcoa_data$PCo1) * 1.1,
             label=cd_label, hjust=1, size = 4, 
             color = 'black', fontface = 'bold') +
     scale_color_manual(values = c('red','#FFA5A5'), 
                       name = expression(bold('GM'['Low'])), 
                       labels = c('AM',
                                  'PM')) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    
    ggprism::theme_prism() +
    labs(x = glue::glue('PCo1 - {round(gm1_position_pcoa_p_var[1],2)}%'),
         y = glue::glue('PCo2 - {round(gm1_position_pcoa_p_var[2],2)}%'),
         caption = glue::glue_col('F = {f_value}\np = {p_value}'),
         title = glue::glue('{str_to_sentence(tissue)}')
         ) +
    guides(
      color = guide_legend(override.aes = list(shape = 15))
    ) +
    theme(
      legend.text = element_text(face = 'bold', color = 'black', size = 14),
      legend.title = element_text(face = 'bold', color = 'black', size = 14),
      plot.caption = element_text(face = 'bold', color = 'black', size = 12),
      plot.title = element_text(face = 'bold', color = 'black', size = 16),
      aspect.ratio = 1
    )
  return(gm1_pos_pcoa)
}

figS4e_f_per_tissue_pcoa <- (gm1_tissue_pcoa('cecum', 15) | gm1_tissue_pcoa('proximal', 16) |
    gm1_tissue_pcoa('mid', 17) | gm1_tissue_pcoa('distal', 18)) /
  (gm4_tissue_pcoa('cecum', 15) | gm4_tissue_pcoa('proximal', 15) |
     gm4_tissue_pcoa('mid', 17) | gm4_tissue_pcoa('distal', 18)) +
  plot_layout(guides = 'collect')
