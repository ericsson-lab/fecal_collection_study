{
  source('code/data_alpha_stats.R')
  library(patchwork)
  library(rstatix)
}

# Order samples
alpha_stats$sample_type <- factor(alpha_stats$sample_type, levels = c('cecum', 'proximal',
                                                                      'mid', 'distal', 'feces'))

alpha_stats_updated <- alpha_stats %>% 
  mutate(time_point = case_when(time_point == '7a' ~ 0,
                               time_point == '8a' ~ 1,
                               time_point == '9a' ~ 2,
                               time_point == '10a' ~ 3,
                               time_point == '1p' ~ 6,
                               time_point == '4p' ~ 9),
         gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High'),
         sex = case_when(sex == 'M' ~ 'Male',
                         sex == 'F' ~ 'Female'))

alpha_stats_updated$gm = factor(alpha_stats_updated$gm, levels = c('GM-Low','GM-High' ))

# Intial alpha diversity figures, pulling t0 from both studies
plot_alpha_bar <- function(stat, y_lab) { 
  alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>%
  ggplot(aes(x = gm, y = {{stat}})) +
  geom_point(aes(color = gm), size = 3,
             position = position_jitterdodge(jitter.width = 0.5,
                                             dodge.width = 0.4)) +
  stat_summary(fun = mean, geom = 'bar',
               width = 0.5, position = position_dodge(0.4),
               fill = NA, color = 'black', linewidth = 1) +
  stat_summary(fun.max = function(x) mean(x) + std_error(x),
               fun.min = function(x) mean(x) - std_error(x),
               geom = 'errorbar', width = 0.1,
               position = position_dodge(0.4),
               linewidth = 1) +
  scale_color_manual(values = c('red', 'blue'), name = 'GM',
                     labels = c(expression(bold('GM'['Low']),
                                           bold('GM'['High'])))) +
  scale_shape_manual(values = c(16,17), name = 'Sex')+
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.1))) +
  scale_x_discrete(labels = c(c(expression(bold('GM'['Low']),
                                           bold('GM'['High']))))) +
  ggprism::theme_prism() +
  labs(y = glue::glue('{y_lab} Index')) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14),
    legend.title = element_text(face = 'bold', color = 'black',
                                size = 14, hjust = 0),
    legend.text = element_text(face = 'bold', color = 'black',
                                size = 14),
    legend.position = 'none',
    aspect.ratio = 2,
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    alpha = guide_legend(order = 3)
  )
}


fig1a_shannon <- plot_alpha_bar(shannon, 'Shannon')
fig1b_chao1 <- plot_alpha_bar(chao1, 'Chao1')


# Test for normality
alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>% 
  shapiro_test(chao1) # fail
alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>% 
  rstatix::shapiro_test(shannon) # pass

# T0 alpha diveristy stats
alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>% 
  rstatix::wilcox_test(chao1 ~ gm) %>% 
  clipr::write_clip()
alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>% 
  rstatix::t_test(shannon ~ gm, detailed = F) %>% 
  clipr::write_clip()


plot_alpha_bar_supplemental <- function(stat, y_lab) { 
  alpha_stats_updated %>% 
    filter(time_point == 0 & 
             sample_type %in% c('distal', 'feces')) %>%
    ggplot(aes(x = gm, y = {{stat}})) +
    geom_point(aes(color = gm, shape = sex), size = 3,
               position = position_jitterdodge(jitter.width = 0.25,
                                               dodge.width = 0.6)) +
    stat_summary(fun = mean, geom = 'bar',
                 width = 0.4, position = position_dodge(0.6),
                 fill = NA, color = 'black', linewidth = 1,
                 aes(group = sex)) +
    stat_summary(fun.max = function(x) mean(x) + std_error(x),
                 fun.min = function(x) mean(x) - std_error(x),
                 geom = 'errorbar', width = 0.1,
                 position = position_dodge(0.6),
                 linewidth = 1, aes(group = sex)) +
    scale_color_manual(values = c('red', 'blue'), name = 'GM',
                       labels = c(expression(bold('GM'['Low']),
                                             bold('GM'['High'])))) +
    scale_shape_manual(values = c(16,17), name = 'Sex')+
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.1))) +
    scale_x_discrete(labels = c(c(expression(bold('GM'['Low']),
                                             bold('GM'['High']))))) +
    ggprism::theme_prism() +
    labs(y = glue::glue('{y_lab} Index')) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14),
      legend.title = element_text(face = 'bold', color = 'black',
                                  size = 14, hjust = 0),
      legend.text = element_text(face = 'bold', color = 'black',
                                 size = 14),
      legend.position = 'none',
      aspect.ratio = 1.75,
    ) +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 2),
      alpha = guide_legend(order = 3)
    )
}

figS1b_shannon <- plot_alpha_bar_supplemental(shannon, 'Shannon')
figS1a_chao1 <- plot_alpha_bar_supplemental(chao1, 'Chao1')


alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>% 
  rstatix::anova_test(chao1 ~ gm * sex) %>% 
  as_tibble() %>% 
  select(Effect, F, p, `p<.05`) %>% 
  clipr::write_clip()

alpha_stats_updated %>% 
  filter(time_point == 0 & 
           sample_type %in% c('distal', 'feces')) %>% 
  rstatix::anova_test(shannon ~ gm * sex) %>% 
  as_tibble() %>% 
  select(Effect, F, p, `p<.05`) %>%  
  clipr::write_clip()

#####
## Room temperature study
#####

#Fig 2A Chao1 across time.
fig2a_chao1_plot_time <- alpha_stats_updated %>% 
  filter(study == 'room_temp') %>% 
  ggplot(aes(x = time_point, y = chao1, color = gm)) +
  stat_summary(fun = mean,  geom = "line", linewidth = 2,
               aes(group =as.factor(gm), color = as.factor(gm),
                   linetype = sex)) +
  stat_summary(aes(group = gm),
               fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x), geom = 'errorbar',
               color = 'black', width = 0.2, linewidth = 1) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 6, 9)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                     limits = c(0, 400)) +
  scale_linetype_manual(values = c(1, 2), name = "Sex") +
  scale_color_manual(values = c('red', 'blue'), name = "GM") +
  labs(x = 'Hours Post Collection',
       y = 'Chao1 Index') +
  ggprism::theme_prism() +
  theme(
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    strip.text = element_text(size = 14),
    legend.position = 'none'
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(override.aes = list(linewidth = 2))
  )

alpha_stats_updated %>% 
  filter(study == 'room_temp') %>% 
  group_by(gm, time_point) %>% 
  count()


# Figure 2B shannon across time
fig2b_shannon_plot_time <- alpha_stats_updated %>% 
  filter(study == 'room_temp') %>% 
  ggplot(aes(x = time_point, y = shannon, color = gm)) +
  stat_summary(fun = mean,  geom = "line", linewidth = 2,
               aes(group =as.factor(gm), color = as.factor(gm),
                   linetype = sex)) +
  stat_summary(aes(group = gm),
               fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x), geom = 'errorbar',
               color = 'black', width = 0.2, linewidth = 1)+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 6, 9)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 5)) +
  scale_linetype_manual(values = c(1, 2), name = "Sex") +
  scale_color_manual(values = c('red', 'blue'), name = "GM",
                     labels = c(expression(bold('GM'['Low']),
                                           bold('GM'['High'])))) +
  labs(x = 'Hours Post Collection',
       y = 'Shannon Index') +
  ggprism::theme_prism() +
  theme(
    legend.text = element_text(face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    strip.text = element_text(size = 14)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(override.aes = list(linewidth = 2))
  )

alpha_stats_updated %>% 
  filter(study == 'room_temp') %>% 
  rstatix::anova_test(chao1 ~ gm * time_point) %>% 
  as_tibble() %>% 
  select(Effect, F, p) %>% 
  clipr::write_clip()

alpha_stats_updated %>% 
  filter(study == 'room_temp') %>% 
  rstatix::anova_test(shannon ~ gm * time_point) %>% 
  as_tibble() %>% 
  select(Effect, F, p) %>% 
  clipr::write_clip()


#####
##### Colon Position
#####
#####

# Order sample position factors
alpha_stats$sample_type <- factor(str_to_sentence(alpha_stats$sample_type), levels = c('Cecum', 'Proximal', 'Mid', 'Distal'))

# format alpha stats for colon position study
alpha_stats_gm <- alpha_stats %>% 
  filter(study == 'colon_position') %>% 
  mutate(gm = case_when(gm == 1 ~ 'GM-Low',
                        gm == 4 ~ 'GM-High'),
         time_point = case_when(time_point == '7a' ~ 'AM',
                                time_point == '4p' ~ 'PM'))
# Order GM factor
alpha_stats_gm$gm <- factor(alpha_stats_gm$gm, levels = c('GM-Low',
                                                          'GM-High'))
# Order time factor
alpha_stats_gm$time_point <- factor(alpha_stats_gm$time_point, levels = c('AM', 'PM'))

# Plot figure 4B 
fig4b_chao1<-alpha_stats_gm %>%
  ggplot(aes(x = sample_type, y = chao1, 
             color = interaction(time_point, gm))) +
  stat_summary(fun = mean,  geom = "line", linewidth = 3,
               aes(group = interaction(gm, time_point),
                   linetype = sex))  +
  stat_summary(aes(group = interaction(gm, sample_type, time_point)),
               fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x), geom = 'errorbar',
               color = 'black', width = 0.2, linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0,400)) +
  scale_color_manual(values = c('red','#FFA5A5', 'blue',  '#9696FF'), 
                     name = 'Timepoint/GM', 
                     labels = c(expression(bold('AM: GM'['Low'])),
                                expression(bold('PM: GM'['Low'])),
                                expression(bold('AM: GM'['High'])),
                                expression(bold('PM: GM'['High'])))) +
  scale_linetype_manual(values = c(1,2), name = 'Sex') +
  labs(x = 'Colon Position',
       y = 'Chao1 Index') +
  ggprism::theme_prism() +
  theme(
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    strip.text = element_text(size = 14),
    aspect.ratio = 2/2
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(override.aes = list(linewidth = 2))
  )

# Plot figure 4C 
fig4c_shannon <- alpha_stats_gm %>%
  mutate(sex = ifelse(sex == 'M', 'Male', 'Female')) %>% 
  ggplot(aes(x = sample_type, y = shannon, 
             color = interaction(time_point, gm))) +
  stat_summary(fun = mean,  geom = "line", linewidth = 3,
               aes(group = interaction(gm, time_point),
                   linetype = sex))  +
  stat_summary(aes(group = interaction(gm, sample_type, time_point)),
               fun.min = function(x) mean(x) - std_error(x),
               fun.max = function(x) mean(x) + std_error(x), geom = 'errorbar',
               color = 'black', width = 0.2, linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, 5)) +
  scale_color_manual(values = c('red','#FFA5A5', 'blue',  '#9696FF'), 
                     name = 'Timepoint/GM', 
                     labels = c(expression(bold('AM: GM'['Low'])),
                                expression(bold('PM: GM'['Low'])),
                                expression(bold('AM: GM'['High'])),
                                expression(bold('PM: GM'['High'])))) +
  scale_linetype_manual(values = c(1,2), name = 'Sex') +
  labs(x = 'Colon Position',
       y = 'Shannon Index') +
  ggprism::theme_prism() +
  theme(
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black', size = 14),
    legend.title = element_text(face = 'bold', color = 'black', size = 14, hjust = 0),
    strip.text = element_text(size = 14),
    aspect.ratio = 2/2
    ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
    linetype = guide_legend(override.aes = list(linewidth = 2))
  )

# Chao1 ~ GM * Sample Type * Time Point
alpha_stats_gm %>% 
  rstatix::anova_test(chao1 ~ gm * sample_type * time_point, effect.size = 'pes') %>% 
  as_tibble() %>% 
  select(Effect, F, p, `p<.05`, pes) %>% 
  clipr::write_clip()

alpha_stats_gm %>% 
  rstatix::tukey_hsd(chao1 ~ gm * sample_type * time_point) %>% 
  as_tibble() %>% 
  write_tsv('stats/colon_position/chao1_tukey_post_hoc.tsv')

# Shannon ~ GM * Sample Type * Time Point
alpha_stats_gm %>% 
  rstatix::anova_test(shannon ~ gm * sample_type * time_point, effect.size = 'pes') %>% 
  as_tibble() %>% 
  select(Effect, F, p, `p<.05`, pes) %>% 
  clipr::write_clip()

alpha_stats_gm %>% 
  rstatix::tukey_hsd(shannon ~ gm * sample_type * time_point) %>% 
  as_tibble() %>% 
  write_tsv('stats/colon_position/shannon_tukey_post_hoc.tsv')


