library(tidyverse)
library(readxl)
library(rstatix)

std_error <- function(x) sd(x)/sqrt(length(x))

colon_data <- read_excel('data/colon_lengths.xlsx')

figS3_fecal_boli <- colon_data %>% 
  ggplot(aes(x = gm, y = fecal_boli)) +
  geom_point(aes(color = interaction(gm, time_of_day)), size = 4,
             position = position_jitterdodge(jitter.width = 0.35,
                                             dodge.width = 0.6)) +
  stat_summary(fun = mean, geom = 'bar', aes(group = time_of_day),
               width = 0.5, position = position_dodge(0.6),
               fill = NA, color = 'black', linewidth = 1) +
  stat_summary(fun.max = function(x) mean(x) + std_error(x),
               fun.min = function(x) mean(x) - std_error(x),
               geom = 'errorbar', width = 0.1,
               position = position_dodge(0.6),
               aes(group = time_of_day), linewidth = 1) +
  scale_color_manual(values = c('red', 'blue', '#FFA5A5', '#9696FF'), name = 'Timepoint/GM',
                     labels = c(expression(bold('AM: GM'['Low'])),
                                expression(bold('AM: GM'['High'])),
                                expression(bold('PM: GM'['Low'])),
                                expression(bold('PM: GM'['High'])))) +
  scale_shape_manual(values = c(16,17), name = 'Sex')+
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.1))) +
  scale_x_discrete(labels = c(expression(bold('GM'['Low'])),
                              expression(bold('GM'['High'])))) +
  ggprism::theme_prism() +
  labs(y = glue::glue('Number of Fecal Boli')) +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_text(hjust = 0.5, face = 'bold', color = 'black',
                                size = 14),
    legend.text = element_text(hjust = 0, face = 'bold', color = 'black',
                               size = 14),
    aspect.ratio = 3/2,
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    alpha = guide_legend(order = 3)
  )

colon_data %>% 
  anova_test(fecal_boli ~ gm * time_of_day) %>% 
  clipr::write_clip()

colon_data %>% 
  tukey_hsd(fecal_boli ~ gm * time_of_day) %>% 
  clipr::write_clip()


  
