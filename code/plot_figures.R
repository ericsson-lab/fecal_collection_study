{
  library(patchwork)
  
  source('code/plot_alpha.R')
  source('code/plot_distances.R')
  source('code/plots_fecal_boli.R')
  source('code/plot_taxonomy.R')
  source('code/plot_pes.R')
  }


# panel A & B for figure 1
fig1b_chao1 + fig1a_shannon + fig1c_t0_pcoa + plot_layout(guides = 'collect')
ggsave('plots/room_temp/fig1_overall_chao1_shannon_pcoa.png', width = 10, height = 4, bg = 'white')


# Panel A & B for supplemental fig 1, add sex
(figS1a_chao1 + figS1b_shannon + figS1C_t0_pcoa_supplemental) / 
  (figS1d_gm1_sex_venn / figS1e_gm4_sex_venn) | figS1f_prevelance_abundance 
+ plot_layout(guides = 'collect')
ggsave('plots/room_temp/supp_fig1_overall_chao1_shannon.png', width = 6, height = 4, bg = 'white')

figS1d_gm1_sex_venn
figS1e_gm4_sex_venn
figS1f_prevelance_abundance

# Build Fig 2A-B using patchwork
fig2a_chao1_plot + theme(legend.position = 'none') +
  fig2b_shannon_plot + theme(
    legend.position = 'right',
    legend.text = element_text(hjust = 0)) &
  plot_annotation(tag_levels = 'A')

fig2c_time_pcoa 
fig2d_t0_distance 
fig2e_prev_t_distance


((chao1_plot_time + shannon_plot_time)/(time_pcoa | t0_distance / prev_t_distance)) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'keep', widths = c(0.5, 0.5), heights = c(0.4, 0.6)) & 
  theme(
    plot.tag = element_text(face = 'bold', color = 'black', size = 18 ),
    legend.title = element_text(hjust = 0),
    legend.text = element_text(hjust = 0)
  )

ggsave('plots/room_temp/fig_2_time.png', width = 11, height = 8, bg = 'white')

## Supplemetnal Fig 2
figS2a_gm1_pcoa
figS2b_gm4_pcoa
figS2c_gm1_bc_mat
figS2d_gm4_bc_mat

(gm1_pcoa + gm4_pcoa) /
  (gm1_bc_mat + gm4_bc_mat) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = 'bold', color = 'black', size = 18 )
  )
ggsave('plots/room_temp/gm1_gm4_bc_pairwise_distances.png', 
       width = 11, height = 8, bg = 'white')
# ggsave('plots/room_temp/gm4_bc_pairwise_distances.png', width = 5, height = 8)
ggsave('plots/room_temp/gm1_gm4_bc_pairwise_distances.png', width = 10, height = 9)


fig3_family
fig3_phylum


fig4e_dist_from_cecum
fig4b_chao1 + theme(legend.position = 'none') + 
  fig4c_shannon + theme(
    legend.position = 'right',
    legend.text = element_text(hjust = 0)) &
ggsave('plots/colon_position/chao1_shannon_main_fig.png', bg = 'white', 
       width = 10, height = 4)

figS3_fecal_boli

figS4a_gm1_pos_pcoa
figS4b_gm4_pos_pcoa
figS4c_gm1_am_pm
figS4d_gm1_am_pm
figS4e_f_per_tissue_pcoa
gm1_pos_pcoa + gm4_pos_pcoa +
  gm1_am_pm + gm4_am_pm +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(0.3, 0.7)) &
  theme(
    plot.tag = element_text(size = 18),
    
  )

library(patchwork)

gm1_pos_pcoa + gm4_pos_pcoa + plot_layout(guides = 'collect', )&
  guides(shape = guide_legend(order = 2),
         color = guide_legend(color = 1))
ggsave('plots/colon_position/gm_low_high_bc_pcoa.png', width = 10, height = 4, bg = 'white')




fig6a_time
fig6b_pos


