setwd("/home/lxk/private/Radial-C/260719_RadialC_K562_FA_EGS_0.5ulMNase19m_ligaRT_5m,30m/output_RadialC_K562_equal_5m_R1r1/layer/")
library(tidyverse)
# library(Cairo)

df <- read.table("paired_layer_count.tab", sep = "\t", header = FALSE)
colnames(df) <- c("layer", "bin_num", "12bp", "18bp", "34bp", "60bp", "91bp")

df <- df %>%
  mutate(
    # 批量对 12bp 到 91bp 这5列除以 bin_num
    across(c("12bp", "18bp", "34bp", "60bp", "91bp"), ~ .x / bin_num)
  ) %>%
  select(-bin_num) %>%  # 计算完成后，丢弃 bin_num 列
  pivot_longer(
    cols = c("12bp", "18bp", "34bp", "60bp", "91bp"),
    names_to = "type",
    values_to = "frequency"      
  )

# df$dis <- factor(df$dis, levels = c("40", "35", "30", "25", "20"))

# CairoPDF("bio_vs_link.pdf", width = 2.5, height = 2, bg = "transparent", family ="Arial")
png("paired_layer_frequency.png", units = "in", width = 2.5, height = 2.5, res = 300, bg = "white")
# pdf("bio_vs_link.pdf", width = 2.0, height = 2.3, bg = "transparent")

ggplot(data = df, mapping = aes(x=layer, y=frequency, color=type)) + 
  geom_line(aes(group=type), size = 1) + 
  geom_point(size = 3) +
  labs(title = "1:1:1:1:1, 5 min") + 
  xlab("Layer") + 
  ylab("Normalized count") + 
  coord_cartesian(
    ylim = c(0, 160)
  )+
  scale_y_continuous(
    # limits = c(30, 120),
    expand = expansion(mult = c(0.00, .01))
  ) +
  scale_x_continuous(
    limits = c(1, 10),
    breaks = seq(1, 10),
    # labels = c("≤ 20",  '≤ 25', '≤ 30', '≤ 35', '≤ 40'),
    expand = expansion(mult = c(0.06, .04))
  ) +
  guides(
    color = guide_legend(
      ncol = 2,
      keywidth = 0.01,
    )
  ) +  
  scale_color_manual(
    limits = c("12bp", "18bp", "34bp", "60bp", "91bp"),
    values = c("12bp"="#4C72B0", "18bp"='#55A868', "34bp"='#C44E52', "60bp"='#CCB974', "91bp"='#798E87'),
    labels = c("12 bp", "18 bp", "34 bp", "60 bp", "91 bp"),
  ) +
  theme(
    axis.text = element_text(color = "black",size = 12), #face = "bold",
    axis.title = element_text(color = "black",size = 12),
    axis.text.x = element_text(color = "black", 
                               size = 12, 
                               # angle = 45, 
                               # hjust = 1,
                               # vjust = 1
    ),
    panel.border = element_blank(), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12,
                               margin = margin(l = 1, unit = "pt")
                               ),
    panel.background = element_blank(), 
    # legend.position = "right",
    legend.position = c(0.33,0.82),
    rect = element_rect(fill = NA, size=0, color = NA),
    plot.background = element_rect(fill = NA, size=0, color = NA),
    legend.background = element_rect(fill = NA, size=0, color = NA),
    legend.key = element_rect(fill = NA, size=0, color = NA),
    # plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    plot.title = element_text(
      hjust = 0.5,          # 水平对齐：0=左, 0.5=中, 1=右
      # vjust = 0.5,          # 垂直对齐
      # size = 14,            # 字体大小
      # face = "bold",        # 字体样式：bold, italic, plain
      # color = "darkblue",   # 字体颜色
      margin = margin(t = 0, b = 10)  # 上下边距
    ),
    legend.key.spacing.x = unit(0.1, "cm")
  ) 

dev.off()
