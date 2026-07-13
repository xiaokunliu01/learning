library(ggplot2)
library(ggseqlogo)
library(dplyr) # 记得加载 dplyr 以使用 slice_sample

setwd("/mnt/disk1/6/lxk/private/DNase-C/dimer/K562/CTCF_EGR1/motif_analysis/")

df <- read.table("EGR_motif_fa.tab", sep = "\t", header = FALSE)
df <- df %>% slice_sample(n = 10000)

# 提取并转换序列
sequences <- toupper(df[,2])

png("EGR_motif_fa.png", units = "in", width = 6, height = 2, res = 300, bg = "white")
# pdf("THRB_0+.pdf", width = 4, height = 2, bg = "transparent")

ggplot() + 
  geom_logo(sequences, method = "bits") +
  scale_x_continuous(
    # limits = c(0,10),
    breaks=seq(1, nchar(sequences[1]), 5),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(0,2),
    expand = c(0,0)
  ) +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    # axis.text.x = element_text(color = "black", size = 10, angle = 0, hjust = 0.5, vjust = 1),
    axis.text = element_text(color = "black", size = 11), #face = "bold",
    axis.title = element_text(color = "black", size = 11),
    legend.position = "none",
    # legend.title = element_blank(),
    # legend.text = element_text(size = 12),
    rect = element_rect(fill = NA, linewidth=0, color = NA),
    plot.background = element_rect(fill = NA, linewidth=0, color = NA),
    legend.background = element_rect(fill = NA, linewidth=0, color = NA),
    legend.key = element_rect(fill = NA, linewidth=0, color = NA),
  )

dev.off()





