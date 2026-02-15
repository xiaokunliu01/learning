setwd("/mnt/disk1/6/lxk/private/MNase-seq/K562/R2/plot_Average_Profile_R1/")

library(tidyverse)
library(pracma)

df <- read.table("MNaseseq_K562_R1_CTCFmotif_profile.tab", sep = "\t", header = FALSE, skip = 2)

transposed_df <- as.data.frame(t(data.frame(df[,-2])))
new_colnames <- transposed_df[1, ]
transposed_df <- as.data.frame(transposed_df[-1, ])
colnames(transposed_df) <- new_colnames

bins <- seq(from = -1000, to = 999, by = 1)
# df <- cbind(bins, transposed_df)

df <- data.frame(
  bin = bins,
  value = as.numeric(transposed_df[,1])
)

x <- df[,1]
y <- df[,2]
# 找波峰（注意：pracma::findpeaks 找的是“正峰”，所以对 y 直接用）
pks <- findpeaks(y, 
                 minpeakheight = 70, 
                 minpeakdistance = 20,   # 峰之间至少隔3个点
                 threshold = 0)       # 峰必须比邻居高至少0.5

# pks 是一个矩阵，列：peak_value, index, start, end
peaks_idx <- pks[, 2]

# 波谷：对 -y 找峰
vlys <- findpeaks(-y, 
                  minpeakheight = -70,
                  minpeakdistance = 20, 
                  threshold = 0)
valleys_idx <- vlys[, 2]


png("r_K562_MNase_CTCFmotif_profile.png", units = "in", width = 5.6, height = 2.8, res = 300, bg = "white")
# pdf("enrichment_HiC_scale_500bp.pdf", width = 3.2, height = 3.9, bg = "transparent")

ggplot(data = df, mapping = aes(x = bin)) + 
  geom_line(aes(y = df[,2], color = "MNase-seq"), linewidth = 1.2) +
  geom_vline(
    # xintercept = c(-972, -882, -779, -697, -596, -511, -419, -312, -244, -137, -81, 82, 137, 254, 325, 416, 506, 592, 686, 775, 876, 962),
    xintercept = c(x[peaks_idx], x[valleys_idx]),
    linetype="dotted",
    ) + 
  annotate("text", 
           # x= c(-882, -697, -511, -312, -137, 137, 325, 506, 686, 876),
           x = x[peaks_idx],
           y = 78,
           # label = c(-882, -697, -511, -312, -137, 137, 325, 506, 686, 876)
           label = x[peaks_idx]
           ) +
  scale_color_manual(
    values = c("MNase-seq" = "#fc9533")
    ) + 
  guides(color=guide_legend(
    # labels = c("DNase-C", expression(paste(italic('in situ'), ' Hi-C')), "Micro-C", "BL-Hi-C", "Hi-TrAC"),
    nrow=1, 
    byrow=TRUE,
    # label.hjust = 1,
    keywidth = 0.45
    )
    ) +
  xlab("Distance from CTCF motif (bp)") + 
  ylab("Enrichment") +
  labs(title = "") +
  # coord_cartesian() +
  scale_x_continuous(
    limits = c(-1000, 1000), 
    expand =  expansion(mult = c(0, .03))
    ) +
  scale_y_continuous(
    limits = c(50, 80),
    # breaks = seq(0,0.25,0.05),
    # name = "DNase-C enrichment",
    # sec.axis = sec_axis(~. *5, name="DNase-seq enrichment"),
    expand = c(0,0)
  ) + 
  theme(
    axis.text = element_text(color = "black",size = 10), #face = "bold",
    axis.title = element_text(color = "black",size = 10),
    panel.border = element_blank(), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.background = element_blank(), 
    legend.position = "top",
    rect = element_rect(fill = NA, linewidth=0, color = NA),
    plot.background = element_rect(fill = NA, linewidth=0, color = NA),
    legend.background = element_rect(fill = NA, linewidth=0, color = NA),
    legend.key = element_rect(fill = NA, linewidth=0, color = NA),
    # plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) 

dev.off()


