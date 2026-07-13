setwd(dir = "/mnt/disk1/6/lxk/private/ChIP-seq/public/K562/EGR1/")

library(ggplot2)

# 准备数据框
df <- data.frame(
  group = c("contain", "no contain"),
  value = c(19602, 64832)
)


# 计算百分比用于标签
df$fraction <- df$value / sum(df$value)
df$percent <- round(df$fraction * 100, 2)
df$label <- paste0(df$group, "\n", df$percent, "%")

png("peak.png", units = "in", width = 2.5, height = 2.5, res = 300, bg = "white")
# pdf("EGR1-EGR1_contact_distance_R2.pdf", width = 4, height = 2.5, bg = "transparent")

# 绘制饼图
ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity") +  # 绘制堆叠条形图
  coord_polar("y", start = 0) +             # 将直角坐标系转换为极坐标系
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5)) +  # 调整标签位置
  scale_fill_manual(
    values = c("contain"="#386CAF",
               # "no contain"="#7FC87F"
               "no contain"="#9f79d3"
               # "no contain"=rgb(1,0.6,0.2)
    ),
    # limits = c("CTCF1_CTCF2", "CTCF2_EGR1", "CTCF1_EGR1"),
    # labels = c("CTCF1-CTCF2", "CTCF2-EGR1", "CTCF1-EGR1"),
  ) +
  guides(fill = guide_legend(
    title = NULL,
    # label.hjust = 0,
    # keywidth = unit(0.4,"cm"),
    # keyheight = unit(0.4, "cm"),
    ncol = 2,
    label.theme = element_text(size = 13),
  )) +
  labs(title = "EGR1 peak") + 
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
  
        # axis.text = element_text(color = "black", size = 15), 
        # face = "bold",
        # axis.title = element_text(color = "black", size = 15),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0.5,0.03),
        # legend.key.width = unit(0.2, "cm"),
        # legend.key.height = unit(0.2, "cm"),
        # legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        rect = element_rect(fill = NA, size=0, color = NA),
        plot.background = element_rect(fill = NA, size=0, color = NA),
        legend.background = element_rect(fill = NA, size=0, color = NA),
        legend.key = element_rect(fill = NA, size=0, color = NA),
        plot.title = element_text(
          hjust = 0.5,          # 水平对齐：0=左, 0.5=中, 1=右
          # vjust = 0.5,          # 垂直对齐
          # size = 14,            # 字体大小
          # face = "bold",        # 字体样式：bold, italic, plain
          # color = "darkblue",   # 字体颜色
          margin = margin(t = 0, b = -10)  # 上下边距
        )
  )

dev.off()
