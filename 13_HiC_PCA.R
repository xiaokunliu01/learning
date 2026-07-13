library(ggplot2)

setwd('/mnt/disk1/6/lxk/private/DNase-C/dimer_paper/fig3/RAD21_degron/revise_R1/PCA/RAD21_CTCF_degron/')

data <- read.table("DMSO4h_dTAG4h_R1R2.tab")

countData <- as.matrix(data)
colnames(countData) <-c("DMSO_R1","DMSO_R2","dTAG_R1","dTAG_R2")


row_var <- apply(countData, 1, var)
nonzero_var_rows <- which(row_var != 0)
filtered_data <- countData[nonzero_var_rows,]
data <- as.data.frame(filtered_data)
gene <- t(data)


pca1 <- prcomp(gene[, -ncol(gene)], center = TRUE, scale. = TRUE)
df1 <- pca1$x
df1 <- as.data.frame(df1)

summ1 <- summary(pca1)
xlab1 <- paste0("PC1 (", round(summ1$importance[2, 1] * 100, 1), "%)")
ylab1 <- paste0("PC2 (", round(summ1$importance[2, 2] * 100, 1), "%)")

name <-c("DMSO","DMSO","dTAG","dTAG")
#
# gene <- cbind(gene, name)
# gene <- as.data.frame(gene)

ggplot(df1) +
  geom_point(aes(x = PC1, y = PC2, color = name), alpha = 1, size = 3) +
  # geom_hline(yintercept = 0, linetype="dashed", colour="#A9A9A9")+
  # geom_vline(xintercept = 0, linetype="dashed", colour="#A9A9A9")+
  labs(
    title = expression(Delta *RAD21*"; "*Delta *CTCF), 
    x=xlab1, 
    y=ylab1
    ) + 
  coord_cartesian(
    xlim = c(-110, 110),
    ylim = c(-110, 110),
  ) +
  scale_x_continuous(
    # limits = c(-100, 100), 
    breaks = seq(-100, 100, 100),
    # expand = expansion(mult = c(0, .04))
  ) +
  scale_y_continuous(
    # limits = c(-100, 100), 
    breaks = seq(-100, 100, 100),
    # expand = c(0,0)
  ) +
  scale_colour_manual(
    values = c("#386CAF", "#fc9533"),
    breaks = c("DMSO", "dTAG"),
    labels = c("DMSO", "dTAG")
  ) +
  theme(
    # panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    # axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text = element_text(color = "black", size = 12), #face = "bold",
    axis.title = element_text(color = "black", size = 12),
    legend.position = "right",
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),
    # rect = element_rect(fill = NA, linewidth = 0, color = NA),
    # plot.background = element_rect(fill = NA, linewidth = 0, color = NA),
    legend.background = element_rect(
      fill = NA,        # 背景填充色
      # color = "black",       # 边框颜色
      # size = 0.5,           # 边框宽度
      linetype = "solid"    # 边框线型
    ),
    legend.margin = margin(0, -5, 0, -5),
    legend.key = element_rect(fill = NA, linewidth = 0, color = NA),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.5, "cm"),
    plot.title = element_text(
      # face = "bold", 
      # color = "black", 
      # size = 16, 
      hjust = 0.5),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1,
      linetype = "solid"  # 或 "dashed", "dotted"
    ),
  )

ggsave("HiC_PCA.pdf", width = 2.6, height = 2.0, dpi = 300)

