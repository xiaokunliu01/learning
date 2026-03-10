# 1. 生成模拟数据（包含多个相近极值）
set.seed(456)
x <- 1:100
y <- sin(x/5) + rnorm(100, sd = 0.3)  # 加入噪声

# # 使用 diff() 和 sign() 找局部极值
# 改进的寻找波峰函数
find_peaks <- function(y, threshold = 0) {
  dy <- diff(y)
  peaks <- which(diff(sign(dy)) == -2) + 1
  
  # 只有当波峰的高度大于 threshold 时才保留
  peaks <- peaks[y[peaks] > threshold]
  
  return(peaks)
}

# 同样改进波谷函数
find_valleys <- function(y, threshold = 0) {
  dy <- diff(y)
  valleys <- which(diff(sign(dy)) == 2) + 1
  
  # 只有当波谷的高度小于 threshold 时才保留
  valleys <- valleys[y[valleys] < threshold]
  
  return(valleys)
}


# 应用
# x <- df[,1]
# y <- df[,2]
peaks <- find_peaks(y, threshold = 0)
valleys <- find_valleys(y, threshold = 0)

# 绘图
plot(x, y, type = "l", main = "Local Extrema")
points(x[peaks], y[peaks], col = "red", pch = 19, cex = 1.5)
points(x[valleys], y[valleys], col = "blue", pch = 19, cex = 1.5)