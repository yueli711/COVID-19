library(ComplexHeatmap)
library(circlize)
#cell56716_PBMC complex_heatmap
setwd("/home/li/covid19/result01/complex_heatmap")
heat <- read.table(file = 'combat_cell56716_PBMC.txt', header = T)
heat_matrix <- as.matrix(heat)
#pdf('2_heat_map-2.pdf')
#Heatmap(heat_matrix, name = "15 pct = 0.9 n_sig:150", col = colorRamp2(c(0, 10, 50), c("blue", "white", "red")), show_row_names = TRUE) #绘图
#dev.new()
col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  #设置绘图的颜色，可以根据数据实际情况调整
show_row_names = TRUE #是否显示每行名称，默认是TRUE。如果数据很多，行名称会堆叠显示，而且字非常小，感觉没必要
column_title = "Column title" #  显示列标题
column_title_side = "top"  #设置列标题的位置，可选"top"或"bottom"
column_title_gp = gpar(fontsize = 6, fontface = "bold")  #更改列文本的字体

Heatmap(heat_matrix,row_title = "Row title", row_title_gp = gpar(fontsize = 6, fontface = "bold")) #设置行标题和字体格式
row_title_gp = gpar(fontsize = 6, fontface = "bold")

Heatmap(heat_matrix, name = "mtcars", column_title = "mock and infected with SARS-CoV-2 samples", 
        column_title_gp = gpar(fontsize = 16, fontface = "bold"), 
        row_title = "Gene Name", row_title_gp = gpar(fontsize = 16, fontface = "bold"),column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 4))

#cell56716_BALF complex_heatmap
heat <- read.table(file = 'combat_cell56716_BALF.txt', header = T)
heat_matrix <- as.matrix(heat)
#pdf('2_heat_map-2.pdf')
#Heatmap(heat_matrix, name = "15 pct = 0.9 n_sig:150", col = colorRamp2(c(0, 10, 50), c("blue", "white", "red")), show_row_names = TRUE) #绘图
#dev.new()
col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  #设置绘图的颜色，可以根据数据实际情况调整
show_row_names = TRUE #是否显示每行名称，默认是TRUE。如果数据很多，行名称会堆叠显示，而且字非常小，感觉没必要
column_title = "Column title" #  显示列标题
column_title_side = "top"  #设置列标题的位置，可选"top"或"bottom"
column_title_gp = gpar(fontsize = 6, fontface = "bold")  #更改列文本的字体

Heatmap(heat_matrix,row_title = "Row title", row_title_gp = gpar(fontsize = 6, fontface = "bold")) #设置行标题和字体格式
row_title_gp = gpar(fontsize = 6, fontface = "bold")

Heatmap(heat_matrix, name = "mtcars", column_title = "mock and infected with SARS-CoV-2 samples", 
        column_title_gp = gpar(fontsize = 16, fontface = "bold"), 
        row_title = "Gene Name", row_title_gp = gpar(fontsize = 16, fontface = "bold"),column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 4))


