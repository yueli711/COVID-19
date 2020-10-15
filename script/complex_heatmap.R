library(ComplexHeatmap)
library(circlize)
#cell56716_PBMC complex_heatmap
setwd("~/complex_heatmap")
heat <- read.table(file = 'combat_cell56716_PBMC.txt', header = T)
heat_matrix <- as.matrix(heat)
#pdf('2_heat_map-2.pdf')
#Heatmap(heat_matrix, name = "15 pct = 0.9 n_sig:150", col = colorRamp2(c(0, 10, 50), c("blue", "white", "red")), show_row_names = TRUE) 
#dev.new()
col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  
show_row_names = TRUE 
column_title = "Column title"
column_title_side = "top"  
column_title_gp = gpar(fontsize = 6, fontface = "bold") 

Heatmap(heat_matrix,row_title = "Row title", row_title_gp = gpar(fontsize = 6, fontface = "bold")) 
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
col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  
show_row_names = TRUE 
column_title = "Column title" 
column_title_side = "top"  
column_title_gp = gpar(fontsize = 6, fontface = "bold")  

Heatmap(heat_matrix,row_title = "Row title", row_title_gp = gpar(fontsize = 6, fontface = "bold")) 
row_title_gp = gpar(fontsize = 6, fontface = "bold")

Heatmap(heat_matrix, name = "mtcars", column_title = "mock and infected with SARS-CoV-2 samples", 
        column_title_gp = gpar(fontsize = 16, fontface = "bold"), 
        row_title = "Gene Name", row_title_gp = gpar(fontsize = 16, fontface = "bold"),column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 4))


