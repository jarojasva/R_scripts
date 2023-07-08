getwd()
setwd()
library(svglite)
library(gplots)
library(RColorBrewer)
data<-read.csv("pyani_ml10562.csv", header=T, sep=",", row.names=1)
head(data)
coul <- colorRampPalette(c("forestgreen","white", "dodgerblue4"))(15)
svglite("ANI.svg", width=10, height=10)
heatmap.2(as.matrix(data),main="ANIm [%]",notecol="black",density.info="none",trace="none",margins=c(18,20),cexRow=0.7,cexCol=0.7,col=coul)
dev.off()



# Cuando hay valores de cero en la tabla de ANI
getwd()
setwd("D:/UNAM/Semestre_9/Art_Genomas_LP/PyANI_Mediterrane/")
library(svglite)
library(gplots)
library(RColorBrewer)
data<-read.csv("pyani_mediterraneo2.csv", header=T, sep=",", row.names=1)
coul <- colorRampPalette(c("white", "white", "white","white","white","white","white","dodgerblue4"))(40)
svglite("ANI.svg", width=10, height=10)
heatmap.2(as.matrix(data),main="ANI [%]",notecol="black",density.info="none",trace="none",margins=c(18,20),cexRow=0.7,cexCol=0.7,col=coul)
dev.off()

#Cuando no se quiere clusterizar
heatmap.2(as.matrix(data),main="ANIm [%]",notecol="black",density.info="none",trace="none",margins=c(18,20),cexRow=0.7,cexCol=0.7,col=coul,Rowv=FALSE, Colv=FALSE)



# COLOCANDO VALORES DENTRO DEL HEATMAP
library(pheatmap)
library(RColorBrewer)
library(svglite)
data<-read.csv("pyani_mediterraneo.csv", header=T, sep=",", row.names=1)
coul <- colorRampPalette(c("white", "gray97", "gray97","gray97","gray97", "steelblue4"))(40)
svglite("ANI_values.svg", width=20, height=10)
pheatmap(data, main = "ANI [%]", fontsize_row = 14, fontsize_col = 14, display_numbers = TRUE, 
	number_color = "black", fontsize_number = 11, color = coul,
	legend_breaks = c(0, 80, 90, 100), border_color = "black")
dev.off()
# number_color <- ifelse(data == 96.00, "white", "black")
# pheatmap(data, main = "ANI [%]", fontsize_row = 12, fontsize_col = 12, display_numbers = TRUE, number_color = number_color, fontsize_number = 12, color = coul)




#Crear un heatmap con gradientes
library(ggplot2)
data<-read.csv("pyani_mediterraneo2.csv", header=T, sep=",", row.names=1)
matriz <- as.matrix(datos)
svglite("ANI.svg", width=10, height=10)
ggplot(data = melt(matriz), aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "lightsteelblue", limits = c(0, 96)) + geom_text(aes(label = value), color = "black") + theme_minimal() + labs(x = "Columnas", y = "Filas", fill = "Valores") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14), axis.text.y = element_text(size = 14))
dev.off()

#con colores intermedios :D
library(ggplot2)
library(scales)
library(reshape)
data<-read.csv("pyani_mediterraneo2.csv", header=T, sep=",", row.names=1)
matriz <- as.matrix(data)
svglite("ANI2.svg", width=20, height=10)
ggplot(data = melt(matriz), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("steelblue", "white", "red"),
                       values = rescale(c(0, 80, 96, 100), to = c(0, 1))) +
  geom_text(aes(label = value), color = "black") +
  theme_minimal() +
  labs(x = "Columnas", y = "Filas", fill = "Valores") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
dev.off()
