#Este script sirve para hacer la grafica de burbujas 
#Autor Hugo Castelán 
args <- commandArgs(TRUE)
file_in <- as.character(args[1])
#Cambiar a coma o por tabulador depede de la tabla 
virus_dsDNA<-read.table(file_in, sep = "\t", header = T)
library(ggplot2)
library(reshape)
prop<-prop.table(data.matrix(datos), 2)
#Modifca la tabla para ggpplo2
df=melt(virus_dsDNA,	variable_name="Families")
#Estos nombre de las columnas 
colnames(df)<-c("Families","metagenomes","abundance")
df[df==0]<-NA
ggplot(df, aes(x =  metagenomes, y = Families, color= Families))+geom_point(aes(size=abundance))+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
