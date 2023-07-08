


#Este script sirve para hacer la grafica de burbujas 
#Autor Hugo Castelán 
args <- commandArgs(TRUE)
file_in <- as.character(args[1])
library(ggplot2)
grafica<-read.table(file=file_in, sep="\t", header=T, row.names=1)
require(ggplot2)
require(reshape2)
prop<-prop.table(data.matrix(grafica), 2)
dat_m <- melt(prop)
colnames(dat_m)<-c("Families", "Metagenome", "Abundance")
 ggplot(dat_m, aes(Metagenome, Abundance, fill =Families)) + 
   geom_bar(stat = "identity")+xlab("Sitios")+ 
   ylab("Relative abundance")+theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_manual(values = c("darkslategrey", "darkseagreen4", "seagreen3", 
                                "darkolivegreen1","forestgreen","mediumspringgreen",
                                "limegreen","magenta3","maroon4","purple","orchid4",
                                "palevioletred4","salmon","peru","sienna1","sienna3",
                                "tan3","tomato","skyblue","turquoise","turquoise3","skyblue4"
                                ,"slateblue4","steelblue","steelblue","deepskyblue3",
                                "dodgerblue4","cornflowerblue","antiquewhite","antiquewhite4"
                                ,"azure3","azure4","gainsboro","gray","gray56","honeydew","honeydew3"
                                ,"lightcyan3","lightcyan4","lightgoldenrod4","lightgoldenrod","khaki"
                                ,"khaki4","gold3","goldenrod1","gold","lightgoldenrod4","midnightblue"
                                ,"navyblue","royalblue","skyblue1","steelblue4","lightcyan1","cadetblue"
                                ,"cadetblue3","aquamarine4","aquamarine"))

   