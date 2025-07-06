##for the feces metabolomics datasets 
cormeta.F<-rcorr(as.matrix(Feces.dataSele1),type = "spearman")
corrplot(cormeta.F$r)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  Cor.matrix<-data.frame( row = rownames(cormat)[row(cormat)[ut]], 
                          column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut] )
}   # 构建相关性对应矩阵
cor.matrixF<-flattenCorrMatrix(cormeta.F$r,cormeta.F$P)

#network visualization  
library(igraph);library(grDevices)
cor.matrixF$Ecolor<-NA
cor.matrixF$Ecolor[which(cor.matrixF$cor<=0)]<-rep("#0099B499",length(which(cor.matrixF$cor<=0)))# "#0099B499"
cor.matrixF$Ecolor[which(cor.matrixF$cor>0)]<-rep("#925E9F99",length(which(cor.matrixF$cor>0)))# "#925E9F99"
cor.matrixF0<-cor.matrixF[which(cor.matrixF$p<0.05),]
cor.matrixF0$cor<-abs(cor.matrixF0$cor)
meta.namesH<-unique(c(cor.matrixF0$row,cor.matrixF0$column))
H.fc<-fc(Heart.data,case = "DM",control = "NS")
H.p<-kruskalTest_p(Heart.data,method = "BH")
meta.colorF<-cbind(H.fc,H.p$p.value)
meta.colorF<-data.frame(rownames(meta.colorF)[which(rownames(meta.colorF) %in% meta.namesH)],
                        meta.colorF[which(rownames(meta.colorF) %in% meta.namesH),])

colnames(meta.colorF)[1:3]<-c("metabolites","FC","pvalue")

meta.colorF$Ncolor<-NA
meta.colorF$Ncolor[which(meta.colorF$FC<0)]<-rep("#00468B99",length(which(meta.colorF$FC<0)))
meta.colorF$Ncolor[which(meta.colorF$FC>0)]<-rep("#ED000099",length(which(meta.colorF$FC>0)))  

meta.colorF<- graph_from_data_frame(cor.matrixF0,directed=FALSE)
## 需要排序
names(V(meta.colorF))

meta.colorF$metabolites<-factor(meta.colorF$metabolites,levels = names(V(meta.colorF)))
meta.colorF<-meta.colorF[order(meta.colorF$metabolites),]

set_vertex_attr(meta.colorF,"name", index = V(meta.colorF), meta.colorF$metabolites)
set_vertex_attr(meta.colorF,"color", index = V(meta.colorF),meta.colorF$Ncolor)
meta.colorF$FC[which(abs(meta.colorF$FC)<0.1)]=0.1
V(meta.colorF)$size<-abs(meta.colorF$FC)*10
V(meta.colorF)$color<-meta.colorF$Ncolor
E(meta.colorF)$width<-cor.matrixF0$cor*5
E(meta.colorF)$color<-cor.matrixF0$Ecolor #E(meta.colorF)$color[E(meta.colorF)$width>0] <- "white"
windows(width=4, height=4)
plot.igraph(meta.colorF,layout=as.matrix(all.layout[,-3]),edge.curved =0.02,edge.arrow.mode ="-",
            rescale =TRUE,asp = 0.8,vertex.frame.color =meta.colorF$Ncolor,vertex.label.color = "black",
            vertex.label.dist =0.3,alpha = 0.3,vertex.label.cex=0.8)
all.layout<-data.frame(layout_with_graphopt(meta.colorF),names(V(meta.colorF)))
all.layout<-data.frame(layout_with_graphopt(meta.colorF),names(V(meta.colorF)))
n1<-which(!(names(V(meta.colorF)) %in% names(V(cor.igrahK))) )
all.layout1<-rbind(all.layout[-n1,]) #相同的变量
all.layout1$names.V.meta.colorF..<-factor(all.layout1$names.V.meta.colorF..,levels =names(V(cor.igrahK))[-n2] )                     #
all.layout1<-all.layout1[order(all.layout1$names.V.meta.colorF..),]
all.layout1.1<-rbind(all.layout1,all.layout[n1,])
all.layout1.1$names.V.meta.colorF..<-as.character(all.layout1.1$names.V.meta.colorF..)
all.layout1.1$names.V.meta.colorF..[29:31]<-c("X3.Methylhistidine","Creatinine","Fumarate")
all.layout1.1$names.V.meta.colorF..<-factor(all.layout1.1$names.V.meta.colorF..,levels = names(V(cor.igrahK)))
all.layout1.1<-all.layout1.1[order(all.layout1.1$names.V.meta.colorF..),]
n2<-which(!(names(V(cor.igrahK)) %in% names(V(meta.colorF))) )
names(V(cor.igrahK))[n2]


