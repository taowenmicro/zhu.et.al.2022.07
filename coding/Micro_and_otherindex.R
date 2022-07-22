
env = read.csv("./data/dataNEW/env.csv")


head(env)

# library(ggClusterNet)
# library(phyloseq)
# library(tidyverse)
result_path <- paste("./","/result_and_plot/",sep = "")
dir_create(result_path)

TFdir.f <- as.data.frame(
  table(tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in% c("Fungi",
                                                                                                            "fungi",
                                                                                                            "K:Fungi",
                                                                                                            "k__Fungi",
                                                                                                            "D__Fungi",
                                                                                                            "d__Fungi",
                                                                                                            "D:Fungi",
                                                                                                            "d::Fungi",
                                                                                                            "k__Fungi")] > 10
if (length(TFdir.f) != 0) {
  print("ITS")
  res1path <- paste(result_path,"/Micro_and_other_index_ITS",sep = "")
}
TFdir.b <- as.data.frame(
  table(tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in% c("Bacteria",
                                                                                                            "bacteria",
                                                                                                            "d:Bacteria",
                                                                                                            "D__Bacteria",
                                                                                                            "D:Bacteria",
                                                                                                            "d__Bacteria",
                                                                                                            "K__Bacteria",
                                                                                                            "k__Bacteria")] > 10
if (length(TFdir.b) != 0) {
  print("16s")
  res1path <- paste(result_path,"/Micro_and_other_index_16s_220113",sep = "")
} else if (AM) {
  print("16s")
  res1path <- paste(result_path,"/Micro_and_other_index_AM",sep = "")
}
dir.create(res1path)

#--微生物和环境因子共排序RDA/CCA #---------
#----------附图环境一因子联合#------
library(EasyStat)
RDApath = paste(res1path,"/RDA_CCA/",sep = "")
dir.create(RDApath)

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/rda-cca.R")
envRDA = env
head(env)
row.names(envRDA) = env$ID
envRDA$ID = NULL
head(envRDA)

result = RDA_CCA(ps = ps_Rlefse,
                 env = envRDA,
                 path = RDApath,
                 chose.env = F
                 )
#提取图片
p1 = result[[1]] + mytheme1 + scale_fill_manual(values = colset1)
p1
# 提取作图数据
dataplot = result[[2]]
# 提取带有标记的图片
p2 = result[[3]]+ mytheme1  + scale_fill_manual(values = colset1)
#提取理化提供群落关系的渐检验结果
aov = result[[4]]

# library("Cairo")

##保存
plotnamea = paste(RDApath,"/RDA_envnew.pdf",sep = "")
ggsave(plotnamea, p1, width = 8, height = 6)
plotnamea4 = paste(RDApath,"/RDA_envnew.jpg",sep = "")
ggsave(plotnamea4, p1, width = 8, height = 6)


filenamea = paste(RDApath,"dataplotnew.txt",sep = "")
write.table(dataplot ,file=filenamea,sep="\t",col.names=NA)

filenamea = paste(RDApath,"aovnew.txt",sep = "")
write.table(aov,file=filenamea,sep="\t",col.names=NA)

plotnamea = paste(RDApath,"/RDA_envlabelnew.pdf",sep = "")
ggsave(plotnamea, p2, width = 18, height = 12)#, device = cairo_pdf, family = "Song"
plotnamea4 = paste(RDApath,"/RDA_envlabelnew.png",sep = "")
ggsave(plotnamea4, p2, width = 18, height = 12)# , device = cairo_pdf, family = "Song"


result = RDA_CCA_explain_percent(ps = ps_Rlefse,
                                 env.dat = envRDA)
out = result[[1]]
wxp = result[[2]]

filenamea = paste(RDApath,"each_env_exp_percent.csv",sep = "")
write.csv(out,filenamea)
filenamea = paste(RDApath,"all_index_explain_percent.csv",sep = "")
write.csv(exp,filenamea)



#--微生物群落和环境数据组合图表#---------
compath = paste(res1path,"/Conbine_env_plot/",sep = "")
dir.create(compath)

otu1 = as.data.frame(t(vegan_otu(ps)))
head(otu1)

tabOTU1 = list(bac = otu1)
# tabOTU1 = list(b = otu)
# tabOTU1 = list(f = otu2)

p0 <- MatCorPlot(env.dat = envRDA,
                 tabOTU = tabOTU1,
                 diag = F,
                 range = 0.1,
                 numpoint = 21,
                 sig = TRUE,
                 siglabel = FALSE,
                 shownum = F,
                 numsymbol = NULL,
                 lacx = "right",
                 lacy = "bottom",
                 p.thur = 0.05,
                 onlysig = T
)
p0

p0 <- p0 + scale_colour_manual(values = c("blue","red")) + 
  scale_fill_distiller(palette="PRGn")
p0

FileName <- paste(compath,"Conbine_envplot", ".pdf", sep = "")
ggsave(FileName, p0,width = 15,height = 10)# , device = cairo_pdf, family = "Song"

FileName <- paste(compath,"Conbine_envplot", ".png", sep = "")
ggsave(FileName, p0,width = 15,height = 10)# , device = cairo_pdf, family = "Song"



#--拆分开来
#--- mantel test
rep = MetalTast (env.dat = envRDA, tabOTU = tabOTU1,distance = "bray",method = "metal")
repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
repP = rep[seq(from=1,to=dim(rep)[2],by=2)]
head( repR)
head( repP)
mat = cbind(repR,repP)
head(mat)

FileName <- paste(compath,"Conbine_envplot_data", ".csv", sep = "")
write.csv(mat,FileName)


result <- Miccorplot(data = envRDA,
                     method.cor = "spearman",
                     cor.p = 0.05,
                     x = F,
                     y = F,
                     diag = T,
                     lacx = "left",
                     lacy = "bottom",
                     sig = T,
                     siglabel = F,
                     shownum = F,
                     numpoint = 21,
                     numsymbol = NULL
)

p1 <- result[[1]]
p1 <- p1 + scale_fill_distiller(palette="PRGn")
p1
FileName <- paste(compath,"envCorplot", ".pdf", sep = "")
ggsave(FileName, p1,width = 15,height = 10)
FileName <- paste(compath,"envCorplot", ".png", sep = "")
ggsave(FileName, p1,width = 15,height = 10)


#--在全部其他指标中影响微生物群落的最重要的指标#---------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/ramdom_env.plot.R")

rampath = paste(res1path,"/Random_env/",sep = "")
dir.create(rampath)


result <- ramdom_env.plot(ps  = ps,env = envRDA )
p <- result[[2]]
p
data = result[[1]]
head(data)
hit <- dim(envRDA)[2]
hit
FileName <- paste(rampath,"ranImportant", ".pdf", sep = "")
ggsave(FileName, p,width = 6,height =hit/5)


FileName <- paste(rampath,"ranImportant", ".csv", sep = "")
write.csv(data,FileName)

# #--在全部其他指标中影响微生物群落的最重要的指标#---------
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro/ramdom_env_micro.heatplot.R")
# rampath = paste(res1path,"/Random_env/",sep = "")
# dir.create(rampath)
# # ps = ps0
# map = sample_data(ps)
# result <- ramdom_env_micro.heatplot( ps = ps, env = env,seed = 1)
# p <- result[[1]] + mytheme2
# dat = result[[2]]
# wid = dim(env)[2]
# hei = length(unique(map$Group))
# hei
# FileName <- paste(rampath,"Randomforest_env_micro_heatmap", ".pdf", sep = "")
# ggsave(FileName, p,width = 5 * hei ,height = wid/4)
# 
# FileName <- paste(rampath,"Randomforest_env_micro_heatmap", ".jpg", sep = "")
# ggsave(FileName, p,width = 5 * hei ,height = wid/4)
# 
# FileName <- paste(rampath,"Random_env_micro_heatmap", ".csv", sep = "")
# write.csv(dat,FileName)




#---特征微生物和环境因子的关系探索#----

#---不同分类环境因子相关图表
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/cor_env_ggcorplot.R")
heatpath = paste(res1path,"/cor_env_heapmap_boplot/",sep = "")
dir.create(heatpath)
library(sna)
library(igraph)

#---提取门水平
jj = 6
tran = TRUE
Top = 20
psdata <- tax_glom_wt(ps = ps,ranks = jj)
if (tran == TRUE) {
  psdata = psdata%>%
    transform_sample_counts(function(x) {x/sum(x)} )
}


otu = otu_table(psdata)
tax = tax_table(psdata)

if (dim(otu)[1] < Top) {
  top10 <- otu[names(sort(rowSums(otu), decreasing = TRUE)[1:dim(otu)[1]]),]
  top10 = t(top10)
} else {
  top10 <- otu[names(sort(rowSums(otu), decreasing = TRUE)[1:Top]),]
  top10 = t(top10)
}


head(top10)



result = cor_env_ggcorplot(
  env1 = envRDA,
  env2 = top10,
  label =  T,
  col_cluster = T,
  row_cluster = T,
  method = "spearman",
  r.threshold= 0,
  p.threshold= 0
)



p1 <- result[[1]] 
p1
p2 <- result[[2]]
p2

hei = dim(env)[2]/5
wid = Top

filename = paste(heatpath,"/",jj,"ggheatmap.pdf",sep = "")
ggsave(filename,p1,width = Top/1.3,height = dim(env)[2]/5)
filename = paste(heatpath,"/",jj,"ggbubble.pdf",sep = "")
ggsave(filename,p2,width = Top/1.3,height = dim(env)[2]/5)

filename = paste(heatpath,"/",jj,"ggheatmap.jpg",sep = "")
ggsave(filename,p1,width = Top/1.3,height = dim(env)[2]/5)
filename = paste(heatpath,"/",jj,"ggbubble.jpg",sep = "")
ggsave(filename,p2,width = Top/1.3,height = dim(env)[2]/5)


#--微生物和环境因子的网络#--------

#--微生物和环境的网络
library(phyloseq)
library(dplyr)
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(ggplot2)
library(ggClusterNet)


Envnetplot<- paste(res1path,"/Env_network",sep = "")

dir.create(Envnetplot)

ps16s = ps0
psITS = NULL

# source("E:\\Shared_Folder\\Function_local\\R_function\\my_R_packages\\ggClusterNet\\R/merge16S_ITS.R")

ps.merge <- merge16S_ITS(ps16s = ps16s,
                         psITS = NULL,
                         N16s = 150)
ps.merge

map = sample_data(ps.merge)
map$Group = "one"
sample_data(ps.merge) <- map

#--环境因子导入

#--导入多个环境变量及其分组
data1  = envRDA
head(data1)
data1$ID = row.names(data1)

data1 <- as.tibble(data1) %>% dplyr::select(ID,everything()) %>% as.data.frame()

colnames(data1)[1]=  "ID"
head(data1)

Gru = data.frame(ID = colnames(data1)[-1],group = "env" )
ps.merge


source("E:/Shared_Folder/Function_local/R_function/micro//corBionetwork2.R")




library(ggpubr)
library(ggClusterNet)
library(ggrepel)
library(dplyr)
result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        r.threshold = 0.6, # 相关阈值
                        p.threshold = 0.05,
                        group = "Group",
                        env = data1, # 环境指标表格
                        envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        path = Envnetplot,# 结果文件存储路径
                        fill = "Phylum", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = TRUE, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = F, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 18,
                        label = TRUE,
                        height = 10
)


p = result[[1]]
p
# 全部样本网络参数比对
data = result[[2]]
plotname1 = paste(Envnetplot,"/network_all.jpg",sep = "")
ggsave(plotname1, p,width = 15,height = 12)
plotname1 = paste(Envnetplot,"/network_all.png",sep = "")
ggsave(plotname1, p,width = 10,height = 8)
plotname1 = paste(Envnetplot,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 15,height = 12)
tablename <- paste(Envnetplot,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)



#---细菌真菌网络#-------------

#--微生物和环境因子的网络#--------

#--微生物和环境的网络
library(phyloseq)
library(dplyr)
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(ggplot2)
library(ggClusterNet)


Envnetplot<- paste(res1path,"/16s_ITS_network",sep = "")

dir.create(Envnetplot)

ps16s = ps
psITS = ps

# source("E:\\Shared_Folder\\Function_local\\R_function\\my_R_packages\\ggClusterNet\\R/merge16S_ITS.R")

ps.merge <- merge16S_ITS(ps16s = ps16s,
                         psITS = psITS,
                         N16s = 100,
                         NITS = 100
                         )
ps.merge

tax = ps.merge %>% vegan_tax() %>%
  as.data.frame()

result <- corBiostripe(ps = ps.merge, r.threshold = 0.8,
                       p.threshold = 0.05
                       )
#--提取相关矩阵
cor = result[[1]]


#--分组信息提取
tax = as.data.frame((vegan_tax(ps.merge)))
netClu = data.frame(ID = row.names(tax),group = tax$filed)
head(netClu)

result2 = PolygonClusterG(cor = cor,nodeGroup =netClu,zoom = 1.5,zoom2 = 1 )
node = result2[[1]]
head(node)
nodes <- node %>%
  inner_join(netClu ,by =c("elements" = "ID") )
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
### 出图
head(edge)
head(nodes)
library(ggrepel)
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,
                                    color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group),pch = 21, size = 5,data = nodes) +
  geom_text_repel(aes(X1, X2,label = elements),data = nodes,size = 1.5) +
  scale_colour_brewer(palette = "Set1") + theme_void()
pnet

p = pnet

# 全部样本网络参数比对
data = edge
plotname1 = paste(Envnetplot,"/network_all.jpg",sep = "")
ggsave(plotname1, p,width = 6,height = 10)
plotname1 = paste(Envnetplot,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 6,height = 10)

tablename <- paste(Envnetplot,"/edge_",".csv",sep = "")
write.csv(data,tablename)

data = nodes
tablename <- paste(Envnetplot,"/nodes_",".csv",sep = "")
write.csv(data,tablename)




