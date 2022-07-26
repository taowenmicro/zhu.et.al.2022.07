library(picante)
library(ape)
library(vegan)
library(FSA)
library(eulerr)
library(grid)
library(gridExtra)
require(minpack.lm)
require(Hmisc)
require(stats4)
library(parallel)


psphy = filter_taxa(ps, function(x) sum(x ) > 1000 , TRUE);psphy
package.amp()

# TFdir.f <- as.data.frame(table(tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in% c("Fungi","K__Fungi","k__Fungi")] > 10
# if (length(TFdir.f) != 0) {
#   print("ITS")
#   res1path <- paste(result_path,"/ITS_env_phylo_processing",sep = "")
# }
# TFdir.b <- as.data.frame(table(tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in% c("Bacteria","d__Bacteria","K__Bacteria","k__Bacteria")] > 10
# if (length(TFdir.b) != 0) {
#   print("16s")
#   res1path <- paste(result_path,"/16S_env_phylo_processing",sep = "")
# }
# dir.create(res1path)

# res1path = "./"


phypath = paste(res1path,"/Phylogenetic_analyse_spacies/",sep = "")
dir.create(phypath)

map = phyloseq::sample_data(ps)
n = map$Group %>% unique() %>%
  length()
n


#--中性模型#-----
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\neutralModel.R")
result = neutralModel(ps = ps,group  = "Group",ncol = 5)
#--合并图表
p1 =  result[[1]]
p1


FileName <- paste(phypath,"1_neutral_modelCul", ".pdf", sep = "")
ggsave(FileName, p1,width = 16,height = 4)
FileName <- paste(phypath,"1_neutral_modelCul", ".png", sep = "")
ggsave(FileName, p1,width = 16,height = 4)


#--系统发育信号#---
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\phyloSignal_and_phySigplot.R")
envphy = env
row.names(env) = env$ID
env$ID = NULL
# #这个函数没有输出只会保存到指定的路径中，因为这个计算很浪费时间#---------
phypath2 = paste(phypath,"/phyloSignal/",sep = "")
dir.create(phypath)
phyloSignal(ps = psphy,
            group  = "Group",
            env = env,
            path = phypath2)

result = phySigPlot(ps = ps,group  = "Group",env = env,path = phypath2)
#
#提取图片
p2 = result[[1]] + mytheme1
p2
#-提取作图数据
data = result[[2]]
head(data)

FileName <- paste(phypath,"2_phySigPlot", ".pdf", sep = "")
ggsave(FileName, p2,width = 15,height = 6)
FileName <- paste(phypath,"2_phySigPlot", ".csv", sep = "")
write.csv(data,FileName)



#---计算零模型#-------
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\nullModel1.R")

result <- nullModel(ps = psphy,
                    group="Group",
                    dist.method =  "bray",
                    gamma.method = "total",
                    transfer = "none",
                    null.model = "ecosphere"
                    )

#--分组零模型运行结果
nullModeltab <- result[[1]]

# 比例
ratiotab <- result[[2]]
#-统计量统计差异
aovtab <- result[[3]]

FileName <- paste(phypath,"3_nullModeltab", ".csv", sep = "")
write.csv(nullModeltab,FileName)

FileName <- paste(phypath,"3_ratiotab", ".csv", sep = "")
write.csv(ratiotab,FileName)

# FileName <- paste(phypath,"3_aovtab", ".csv", sep = "")
# write.csv(aovtab,FileName)




#--BNTI#----
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\bNTICul.R")

result = bNTICul(ps = psphy,group  = "Group",num = 10,thread = 1)

bNTI = result[[1]]
head(bNTI)


filename = paste(phypath,"/4_bNTI.csv",sep = "")
write.csv(bNTI, filename)


#--计算RCbray#-----------
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\RCbary.R")

result = RCbary(ps = psphy ,group  = "Group",num = 10,thread = 1)

RCbary = result[[1]]
head(RCbary)
filename = paste(phypath,"/5_RCb.csv",sep = "")
write.csv(RCbary,filename)

#--BetaNTI和RCbray联合出图#---------
# phypath = "./result_and_plot/16S_env_phylo_processing/Phylogenetic_analyse_spacies/"
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\bNTIRCPlot.R")

bNTI = read.csv(paste(phypath,"/4_bNTI.csv",sep = ""),row.names = 1)
head(bNTI)
# RCbray 数据读入，修改列名
RCb = read.csv(paste(phypath,"/5_RCb.csv",sep = ""),row.names = 1) %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)
head(RCb)

result = bNTIRCPlot(ps = psphy ,RCb  = RCb,bNTI = bNTI,group  = "Group")

#--bNTI出图片
p3 <- result[[1]] + mytheme1
p3

#RCbary可视化
p4 <- result[[2]] + mytheme1
p4
#组合图片BNTI，RCbray
p5 <- result[[3]]
p5
plotdata = result[[4]]
head(plotdata)




filename = paste(phypath,"/6_bNTI_RCbray.csv",sep = "")
write.csv(plotdata,filename)


plotdata = result[[5]]
head(plotdata)

filename = paste(phypath,"/6_RCbray_percent.csv",sep = "")
write.csv(plotdata,filename)


FileName <- paste(phypath,"6_bNTI", ".pdf", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".pdf", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".pdf", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

FileName <- paste(phypath,"6_bNTI", ".png", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".png", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".png", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

#---环境因子和BetaNTI相关#---------
source("F:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\EnvCorbNTI.R")

#-导入bNTI函数
bNTIRC = read.csv(paste(phypath,"/6_bNTI_RCbray.csv",sep = ""),row.names = 1)
head(bNTIRC)

map = sample_data(psphy)
head(map)
plot = EnvCorbNTI(ps = psphy,
                  bNTIRC = bNTIRC,
                  group  = "Group",
                  env = envRDA
                  )

## 提取相关分析结果，总图
p6 <- plot[[1]]
#提取单个
# plot[[2]][1]

FileName <- paste(phypath,"7_env_corWithBNTI", ".pdf", sep = "")
ggsave(FileName, p6,width = 16,height = 16)

FileName <- paste(phypath,"7_env_corWithBNTI", ".png", sep = "")
ggsave(FileName, p6,width = 16,height = 16)

