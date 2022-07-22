library(tidyverse)

#-------------16S#----------
ps0 = readRDS("./ps_16s.rds")
# map = read.delim("./map.txt")
# head(map)
# row.names(map) = map$ID
# sample_data(ps0) = map
# saveRDS(ps0,"./ps_16s.rds")
#-主题--颜色等------#-------
source("F:\\Shared_Folder\\Function_local\\R_function\\micro//total_amplicon.R")
#---扩增子环境布置

res = theme_my()
mytheme1 = res[[1]];mytheme2 = res[[2]]; colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
result<- dir.amp(ps0 = ps0)
res1path = result[[1]];res1path
id = result[[2]];id
# 微生物组数据输入#------

#--样本筛选-根据分组和ID等map文件中的信息#-------
# ps_sub <- subset_samples(ps0,!ID %in% c("sample1"));ps_sub 
# 序列筛选-根际七大门类信息或者OTU的ID信息#---------

#--样本筛选-根据分组和ID等map文件中的信息#-------
# ps_sub <- subset_samples(ps0,!ID %in% c("sample1"));ps_sub 
# 序列筛选-根际七大门类信息或者OTU的ID信息#---------
ps0 <- ps0 %>%
  # subset_taxa(
  #   Kingdom == id
  # ) %>% 
  phyloseq::filter_taxa(function(x) sum(x ) > 0 , TRUE);ps0

#--最终确定的phyloseq对象定义为ps#--------
ps = ps0

# ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)
#--提取有多少个分组#-----------
gnum = phyloseq::sample_data(ps)$Group %>%unique() %>% length()
gnum

# 设定排序顺序
axis_order = phyloseq::sample_data(ps)$Group %>%unique()

#--默认进化树展示的微生物数量
Top_micro = 150

# 堆叠柱状图展示TOp前是个门,j 展示的水平为门水平#-----
Top = 12
# rank.names(ps0)
jj = j = "Phylum"

# -类韦恩图的二分网络设置过滤阈值#-------
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

# 差异分析 edger设定分组#------
# group1 = c("Gro1","Gro2")
# b= data.frame(group1)
b = NULL# 如果每两个组之间都做差异，那就指定

# 设置CK，用于双向柱状图绘制#------
CK = unique(phyloseq::sample_data(ps)$Group)[1]

# 热图和气泡图差异系数最大的前多少个OTU#--------
heatnum　=　30


# 用python做lefse#-----
lefsenum = 0
ps_lefse <- ps %>%
  phyloseq::subset_taxa(
    # Kingdom == "Fungi"
    Kingdom == id
    # Genus  == "Genus1"
    # Species %in%c("species1") 
    # row.names(tax_table(ps0))%in%c("OTU1")
  )


ps_lefse = ggClusterNet::filter_OTU_ps(ps = ps_lefse,Top = 400)

# #文件预处理
# format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
# # 注意这里 –c用来指定分组信息-u 1指定样品信息
# #文件分析,这里-l设置LDA阈值，默认为2，我们使用4 会更加严格
# ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res  -l 2
# #柱状图绘制
# plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
# #树状图绘制
# plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
# #做每个差异的柱状图
# mkdir biomarkers_raw_images
# plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/


#--R 语言做lefse法分析-过滤#----------
# ps_sub = filter_taxa(ps0, function(x) sum(x ) > 1000 , TRUE);ps_sub
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)


#---机器学习部分#------
# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
ROC = FALSE
# 是否做交叉检验
rfcv = FALSE
# 选择多少个重要变量
optimal = 40

# 网络
# 过滤多少丰度的做网络
N = 200
zipi = FALSE

#-------功能预测#----


if (is.null(ps0@refseq)) {
  Tax4Fun = FALSE
} else if(!is.null(ps0@refseq)){
  Tax4Fun = TRUE
}



if (Tax4Fun) {
  dir.create("data")
  otu = ps0 %>% 
    ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps0 %>% 
    ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  ps0 = ps0
  
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}

#-------------ITS#----------
ps0 = readRDS("./ITS//ps_ITS.rds")
tax_table(ps0)
# map = read.delim("./map.txt")
# head(map)
# row.names(map) = map$ID
# sample_data(ps0) = map
# saveRDS(ps0,"./ps_ITS.rds")
#-主题--颜色等------#-------
source("F:\\Shared_Folder\\Function_local\\R_function\\micro//total_amplicon.R")
#---扩增子环境布置
package.amp()
result<- dir.amp(ps0)
res1path = result[[1]];res1path
res1path = "./result_and_plot/Base_diversity_ITS/"
id = result[[2]];id

# 微生物组数据输入#------

#--样本筛选-根据分组和ID等map文件中的信息#-------
# ps_sub <- subset_samples(ps0,!ID %in% c("sample1"));ps_sub 
# 序列筛选-根际七大门类信息或者OTU的ID信息#---------
ps0 <- ps0 %>%
  # subset_taxa(
  #   Kingdom == id
  # ) %>% 
  filter_taxa(function(x) sum(x ) > 0 , TRUE);ps0

#--最终确定的phyloseq对象定义为ps#--------
ps = ps0



# ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)
#--提取有多少个分组#-----------
gnum = sample_data(ps)$Group %>%unique() %>% length()
gnum

# 设定排序顺序
axis_order = c("Soil","Leaf", "Waxberry","D","F")
sample_data(ps)$Group %>%unique()

#--默认进化树展示的微生物数量
Top_micro = 150

# 堆叠柱状图展示TOp前是个门,j 展示的水平为门水平#-----
Top = 12
# rank.names(ps0)
jj = j = "Phylum"

# -类韦恩图的二分网络设置过滤阈值#-------
ps_biost = filter_OTU_ps(ps = ps,Top = 500)

# 差异分析 edger设定分组#------
# group1 = c("Gro1","Gro2")
# b= data.frame(group1)
b = NULL# 如果每两个组之间都做差异，那就指定

# 设置CK，用于双向柱状图绘制#------
CK = unique(sample_data(ps)$Group)[1]

# 热图和气泡图差异系数最大的前多少个OTU#--------
heatnum　=　30


# 用python做lefse#-----
lefsenum = 0
ps_lefse <- ps %>%
  subset_taxa(
    # Kingdom == "Fungi"
    Kingdom == id
    # Genus  == "Genus1"
    # Species %in%c("species1") 
    # row.names(tax_table(ps0))%in%c("OTU1")
  )


ps_lefse = filter_OTU_ps(ps = ps_lefse,Top = 400)

# #文件预处理
# format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
# # 注意这里 –c用来指定分组信息-u 1指定样品信息
# #文件分析,这里-l设置LDA阈值，默认为2，我们使用4 会更加严格
# ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res  -l 2
# #柱状图绘制
# plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
# #树状图绘制
# plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
# #做每个差异的柱状图
# mkdir biomarkers_raw_images
# plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/


#--R 语言做lefse法分析-过滤#----------
# ps_sub = filter_taxa(ps0, function(x) sum(x ) > 1000 , TRUE);ps_sub
ps_Rlefse = filter_OTU_ps(ps = ps,Top = 400)


#---机器学习部分#------
# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
ROC = FALSE
# 是否做交叉检验
rfcv = FALSE
# 选择多少个重要变量
optimal = 40

# 网络
# 过滤多少丰度的做网络
N = 200
zipi = FALSE

#-------功能预测#----


if (is.null(ps0@refseq)) {
  Tax4Fun = FALSE
} else if(!is.null(ps0@refseq)){
  Tax4Fun = TRUE
}


if (Tax4Fun) {
  dir.create("data")
  otu = ps0 %>% 
    filter_OTU_ps(1000) %>%
    vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps0 %>% 
    filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  library(Biostrings)
  writeXStringSet(rep,"./data/otu.fa")
  ps0 = ps0
  
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}






