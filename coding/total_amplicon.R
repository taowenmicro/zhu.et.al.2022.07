
# dir.amp(ps0)
# theme.amp(theme = T)
# package.amp()
# color.amp()
#-扩增子分析需要导入的R包
package.amp <- function(){
  #导入R包#-------
  library(phyloseq)
  library(tidyverse)
  library(ggClusterNet)
  library(EasyStat)
  library(fs)
  library(ggthemes)
  library(RColorBrewer)#调色板调用包
  
  
  
  # library(ggplot2)
  # library(dplyr)
  library(magrittr)
  # library(RColorBrewer)
  # devtools::install_github("YuLab-SMU/treeio") # 安装1.7以上版本的才能支持MicrobiotaProcess
  # devtools::install_github("YuLab-SMU/MicrobiotaProcess")
  # BiocManager::install("treeio")
  library(MicrobiotaProcess)
  # library(tibble)
  library(ggsignif)
  library(ggtree)
  library(ggtreeExtra)
  # library(ggplot2)
  library(ggstar)
  library(MicrobiotaProcess)
  library(ggnewscale)
  library(grid)
  
}
package.amp()
#-1 直接在当前文件夹创建一个名为：result_and_plotde 的文件夹
#-2 根据给出的phyloseq对象，使用tax表格得出这是细菌还是真菌数据，然后建立对应的二级文件夹
dir.amp <- function(ps0){
  
  # # 建立结构保存一级目录#--------
  result_path <- paste("./","/result_and_plot/",sep = "")
  dir_create(result_path)
  #---构建结果保存文件夹#---------
  tax.1 = c("Fungi",
    "fungi",
    "K__Fungi",
    "k__Fungi",
    "d__Fungi",
    "d__fungi")
  
  TFdir.f <- as.data.frame(table(
    tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in%
                               tax.1] > 10
  
  if (length(TFdir.f) != 0) {
    print("ITS")
    res1path <- paste(result_path,"/Base_diversity_ITS",sep = "")
    id = tax.1
  }
  
  tax.2 = c("Bacteria",
            "K__Bacteria",
            "k__Bacteria",
            "d__Bacteria",
            "k:Bacteria",
            "bacteria")
  
  
  TFdir.b <- as.data.frame(table(tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in%
                                                            tax.2 ] > 10
  
  if (length(TFdir.b) != 0) {
    print("16s")
    res1path <- paste(result_path,"/Base_diversity_16s",sep = "")
    id = tax.2
  }
  dir.create(res1path)
  
  return(list(res1path,id))
}

#-包括几个常用的主题

    # 设置主题#-----
    mytheme1 = theme_bw() + theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      legend.position="right",
      
      legend.title = element_blank(),
      legend.background=element_blank(),
      legend.key=element_blank(),
      # legend.text= element_text(size=7),
      # text=element_text(),
      # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
    )
    
    mytheme2 = theme_bw() + theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      legend.position="right",
      
      legend.title = element_blank(),
      legend.background=element_blank(),
      legend.key=element_blank(),
      # legend.text= element_text(size=7),
      # text=element_text(),
      # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14,angle = 90),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
    )

#--设定颜色
  gnum <- unique(sample_data(ps)$Group) %>% length()
  gnum 
  
  if (gnum <= 9) {
    #设定颜色#------------
    #调用所有这个包中的调色板
    RColorBrewer::display.brewer.all()
    #提取特定个数的调色板颜色，会出图显示
    # RColorBrewer::display.brewer.pal(9,"Set1")
    colset1 <- c(brewer.pal(9,"Set1"))
    colset2 <- brewer.pal(12,"Paired")
    colset3 <- c(brewer.pal(12,"Set1"),brewer.pal(9,"Pastel1"))
    colset4 = colset3
  }
  
  
  if (gnum > 9) {
    #设定颜色#------------
    #调用所有这个包中的调色板
    RColorBrewer::display.brewer.all()
    #提取特定个数的调色板颜色，会出图显示
    # RColorBrewer::display.brewer.pal(9,"Set1")
    colset1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(gnum)
    colset2 <- brewer.pal(12,"Paired")
    colset3 <- c(brewer.pal(12,"Set1"),brewer.pal(9,"Pastel1"))
    colset4 = colset3
  }
  






