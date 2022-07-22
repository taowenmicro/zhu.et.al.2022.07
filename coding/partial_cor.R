

# # sample#--------
# library(ggClusterNet)
# 
# data("env1")
# head(env1)
# head(env1)
# res = partial_cor(
#   data = env1,
#   method = "spearman",
#   x = c("env1","env7"),# 偏相关变量
#   y = "env2", # 偏相关变量
#   z = colnames(env1)[3:6])
# 
# res[[1]]

partial_cor = function(
  data = env1,
  method = "spearman",
  x = c("env1","env7"),# 偏相关变量
  y = "env2", # 偏相关变量
  z = colnames(env1)[3:6]
  
){
  a = c()
  b = c()
  c = c()
  R = c()
  P = c()
  
  if (length(x) == 1) {
    for (i in 1:length(z)) {
      # spearman 偏相关系数
      tem1 = ppcor::pcor.test(x = data[[x]], y = data[[y]], z = data[[z[i]]], method = method)
      a[i] = x
      b[i] = y
      c[i] = z[i]
      R[i] = tem1$estimate
      P[i] = tem1$p.value
    }
  }
  
  
  
  
  if (length(x) > 1) {
    for (j in 1:length(x)) {
      for (i in 1:length(z)) {
        # spearman 偏相关系数
        tem1 = ppcor::pcor.test(x = data[[x[j]]], y = data[[y]], z = data[[z[i]]], method = method)
        a[i] = x[j]
        b[i] = y
        c[i] = z[i]
        R[i] = tem1$estimate
        P[i] = tem1$p.value
        pval <- P
        star=ifelse(pval>0.1,"NS",ifelse(pval>0.05,".",ifelse(pval>0.01,"*",ifelse(pval>0.001,"**","***"))))
        star
        tem2 = data.frame(a,b,c,R,P,sig = star)
      }
      
      if (j == 1) {
        tem3 = tem2
      } else {
        tem3 = rbind(tem3,tem2)
      }
      
    }
    
  }
  
  
  
  p <- ggplot(tem3, aes(x = c, y = a)) +
    geom_tile(aes(fill = R))+  #绘制热图
    scale_fill_gradientn(colors = c('blue', 'grey95', 'red'), limit = c(-1, 1)) +  
    geom_text(aes(label = sig), size = 2.8) +  
    facet_grid(~b) +  #分面图
    theme(axis.line = element_blank(), 
          axis.ticks = element_blank(), axis.text = element_text(color = 'black')) +  
    labs(x = '\nPartial correlation control', y = '', fill = "(Partial)\ncorrelations\n\n r") + 
    scale_x_discrete(expand = c(0, 0)) +  
    scale_y_discrete(expand = c(0, 0))
  
  return(list(p,tem3))
  
}

# emf_cor <- Hmisc::rcorr(as.matrix(data), type = method)
# emf_cor$r  
# emf_cor$P  







