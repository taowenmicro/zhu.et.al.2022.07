
# 功能预测--Tax4Fun2#-----

#安装 tax4fun2
# install.packages(pkgs = 'E:/Shared_Folder/Database//Tax4Fun2/Tax4Fun2_1.1.5.tar.gz', repos = NULL, source = TRUE)
# buildReferenceData(path_to_working_directory = './Tax4Fun2/', use_force = FALSE, install_suggested_packages = TRUE)
# buildDependencies(path_to_reference_data = './Tax4Fun2/Tax4Fun2_ReferenceData_v2',
#                   use_force = T, install_suggested_packages = TRUE)

funcpath = paste(otupath,"/Tax4Fun2/",sep = "")
dir.create(funcpath)

path_to_reference_data = "E:/Shared_Folder/Database/Tax4Fun2/Tax4Fun2_ReferenceData_v2"
otudir = uncpath
#加载
library(Tax4Fun2)
#物种注释
#指定 OTU 代表序列、Tax4Fun2 库的位置、参考数据库版本、序列比对（blastn）线程数等
runRefBlast(path_to_otus = './data/otu.fa', 
            path_to_reference_data = path_to_reference_data, 
            path_to_temp_folder = otudir, database_mode = 'Ref100NR', 
            use_force = TRUE, num_threads = 4)

#预测群落功能
#指定 OTU 丰度表、Tax4Fun2 库的位置、参考数据库版本、上步的物种注释结果路径等
makeFunctionalPrediction(path_to_otu_table = './data/otu.txt',
                         path_to_reference_data = path_to_reference_data, 
                         path_to_temp_folder = otudir, database_mode = 'Ref100NR', 
                         normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)

