library(grid)
library(ggplot2)
library(dplyr)

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

df1 <- read.csv("./data/HepG2_Smartseq2_TF_RNA1000_RBP_RNA1000.csv")
df1$Reg_Type <- factor(df1$Reg_Type,levels = c("TF","RBP"))
df1$Method <- factor(df1$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df1$Dataset <- "HepG2_Smartseq2"

df2 <- read.csv("./data/HepG2_DNBelab_TF_RNA1000_RBP_RNA1000.csv")
df2$Reg_Type <- factor(df2$Reg_Type,levels = c("TF","RBP"))
df2$Method <- factor(df2$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df2$Dataset <- "HepG2_10x"

df3 <- read.csv("./data/HepG2_9CL_SCAN_TF_RNA1000_RBP_RNA1000.csv")
df3$Reg_Type <- factor(df3$Reg_Type,levels = c("TF","RBP"))
df3$Method <- factor(df3$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df3$Dataset <- "HepG2_SCANseq2_9CL"

df4 <- read.csv("./data/K562_CELseq_TF_RNA1000_RBP_RNA1000.csv")
df4$Reg_Type <- factor(df4$Reg_Type,levels = c("TF","RBP"))
df4$Method <- factor(df4$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df4$Dataset <- "K562_CELseq"

df5 <- read.csv("./data/K562_STORMseq_TF_RNA1000_RBP_RNA1000.csv")
df5$Reg_Type <- factor(df5$Reg_Type,levels = c("TF","RBP"))
df5$Method <- factor(df5$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df5$Dataset <- "K562_STORMseq"

df6 <- read.csv("./data/K562_Smartseq3_TF_RNA1000_RBP_RNA1000.csv")
df6$Reg_Type <- factor(df6$Reg_Type,levels = c("TF","RBP"))
df6$Method <- factor(df6$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df6$Dataset <- "K562_Smartseq3"

df7 <- read.csv("./data/K562_9CL_SCAN_TF_RNA1000_RBP_RNA1000.csv")
df7$Reg_Type <- factor(df7$Reg_Type,levels = c("TF","RBP"))
df7$Method <- factor(df7$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df7$Dataset <- "K562_SCANseq2_9CL"

df8 <- read.csv("./data/K562_UMI200_SCAN_TF_RNA1000_RBP_RNA1000.csv")
df8$Reg_Type <- factor(df8$Reg_Type,levels = c("TF","RBP"))
df8$Method <- factor(df8$Method,levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM_CT","ARACNe_FDR"))
df8$Dataset <- "K562_SCANseq2_UMI200"

df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- rbind(df,df4)
df <- rbind(df,df5)
df <- rbind(df,df6)
df <- rbind(df,df7)
df <- rbind(df,df8)

write.csv(df,"TF_RBP_1000.csv")
df$Dataset <- factor(df$Dataset,levels = c("HepG2_Smartseq2","HepG2_DNBelab","HepG2_SCANseq2_9CL","K562_CELseq","K562_STORMseq","K562_Smartseq3",
                                           "K562_SCANseq2_9CL","K562_SCANseq2_UMI200"))

df$Method <- factor(df$Method,
                    levels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM","ARACNe_FDR"),
                    labels = c("PIDC", "GRNBOOST2","SINCERITIES","TENET","TENET_A","TENET_B","DeePSEM","ARACNe"))

p <- ggplot() + geom_col(data = df, aes(x = Method, y = EPR, fill = Reg_Type), position = "dodge") +
  scale_fill_manual(
    values = c("TF" = "#a1e9f0", 
               "RBP"    = "#d9b1f0"))+
  facet_grid(Dataset ~ Method, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_blank())+#   element_text(angle = 45, hjust=1)) + 
  geom_hline(yintercept = 1,linetype='dashed')
p
ggsave("barplot_TF_1000.pdf",width = 12,height = 14)