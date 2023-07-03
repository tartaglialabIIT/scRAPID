# This scripts prepares the data for running ARACNe-AP as a tsv table
library(data.table)

# Set your working directory
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

datasets <- c("HepG2_Smartseq2_RBP_RNA500processed_HepG2_Smartseq2_RBP_RNA500",
              "HepG2_9CL_SCAN_RBP_RNA500processed_HepG2_9CL_SCAN_RBP_RNA500",
              "K562_CELseq_RBP_RNA500processed_K562_CELseq_RBP_RNA500",
              "K562_STORMseq_RBP_RNA500processed_K562_STORMseq_RBP_RNA500",
              "K562_Smartseq3_RBP_RNA500processed_K562_Smartseq3_RBP_RNA500",
              "K562_9CL_SCAN_TF_RNA500processed_K562_9CL_SCAN_TF_RNA500",
              "K562_9CL_SCAN_RBP_RNA500processed_K562_9CL_SCAN_RBP_RNA500",
              "K562_UMI200_SCAN_RBP_RNA500processed_K562_UMI200_SCAN_RBP_RNA500")

identities <- c("HepG2_Smartseq2_RBP_RNA500",
                "HepG2_9CL_SCAN_RBP_RNA500",
                "K562_CELseq_RBP_RNA500",
                "K562_STORMseq_RBP_RNA500",
                "K562_Smartseq3_RBP_RNA500",
                "K562_9CL_SCAN_RBP_RNA500",
                "K562_UMI200_RBP_RNA500")

i <- 1

dir.create("input_tables")

for(d in datasets){
  fname <- paste("./rds/",d,sep = "")
  fname <- paste(fname, ".rds", sep = "")
  fname2 <- paste("./input_tables/ARACNE-RUN-",identities[i],sep = "")
  fname2 <- paste(fname2,"_ARACNe-table.tsv",sep = "")
  mydata <- readRDS(fname)
  data_to_write_out <- as.data.frame(as.matrix(GetAssayData(mydata,slot='counts')))
  fwrite(x = data_to_write_out, row.names = TRUE, sep="\t", file  = fname2)
  i <- i + 1
}

# For the HepG2 DNBelab datasets we computed metacells, so we just need the gene names 
# from the gene selection step to create the table for running ARACNe
HepG2DNBelab <- readRDS("../HepG2_DNBelab_METACELLS/mats_pisces-HepG2_DNBelab.rds")

gnames.RBPRNA500 <- read.csv("../HepG2_DNBelab_METACELLS/gnamesHepG2_DNBelab_RBP_RNA500.txt",
                             sep = '\t',header = FALSE)
HepG2DNBelab.RBPRNA500 <- HepG2DNBelab[unlist(gnames.RBPRNA500),]


data_to_write_out2 <- as.data.frame(as.matrix(HepG2DNBelab.RBPRNA500))
fwrite(x = data_to_write_out2, row.names = TRUE, sep="\t", file ="./input_tables/ARACNE-RUN-HepG2_DNBelab_RBP_RNA500_ARACNe-table.tsv")