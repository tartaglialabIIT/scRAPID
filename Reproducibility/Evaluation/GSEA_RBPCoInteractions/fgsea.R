library(fgsea)
library(data.table)
library(ggplot2)

base_folder <- '/mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/COINTER_RBPs/'

datasets<- c('HEK293T')
strings <- c('eCLIP_HepG2_single_FC_1_P_3_canonical')
sets <- c('RBP_RNA1000','RBP_RNA2000','RBP_RNA3000')

methods <- c('PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR')

gmt.file <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_LAST/Bioplex_HEK293T.gmt"
#gmt.file <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_RBPs/Lang.gmt"

pathways <- gmtPathways(gmt.file)

wd <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_LAST/HEK293T_Jaccard/"
setwd(wd)

for(i in 1:length(datasets)){
  print(datasets[i])
  res_folder <- paste(base_folder,strings[i],sep = '')
  res_folder <- paste(res_folder,datasets[i],sep = '')
  res_folder <- paste(res_folder,'_',sep = '')
  for(j in 1:length(sets)){
    print(sets[i])
    res_folder2 <- paste(res_folder,sets[j],sep = '')
    res_folder2 <- paste(res_folder2,'/',sep = '')
    #res_folder2 <- paste(res_folder2,'_cfilt/',sep = '')
    for(k in 1:length(methods)){
      print(methods[k])
      rnk.file <-paste(res_folder2,'Complexes_',sep='')
      rnk.file <-paste(rnk.file,methods[k],sep = '')
      rnk.file <-paste(rnk.file,'.rnk',sep = '')
      
      output.file <- paste(res_folder2,'fgsea_res',sep='')
      output.file <- paste(output.file,methods[k],sep='')
      output.file <- paste(output.file,'.tsv',sep='')
      
      if(file.exists(rnk.file)){
        ranks <- read.table(rnk.file,
                            header=TRUE, colClasses = c("character", "numeric"))
        ranks <- setNames(ranks$Jaccard, ranks$Edges)
        fgseaRes <- fgsea(pathways, ranks)
        fname <- paste(datasets[i],sets[j],sep='_')
	fname <- paste(fname,methods[k],sep='_')
	#fname <- paste(fname,'cfilt',sep='_') 
	fname <- paste(fname,'.pdf',sep='')
	print(fname)
	pdf(fname,width=6,height=3)
	#plt <- plotEnrichment(pathways[["Bioplex"]],
        #       ranks) + labs(title="Bioplex")
	plt <- plotEnrichment(pathways[["Bioplex_HEK293T"]],ranks) + labs(title="Bioplex_HEK293T")
	print(plt)
	dev.off()
	print(fgseaRes)
        print(typeof(fgseaRes))
        fwrite(fgseaRes, file=output.file, sep="\t", sep2=c("", " ", ""))
      }

    }
  }
}

datasets<- c('HEK293T_smartseq3')
strings <- c('eCLIP_HepG2_single_FC_1_P_3_canonical')
sets <- c('RBP_RNA500','RBP_RNA1000','RBP_RNA2000','RBP_RNA3000')

methods <- c('PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR')

gmt.file <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_LAST/Bioplex_HEK293T.gmt"
#gmt.file <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_RBPs/Lang.gmt"

pathways <- gmtPathways(gmt.file)

wd <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_LAST/HEK293T_smartseq3_Jaccard/"
setwd(wd)

for(i in 1:length(datasets)){
  print(datasets[i])
  res_folder <- paste(base_folder,strings[i],sep = '')
  res_folder <- paste(res_folder,datasets[i],sep = '')
  res_folder <- paste(res_folder,'_',sep = '')
  for(j in 1:length(sets)){
    print(sets[i])
    res_folder2 <- paste(res_folder,sets[j],sep = '')
    res_folder2 <- paste(res_folder2,'/',sep = '')
    #res_folder2 <- paste(res_folder2,'_cfilt/',sep = '')
    for(k in 1:length(methods)){
      print(methods[k])
      rnk.file <-paste(res_folder2,'Complexes_',sep='')
      rnk.file <-paste(rnk.file,methods[k],sep = '')
      rnk.file <-paste(rnk.file,'.rnk',sep = '')
      
      output.file <- paste(res_folder2,'fgsea_res',sep='')
      output.file <- paste(output.file,methods[k],sep='')
      output.file <- paste(output.file,'.tsv',sep='')
      
      if(file.exists(rnk.file)){
        ranks <- read.table(rnk.file,
                            header=TRUE, colClasses = c("character", "numeric"))
        ranks <- setNames(ranks$Jaccard, ranks$Edges)
        fgseaRes <- fgsea(pathways, ranks)
        fname <- paste(datasets[i],sets[j],sep='_')
        fname <- paste(fname,methods[k],sep='_')
        #fname <- paste(fname,'cfilt',sep='_') 
        fname <- paste(fname,'.pdf',sep='')
        print(fname)
        pdf(fname,width=6,height=3)
        #plt <- plotEnrichment(pathways[["Bioplex"]],
        #       ranks) + labs(title="Bioplex")
        plt <- plotEnrichment(pathways[["Bioplex_HEK293T"]],ranks) + labs(title="Bioplex_HEK293T")
        print(plt)
        dev.off()
        print(fgseaRes)
        print(typeof(fgseaRes))
        fwrite(fgseaRes, file=output.file, sep="\t", sep2=c("", " ", ""))
      }
      
    }
  }
}

datasets<- c('HCT116')
strings <- c('eCLIP_HepG2_single_FC_1_P_3_canonical')
sets <- c('RBP_RNA500','RBP_RNA1000','RBP_RNA2000','RBP_RNA3000')

methods <- c('PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR')

gmt.file <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_LAST/Bioplex_HCT116.gmt"
#gmt.file <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_RBPs/Lang.gmt"

pathways <- gmtPathways(gmt.file)

wd <- "/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_LAST/HCT116_Jaccard/"
setwd(wd)

for(i in 1:length(datasets)){
  print(datasets[i])
  res_folder <- paste(base_folder,strings[i],sep = '')
  res_folder <- paste(res_folder,datasets[i],sep = '')
  res_folder <- paste(res_folder,'_',sep = '')
  for(j in 1:length(sets)){
    print(sets[i])
    res_folder2 <- paste(res_folder,sets[j],sep = '')
    res_folder2 <- paste(res_folder2,'/',sep = '')
    #res_folder2 <- paste(res_folder2,'_cfilt/',sep = '')
    for(k in 1:length(methods)){
      print(methods[k])
      rnk.file <-paste(res_folder2,'Complexes_',sep='')
      rnk.file <-paste(rnk.file,methods[k],sep = '')
      rnk.file <-paste(rnk.file,'.rnk',sep = '')
      
      output.file <- paste(res_folder2,'fgsea_res',sep='')
      output.file <- paste(output.file,methods[k],sep='')
      output.file <- paste(output.file,'.tsv',sep='')
      
      if(file.exists(rnk.file)){
        ranks <- read.table(rnk.file,
                            header=TRUE, colClasses = c("character", "numeric"))
        ranks <- setNames(ranks$Jaccard, ranks$Edges)
        fgseaRes <- fgsea(pathways, ranks)
        fname <- paste(datasets[i],sets[j],sep='_')
        fname <- paste(fname,methods[k],sep='_')
        #fname <- paste(fname,'cfilt',sep='_') 
        fname <- paste(fname,'.pdf',sep='')
        print(fname)
        pdf(fname,width=6,height=3)
        #plt <- plotEnrichment(pathways[["Bioplex"]],
        #       ranks) + labs(title="Bioplex")
        plt <- plotEnrichment(pathways[["Bioplex_HCT116"]],ranks) + labs(title="Bioplex_HCT116")
        print(plt)
        dev.off()
        print(fgseaRes)
        print(typeof(fgseaRes))
        fwrite(fgseaRes, file=output.file, sep="\t", sep2=c("", " ", ""))
      }
      
    }
  }
}
