library(topGO)
library(pgirmess)
setwd("~/Studium/Master/Projektmodul Turck/GO-Analyse")

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#############################TopGO#########################

geneID2GO <- (readMappings("GO-Annotation TopGO.map"))

####Gene files for TopGO: scan Gene txt files. Create vector with p-values and assign gene names to vector members 

Genes=scan(
  "genes.txt",
  quiet=TRUE,
  what=""
)

TargetGeneClusters=list.files(pattern="\\group_telobox\\.txt$", recursive=TRUE)

for (i in 1:length(TargetGeneClusters)) { 
  
  TargetGenes =scan(
    TargetGeneClusters[[i]],
    quiet=TRUE,
    what=""
  )
  
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[4])

###Gene files with p-values

#pvalue_targetGenes = rep(c(0.001),times=length(TargetGenes))
#names(pvalue_targetGenes)=c(TargetGenes)
#pvalue_allGenes = rep(c(1),times=19213)
#names(pvalue_allGenes)=c(Genes)   ###assign names to vector members

####predefined list of interesting genes
  geneNames <- names(geneID2GO)
  myInterestingGenes <- TargetGenes
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames

##store all data in sampleGOdata to faciliate access to identifiers, annotation 
##and basic data statistics

  sampleMFGOdata <- new("topGOdata", ontology = "MF",
                      allGenes = geneList, annot = annFUN.gene2GO, 
                      gene2GO=geneID2GO)
  
  sampleCCGOdata <- new("topGOdata", ontology = "CC",
                        allGenes = geneList, annot = annFUN.gene2GO, 
                        gene2GO=geneID2GO)
  
  sampleBPGOdata <- new("topGOdata", ontology = "BP",
                        allGenes = geneList, annot = annFUN.gene2GO, 
                        gene2GO=geneID2GO)

###########Fisher-Test
  resultFisMF <- runTest(sampleMFGOdata, algorithm = "classic", statistic = "fisher")
  resultFisCC <- runTest(sampleCCGOdata, algorithm = "classic", statistic = "fisher")
  resultFisBP <- runTest(sampleBPGOdata, algorithm = "classic", statistic = "fisher")
  
  ######## CreateTable######
  
  
  concepts <- list("BP"=resultFisBP, "CC"=resultFisCC, "MF"=resultFisMF)
  
  
  for (i in 1:length(concepts)) { 
  
    Bioconcept = concepts[[i]]
    ABR = names(concepts[i])
    
    if  ( ABR ==  "MF" ) {
      GOdata = sampleMFGOdata
    } else if ( ABR ==  "CC" ) {
      GOdata = sampleCCGOdata
    } else if ( ABR ==  "BP" ) {
      GOdata = sampleBPGOdata
    }
    
     ###Calculate FDR
    allGO = usedGO(object = GOdata)
    allRes <-GenTable(GOdata, classic=Bioconcept, ranksOf="classic", topNodes = length(allGO))
    allRes$classic[allRes$classic == "< 1e-30"] <- "1e-30"
    fdr=p.adjust(allRes$classic, method="fdr", n=length(allRes$classic))
    allRes$FDR=fdr
    
    ###alle GO-Terms, die signifikant sind (FDR corrected term below 0.05)
    important=allRes[which(allRes$FDR < 0.01),]
    
    #######create columns for table
    #######create columns for table
    GO_biological_process=paste0(important$Term,"(",important$GO.ID,")")
    GO_biological_process=gsub(" ", "_", GO_biological_process, fixed = TRUE)
    
    
    Fold_enrichment = important$Significant/important$Expected
    
    
    FDR=important$FDR
    Gene_number=important$Significant
    Group=rep(teloclustername,times=nrow(important))
    
    ############export table
    if (length(GO_biological_process) != length(Gene_number)) { 
      next 
    } else { 
      table=data.frame(GO_biological_process,	Gene_number, Fold_enrichment,	FDR)
    }
    # install.packages(pgirmess)
    
    write.delim(table, file = paste0(ABR,"_",teloclustername,".txt"), row.names = FALSE, quote = FALSE, sep = "\t", )
  
    ####export group table
    grouptable=data.frame(Group,GO_biological_process,	Gene_number,	Fold_enrichment,	FDR)
    
    # assign flexible variable with grouptable within loop
    assign(paste0(ABR,"_group_",teloclustername), grouptable)
  }
}

group_MF_telolike = rbind(MF_group_telolike_cluster1, MF_group_telolike_cluster2, MF_group_telolike_cluster4)

group_CC_telolike =rbind(CC_group_telolike_cluster1,CC_group_telolike_cluster2,CC_group_telolike_cluster4)

group_BP_telolike=rbind(BP_group_telolike_cluster1, BP_group_telolike_cluster4)

group_MF_telobox = rbind(MF_group_telobox_cluster1, MF_group_telobox_cluster2, MF_group_telobox_cluster4)

group_CC_telobox =rbind (CC_group_telobox_cluster1,CC_group_telobox_cluster2,CC_group_telobox_cluster4)

group_BP_telobox = rbind(BP_group_telobox_cluster1, BP_group_telobox_cluster4)


groupfiles = list("group_telobox_MF"=group_MF_telobox, 
              "group_telobox_CC"=group_CC_telobox,
              "group_telobox_BP"=group_BP_telobox,
              "group_telolike_MF"=group_MF_telolike,
              "group_telolike_CC"=group_CC_telolike,
              "group_telolike_BP"=group_BP_telolike)

for (a in 1:length(groupfiles)) {
  fullgrouptable = groupfiles[[a]]
  filename = names(groupfiles[a])

  write.delim(fullgrouptable, file = paste0(filename,".txt"), row.names = FALSE, quote = FALSE, sep = "\t" )
}



########### ViSEAGO ############################################################################################

BiocManager::install("ViSEAGO") # install package from Bioconductor
# load genes background
background<-scan( "genes.txt",quiet=TRUE, what="")
# load gene selection
selection<-scan( "telolike_cluster1.txt", quiet=TRUE,what="")
# connect to Custom file
Custom <- scan( "Go-Annotation.txt", quiet=TRUE, what="")
Custom2<-ViSEAGO::Custom2GO(system.file("extdata/Go-Annotation.txt",package = "ViSEAGO"))
select(GO.db,columns=columns(GO.db),keys=keys(GO.db))
# load GO annotations from Custom
myGENE2GO<-scan("Go-Annotation.txt", quiet=TRUE, what="")
BP<-ViSEAGO::create_topGOdata(geneSel=selection, allGenes=background, gene2GO=myGENE2GO,ont="BP", nodeSize=5)
