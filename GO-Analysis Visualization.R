#################################
# Plot gene ontology enrichment #
#################################

# Requires the package 'ggplot2' (needs to be installed first)
# Load the ggplot2 package
library(ggplot2)

# set the working directory where the tables to use are located
setwd("~/Studium/Master/Projektmodul Turck/GO-Analyse")

##find relevant telobox and telolike files for each biological process 
GOfiles=list.files(pattern="^group_telo[a-z]{3,4}_[A-Z]{2}.txt$", recursive=TRUE)

## loop over all files, save GOdomain and respective cis motif (Telobox/ Telolike) in variable for future naming
for (i in 1:length(GOfiles)) {
  splitname=strsplit(GOfiles[[i]], '[-__|.]')[[1]]
  GOdomain=splitname[3]
  Telo=splitname[2]
  
  # Import the table containing the enriched GO terms by groups
  GO_gp <- read.table(GOfiles[i],header=T,stringsAsFactors = T, fill = TRUE )
  
  # List objects and their structure contained in the dataframe 'GO_gp'
  ls.str(GO_gp)
  
  # Transform the column 'Gene_number' into a numeric variable
  GO_gp$Gene_number <- as.numeric(GO_gp$Gene_number)
  
  # Replace all the "_" by a space in the column containing the GO terms
  GO_gp$GO_biological_process <- chartr("_", " ", GO_gp$GO_biological_process)
  
  # Transform the column 'GO_biological_process' into factors
  GO_gp$GO_biological_process<-as.factor(GO_gp$GO_biological_process)
  
  # Transform FDR values by -log10('FDR values')
  GO_gp$'|log10(FDR)|' <- -(log10(GO_gp$FDR))
  
  # Change factor order
  #GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","WT_down","mutant_up","mutant_down"))
  GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))
  
  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(x = GO_biological_process, y = Fold_enrichment)) +
    geom_hline(yintercept = 1, linetype="dashed", 
               color = "azure4", size=.5)+
    geom_point(data=GO_gp,aes(x=GO_biological_process, y=Fold_enrichment,size = Gene_number, colour = `|log10(FDR)|`), alpha=.7)+
    scale_color_gradient(low="green",high="red",limits=c(0, NA))+
    coord_flip()+
    theme_bw()+
    theme(axis.ticks.length=unit(-0.1, "cm"),
          axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
          axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
          axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank(),
          legend.title.align=0.5)+
    xlab(paste(GOdomain))+
    ylab("Fold enrichment")+
    labs(color="-log10(FDR)", size="Number\nof genes")+
    facet_wrap(~Group,ncol=4)
  guides(y = guide_legend(order=2),
         colour = guide_colourbar(order=1))
  
  #export plot in png file
  ggsave (filename=paste(Telo,"_",GOdomain,".png"),width=8, height=6)
}
