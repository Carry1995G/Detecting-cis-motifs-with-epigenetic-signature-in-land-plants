#############################################################
## Plot cis motif hits per window distrubtion and coverage###
#############################################################


# Requires the package 'ggplot2' (needs to be installed first)
# Load the ggplot2 package
library(ggplot2)

# set working directory to where files used are located 
setwd("~/Studium/Master/Projektmodul Turck/GO-Analyse")


##############################

### convert hits per window BED files to data frame 

MpTak_bed <- read.table("MpTak1v5_ABswitch_hitsperwindow.bed",
                       header = FALSE, sep="\t", 
                       stringsAsFactors=FALSE, 
                       quote=""
)

Blacklist <-MpTak_bed[(MpTak_bed$V4>5),]

##Hits 
table (MpTak_bed$V4>0) #Anzahl an Positionen, wo Hits verzeichnet sind
max(MpTak_bed$V4) #Max Anzahl Hits
min(MpTak_bed$V4) #Min Anzahl Hits


###### Hit Distribution analysis 

##removing 0 hits from data
MpTak_hits <- MpTak_bed[!(MpTak_bed$V4==0),]  

## total hit distribution
plot= ggplot(MpTak_hits,aes(x=V4)) +xlab("Telobox hits")+labs(title=paste("Telobox window hits distribution  in Marchantia polymorpha")) + geom_histogram() 
plot


## hit distrubtion of each chromosome: Loop over all chromosomes separately and save in png file
Chromosomes=(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chrV"))
for (chrom in Chromosomes)  {
  ChromMpTak_set = MpTak_hits[which(MpTak_hits$V1 == chrom),]
  
  plot=ggplot(ChromMpTak_set,aes(x=V4)) + geom_histogram() + xlab("Telobox hits") +labs(title=paste("Distribution of telobox hits in 200 bp windows in", chrom))
  
  ggsave(plot, file=sprintf("MpTak Telobox hits distribution in %s.png", chrom), width = 14, height = 10, units = "cm")
}

## hit distrubtion of scaffold, saved in png file
ScaffoldMpTak_set = MpTak_hits[c(grep("unplaced-scaffold_...", MpTak_hits$V1)),]
scPlot=ggplot(ScaffoldMpTak_set,aes(x=V4)) + geom_histogram() + xlab("Telobox hits") +labs(title=paste("Distribution of telobox hits in 200 bp windows in unplaced scaffolds"))
ggsave(scPlot, file="MpTak Telobox hits distribution in unplaced scaffolds.png", width = 14, height = 10, units = "cm")


## Coverage plot of entire genome  
plot(MpTak_bed$V4,
     type='h',
     col='blue',
     ylab='Coverage',
     xlab='Region',
     xaxt='n',
     main='telobox MpTak v1.5'
)

### Saving telobox coverage plots of different chromosomes in png files

Chromosomes=(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chrV",))

for (chrom in Chromosomes)  {
  # PNG device
  png(paste('Telobox MpTak v1.5',chrom,'.png'))
  
  # Code
  plot(MpTak_bed$V4[which(MpTak_bed$V1 == chrom)],
       type='h',
       col='blue',
       ylab='Coverage',
       xlab='Region',
       xaxt='n',
       main=paste('telobox MpTak v1.5', chrom)
  )
  # Close device
  dev.off()
}

##saving unplaced scaffolds

grep("unplaced-scaffold", MpTak_bed$V1)

png('Telobox MpTak v1.5 unplaced scaffolds.png')
plot(MpTak_bed$V4[grep("unplaced-scaffold", MpTak_bed$V1)],
     type='h',
     col='blue',
     ylab='Coverage',
     xlab='Region',
     xaxt='n',
     main=paste('telobox MpTak v1.5 unplaced-scaffolds'))
dev.off()






## Schleife für das Speichern von Plots mit Regions with telobox repeats 
Chromosomen=seq(1,27)
  
for (chrom in Chromosomes)  {
  # PNG device
  png(paste('Telobox MpTak v1.5',chrom,'highrepeatsregion.png'))
  
  #prepare plot range
  data=MpTak_bed$V4[which(MpTak_bed$V1 == chrom)]
  max=max(MpTak_bed$V4[which(MpTak_bed$V1 == chrom)])
  Index=which(MpTak_bed$V4[which(MpTak_bed$V1 == chrom)] == max)
  if (Index < 11 ){
    Seq = seq(1,10,1)
  }
  elsif {
    Seq = seq((Index-10),(Index+10),1)
  }
  
  #Plot Code
  plot(data[Seq],
       type='h',
       col='blue',
       ylab='Coverage',
       xlab='Region',
       xaxt='n',
       main=paste('telobox MpTak v1.5', chrom)
  )
  
  # Close device
  dev.off()
}