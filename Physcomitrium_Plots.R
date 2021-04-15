Phys_bed <- read.table("Physcomitriumpatens_counts_hitsperwindow.bed",
                       header = FALSE, sep="\t", 
                       stringsAsFactors=FALSE, 
                       quote=""
)

Phys_hits <- Phys_bed[!(Phys_bed$V4==0),]

Phys_hits


###### Hit Distribution analysis 

library(ggplot2)

plot= ggplot(Phys_hits,aes(x=V4)) +xlab("Telobox hits")+labs(title=paste("Telobox window hits distribution in Physcomitrella")) + geom_histogram() 
plot
ggsave(plot, file="Physcomitrella hits distribution.png") #, width = 20, height = 10, units = "cm")

Chromosomen=seq(1,27)

for (chrom in Chromosomen)  {
  ChromPhys_set = Phys_hits[which(Phys_hits$V1 == chrom),]
  
  plot=ggplot(ChromPhys_set,aes(x=V4)) + geom_histogram() + xlab("Telobox hits") +labs(title=paste("Distribution of telobox hits in 200 bp windows in", chrom))
  
  ggsave(plot, file=sprintf("Physcomitrella hits distribution in chrom %s.png", chrom), width = 14, height = 10, units = "cm")
}

ScaffoldPhys_set = Phys_hits[c(grep("scaffold_...", Phys_hits$V1)),]
scPlot=ggplot(ScaffoldPhys_set,aes(x=V4)) + geom_histogram() + xlab("Telobox hits") +labs(title=paste("Distribution of telobox hits in 200 bp windows in unplaced scaffolds"))
ggsave(scPlot, file="Physcomitrella hits distribution in unplaced scaffolds.png", width = 14, height = 10, units = "cm")

#########Physcomitrium patens analysis 

head(Phys_bed)
table (Phys_bed$V4>0)
max(Phys_bed$V4)
min(Phys_bed$V4)

#Plot des ganzen Genoms

plot(Phys_bed$V4,
     type='h',
     col='blue',
     ylab='Coverage',
     xlab='Region',
     xaxt='n',
     main='telobox Physcomitrium'
)

### Saving telobox Coverage plots of different chromosomes in png files

Chromosomen=seq(1,27)

for (chrom in Chromosomen)  {
  # PNG device
  png(paste('Telobox Physcomitrella',chrom,'.png'))
  
  # Code
  plot(Phys_bed$V4[which(Phys_bed$V1 == chrom)],
       type='h',
       col='blue',
       ylab='Coverage',
       xlab='Region',
       xaxt='n',
       main=paste('telobox Physcomitrella chrom', chrom)
  )
  # Close device
  dev.off()
}

##saving unplaced scaffolds

grep("scaffold", Phys_bed$V1)

png('Telobox Physcomitrium scaffolds.png')
plot(MpTak_bed$V4[grep("scaffold", Phys_bed$V1)],
     type='h',
     col='blue',
     ylab='Coverage',
     xlab='Region',
     xaxt='n',
     main=paste('telobox Physcomitrium scaffolds'))
dev.off()




remove(Index)

Index
