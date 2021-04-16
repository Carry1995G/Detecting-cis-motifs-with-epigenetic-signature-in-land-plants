#!/bin/bash


#SRR9974612: H3K27me3_Tak1_rep3
#SRR9974623: H3K27me3_Tak1_rep2
#SRR9974624:  H3K27me3_Tak1_rep1
#SRR9974632: H3_Tak1_rep1

################################################bamcoverage to bigwig

#cd ./BWA-align/BAM-files/


### did not work for sorted full Blacklist

#for i in *.bam                                             
#do                                                              
#samtools index $i
#bamCoverage -b $i -p 70 --normalizeUsing RPKM --ignoreDuplicates -bl /netscratch/dep_coupland/grp_turck/carina/Blacklists/sorted*.bed -o ../../Heatmaps/${i/.bam/fullblacklist.bw}
#done

#cd ../../Heatmaps 

############

#cd Heatmaps

#mv SRR9974612fullblacklist.bw H3K27me3_Tak1_rep3_fullblacklist.bw
#mv SRR9974623fullblacklist.bw H3K27me3_Tak1_rep2_fullblacklist.bw
#mv SRR9974624fullblacklist.bw H3K27me3_Tak1_rep1_fullblacklist.bw
#mv SRR9974632fullblacklist.bw H3_Tak1_rep1_fullblacklist.bw

##########################################################




###############################################################################################################################
#plot ci-element coverage for interegenic regions
cd Heatmaps
computeMatrix reference-point -S H3K27me3_Tak1_rep3_blacklist.bw \
				H3_Tak1_rep1_blacklist.bw \
				-R /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/testrep1.intergenic.H3K27me3.bed \
				/netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/shuffled-intergenic-testrep1.intergenic.H3K27me3.bed\
				/netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/shuffled-notself-testrep1.intergenic.H3K27me3.bed\
                --referencePoint center \
				-b 1000 -a 1000 \
				--missingDataAsZero \
				-bl /netscratch/dep_coupland/grp_turck/carina/Blacklists/cut*.bed \
				--skipZeros \
				-p 2 \
				-o Histone_intergenicK27+_fullblacklist_misDas0_skip0_CDS.gz

plotHeatmap 	-m /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/Heatmaps/Histone_intergenicK27+_fullblacklist_misDas0_skip0_CDS.gz \
		-o  Intergenic_pretty.png \
		--samplesLabel "H3K27me3" "H3" \
		--sortUsing region_length \
		--legendLocation "best"

computeMatrix reference-point -S /netscratch/dep_coupland/grp_turck/carina/BIGWIG/sorted_telobox_MpTak1v5.1blacklist.bw \
				/netscratch/dep_coupland/grp_turck/carina/BedtoolsShuffle/MpTak1v5_telobox_shuffle_blacklist.bw \
				-R /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/testrep1.intergenic.H3K27me3.bed \
				/netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/shuffled-intergenic-testrep1.intergenic.H3K27me3.bed\
				/netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/shuffled-notself-testrep1.intergenic.H3K27me3.bed\
                --referencePoint center \
				-b 1000 -a 1000 \
				--missingDataAsZero \
				-bl /netscratch/dep_coupland/grp_turck/carina/Blacklists/cut*.bed \
				--skipZeros \
				-p 2 \
				-o telobox_intergenicK27+_fullblacklist_misDas0_skip0_CDS.gz

plotHeatmap 	-m /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/Heatmaps/telobox_intergenicK27+_fullblacklist_misDas0_skip0_CDS.gz \
		-o  Intergenic_pretty2.png \
		--samplesLabel "telobox" "shuffled" \
				--sortUsing region_length \
		--legendLocation "best"

