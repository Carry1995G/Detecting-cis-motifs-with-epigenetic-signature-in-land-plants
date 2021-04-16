#!/bin/bash

#For chromatin_daten kann man auch bedtools intersect benutzen, allerdings muss der annotationsfile im bed format sein also weder gtf noch gff file sind o.k. 
#Kann man entweder umformatieren oder per Linux sich die Spalten holen, die man braucht.

#mkdir FeatureOverlap

#cd /netscratch/dep_coupland/grp_turck/carina/Downloads


#Annotationfile im BED format
#gtf2bed < Mp*.gtf > MpTak1v5_GTFtoBED.bed


#would like a list of H3K27me3 positive regions that are interegenic, not overlapping with any annotationsfile
#cd /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/EPIC/
#for i in ./*
#do
#bedtools intersect -a $i -b ../../MpTak1v5_GTFtoBED.bed -v | sort -k1,1 -k2,2 - > /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/$i.intergenic.H3K27me3.bed
#done

#I need a contol, which could be a bed file of randomly shouffled positions and a bed file of randomly shuffled positions also not overlapping with any transcrip annotationsfile
cd /netscratch/dep_coupland/grp_turck/carina/Downloads/Marchantia-chromatin/FeatureOverlap/

for i in *.intergenic.H3K27me3.bed
do
bedtools shuffle -excl $i -i $i -g /netscratch/dep_coupland/grp_turck/carina/MpTak1v5_chrom.txt | sort -k1,1 -k2,2 - > shuffled-notself-$i
bedtools shuffle -excl /netscratch/dep_coupland/grp_turck/carina/Downloads/MpTak1v5_GTFtoBED.bed -i $i -g /netscratch/dep_coupland/grp_turck/carina/MpTak1v5_chrom.txt | sort -k1,1 -k2,2 - > shuffled-intergenic-$i
done


#durch doe -wo function bleibt der -a Eintrag erhalten, allerdings nur für die Annotationen, die mindestens zu 50% mit H3K27me3 hits abgedeckt sind.
#Funktioniert gut für "grosse" H3K27me3 Regionen, nicht so gut für in kleine Blöcke formatieret Regionen

#Das gleiche intersect kann auch genommen werden um die H3K27me3 negativen Einträge zu behalten, allerdings mit anderer Option -v statt -wa oder -wo

#den PerlTeil hat mir jemand gesagt, der Perl als Muttersprache hatte. Hole den Gennamen, der nach ID=bis zum semicolon steht raus, dadurch sollte man eine Liste bekommen.
#um das ganze nachher in deeptools zu pluggen, wäre es vielleicht aber nicht schlecht, auch gleich einen korrekten bed file zu generieren


#########################################################################################
#for i in Annotate*.txt
#do
#cut -f4 $i | grep -o '\(A[TG0-9]*\)' | cut -c-9 | sed "1 d" |sort | uniq >$i.AGIs.txt
#done

#for i in Annotate*FDR.01.txt
#do 
#awk -F"\t" 'BEGIN {OFS="\t"}{if (/tss/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"TSS",$4;
#	else if (/5*UTR/)  print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"5UTR",$4;
#	else if (/3*UTR/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"3UTR",$4 ;
#	else if (/TTS/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"TTS",$4 ;
#	else if (/exon/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"genebody",$4;
#	else if (/intron/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"genebody",$4 ;
#	else if (/1kb-promoter/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"1kb-promoter",$4 ;
#	else if (/3kb-promoter/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"3kb-promoter",$4 ;
#	else print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"Intergenic",$4 ;}' $i >$i.final


#wc -l $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="genebody" { count++ } END { print "genebody" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="3UTR" { count++ } END { print "3UTR" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="TSS" { count++ } END { print "TSS" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="5UTR" { count++ } END { print "5UTR" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="TTS" { count++ } END { print "TTS" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="1kb-promoter" { count++ } END { print "1kb-promoter" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="3kb-promoter" { count++ } END { print "3kb-promoter" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="Intergenic" { count++ } END { print "Intergenic" , count }' $i.final >>Homer.analysis.txt
#done