###Author: Franziska Turck
##download Marchantia chromatin data from Berger lab

#link_to_publication: https://www.cell.com/current-biology/pdf/S0960-9822(19)31610-0.pdf


#get the raw data from SRA
mkdir fastq
cd fastq

for i in SRR9974612 SRR9974623 SRR9974624 SRR9974632
do 
fastq-dump --skip-technical --split-3 $i
done
cd -


