#!/bin/bash -l

#Aligns using Bowties with HG19 and then sorts and creates a BAM from the SAM orginally created.

#output and error files and job naming in pipeline script

#notify changes via email and use all changes
#SBATCH --mail-user [YOUR EMAIL ADRESS]
#SBATCH --mail-type=ALL

#name the job
#SBATCH -J bowtie_cell_lines

#output and error files

#run time
#SBATCH -t 24:0:0

#SBATCH -p tier3 -c 24 #Make sure this works with your SLURM system.

#memory requiriments in MB 
#may have to change based on what happens with run
#SBATCH --mem=80000


#What version of samtools/bowtie will depend on what is loaded on your version on SLURM/SPACK.

spack load samtools@1.16.1
spack load bowtie2@2.5.1

#ALL FOLDERS NEED TO BE MADE PRIOR TO RUNNING THE SCRIPT!!!!
#You can of course edit script for your own file system preference


#location of the transcriptome made from the GTF file
REFERENCE=[LOCATION ON YOUR SLURM DRIVE]/HG19/bowtie/hg19 
#location of the GTF file
GTF=[LOCATION ON YOUR SLURM DRIVE]/HG19/$2 #$2 is the second input after the script
#output directory
#INPUT directory
INPUT=[LOCATION ON YOUR SLURM DRIVE]/DE_REs/other_cell_line
#output directory
SORTED=[LOCATION ON YOUR SLURM DRIVE]/DE_REs/other_cell_line/BAM_files/BAM_sorted
OUTPUT=[LOCATION ON YOUR SLURM DRIVE]/other_cell_line/BAM_files


#Run while loop for each name in file
while read FILE_NAME;
do 
bowtie2 -p 24 -x $REFERENCE -k 100 -S "$OUTPUT/$FILE_NAME.sam" -U "$INPUT/$FILE_NAME.fastq.gz" 

samtools view -@ 6 -S -b "$OUTPUT/$FILE_NAME.sam" > "$OUTPUT/$FILE_NAME.bam"
samtools sort -@ 6 -o $SORTED/$FILE_NAME.bam $OUTPUT/$FILE_NAME.bam
samtools index -@ 6 $SORTED/$FILE_NAME.bam

rm -v "$OUTPUT/$FILE_NAME.sam"
rm -v "$OUTPUT/$FILE_NAME.bam"
mv $SORTED/$FILE_NAME.bam.bai $SORTED/$FILE_NAME.bai

done < [LOCATION OF NAME ON YOUR SLURM DRIVE]/$1 #$1 is the first input after the script.


