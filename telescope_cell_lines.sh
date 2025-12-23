#!/bin/bash -l

#Counts the TEs using Telescope and Sort and rename the resulting file.

#notify changes via email and use all changes
#SBATCH --mail-user [YOUR EMAIL ADRESS]
#SBATCH --mail-type=ALL

#name the job
#SBATCH -J telescope


#run time
#SBATCH -t 72:0:0

#SBATCH -p tier3 -c 1

#memory requiriments in MB 
#may have to change based on what happens with run
#SBATCH --mem=250000


#these will depend on whether your IT department was willing to install telescope or you and to install it on your own instance.
#you make be using a different version on samtools depending.
spack load py-telescope-ngs/g4vuv7e
module load samtools/1.3.1

#ALL FOLDERS NEED TO BE MADE PRIOR TO RUNNING THE SCRIPT!!!!
#You can of course edit script for your own file system preference but any changes made in the Bowtie2 script need to carried forward.

INPUT=[LOCATION ON YOUR SLURM DRIVE]/DE_REs/other_cell_line/BAM_files/BAM_sorted
OUTPUT=[LOCATION ON YOUR SLURM DRIVE]/DE_REs/telescope/other_cell_lines
GTF=[LOCATION ON YOUR SLURM DRIVE]/HG19_Repeatmakser_update.gtf
UPDATE=[LOCATION ON YOUR SLURM DRIVE]/DE_REs/other_cell_line/BAM_files/Updated_BAM
SORT=[LOCATION ON YOUR SLURM DRIVE]/DE_REs/other_cell_line/BAM_files/Updated_BAM/BAM_sorted


while read FILE_NAME;
do

telescope assign $INPUT/$FILE_NAME.bam $GTF --outdir $OUTPUT/$2 --updated_sam #$2 is the second input after the script

mv $OUTPUT/$2/"telescope-telescope_report.tsv" $OUTPUT/$FILE_NAME.tsv
mv $OUTPUT/$2/"telescope-updated.bam" $UPDATE/$FILE_NAME.bam


        samtools sort -o $SORT/$FILE_NAME.bam $UPDATE/$FILE_NAME.bam
        samtools index $SORT/$FILE_NAME.bam

        
        rm -v "$UPDATE/$FILE_NAME.bam"
        mv $SORT/$FILE_NAME.bam.bai $SORT/$FILE_NAME.bai 


done < /home/jmfsbi/scripts/$1 #$1 is the first input after the script.
