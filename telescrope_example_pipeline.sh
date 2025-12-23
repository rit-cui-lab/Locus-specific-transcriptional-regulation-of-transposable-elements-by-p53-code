#!/bin/bash
# Pipeline for full anyalisis of Bowtie and Telescrope using SLURM


BOWTIE=[LOCATION TO YOUR SLURM FOLDER]/BAM_files #This needs to be filled and made befor running

#The files names without the extension for the samples being anyalized seperated into several different anyalisis files.
#The names are the example files given on GitHub and used in this anyalisis.

BUND1="cell_names_1.txt"
BUND2="cell_names_2.txt"
BUND3="cell_names_3.txt"
BUND4="cell_names_4.txt"
BUND5="cell_names_5.txt"

GTF1="hg19_1.gtf" #This is provided on GitHhub.
TELESCOPE=[LOCATION TO YOUR SLURM FOLDER}/telescope/ #This needs to be filled and made befor running

#You can chanage the name of the folders if you want to change thier names
#make thses folders before running the script!

FOLDER1="run_1"
FOLDER2="run_2"
FOLDER3="run_3"
FOLDER4="run_4"
FOLDER5="run_5"



COMMAND='sbatch --output='$BOWTIE'/bowtie1.output --error='$BOWTIE'/bowtie1.output cell_line_bowtie2.sh '$BUND1' '$GTF1''
ID1=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --output='$BOWTIE'/bowtie2.output --error='$BOWTIE'/bowtie2.output cell_line_bowtie2.sh '$BUND2' '$GTF1''
ID2=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --output='$BOWTIE'/bowtie3.output --error='$BOWTIE'/bowtie3.output cell_line_bowtie2.sh '$BUND3' '$GTF1''
ID3=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --output='$BOWTIE'/bowtie4.output --error='$BOWTIE'/bowtie4.output cell_line_bowtie2.sh '$BUND4' '$GTF1''
ID4=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --output='$BOWTIE'/bowtie5.output --error='$BOWTIE'/bowtie5.output cell_line_bowtie2.sh '$BUND5' '$GTF1''
ID5=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --dependency=afterok:'$ID1' --output='$TELESCOPE'/telescope1.output --error='$TELESCOPE'/telescope1.output telescope_cell_lines.sh '$BUND1' '$FOLDER1''
ID=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --dependency=afterok:'$ID2' --output='$TELESCOPE'/telescope2.output --error='$TELESCOPE'/telescope2.output telescope_cell_lines.sh '$BUND2' '$FOLDER2''
ID=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --dependency=afterok:'$ID3' --output='$TELESCOPE'/telescope3.output --error='$TELESCOPE'/telescope3.output telescope_cell_lines.sh '$BUND3' '$FOLDER3''
ID=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --dependency=afterok:'$ID4' --output='$TELESCOPE'/telescope4.output --error='$TELESCOPE'/telescope4.output telescope_cell_lines.sh '$BUND4' '$FOLDER4''
ID=$($COMMAND | awk ' { print $4 }')
COMMAND='sbatch --dependency=afterok:'$ID5' --output='$TELESCOPE'/telescope5.output --error='$TELESCOPE'/telescope5.output telescope_cell_lines.sh '$BUND5' '$FOLDER5''
ID=$($COMMAND | awk ' { print $4 }')
