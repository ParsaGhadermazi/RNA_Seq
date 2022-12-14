Parsa Ghadermazi

This is my notebook for RNA_Seq project.


# STEP1 (Building the folder layouts and downloading the genome)
------------------------------------------
## Downloading the raw reads |  12/3/2022

I picked the following paper for my RNA_Seq project

|Paper Name|SRA accession|
|-----|-------------|
|Dong Y, Hu J, Fan L, and Chen Q. (2016) RNA-seq-based transcriptomic and metabolomic analysis reveal stress responses and programmed cell death induced by acetic acid in Saccharomyces cerevisiae. Scientific Reports."|SRP075510|

First, I made the necessary folder layout:

```

cd /scratch/summit/$USER
mkdir RNA_Project
cd RNA_Project 

```
### Making the directory for the raw_data and scripts

```
mkdir raw_data
mkdir scripts

```

### Downloading the accession files from SRA

[Link to the SRA project](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP075510&o=acc_s%3Aa)

Downloaded the accessions from the link above and moved it to RNA_Project by jupyterhub in a file called accssions.txt to the script directory

There are two ways to download the raw reads from SRA. One is the script introduced in class I call it:

SRA_fetch_sync.sbatch

```
#!/usr/bin/env bash
 
 
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=06:00:00
#SBATCH --partition=shas
#SBATCH --mail-type=end
#SBATCH --mail-user=tstark@colostate.edu
#SBATCH --output=log-download-%j.out
 
line=accessions.txt
 
while read line
do
   echo -e $line
   echo "fasterq-dump -e $SLURM_NTASKS --split-files $line -O ../raw_data/"
   time fasterq-dump -e $SLURM_NTASKS --split-files $line -O ../raw_data/
   echo "vdb-validate $line"
   vdb-validate $line
done < $1

```

And another one provided by David, which I call it:

SRA_fetch_async.sbatch

```
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:20:00
#SBATCH --partition=shas
#SBATCH --output=log-download-%A.%a.out
# Execute code with: $ sbatch --array=0-18 download_array.sbatch
 
ARR=(SRR5832182 SRR5832183 SRR5832184 SRR5832185 SRR5832186 SRR5832187 SRR5832188 SRR5832189 SRR5832190 SRR5832191 SRR5832192 SRR5832193 SRR5832194 SRR5832195 SRR5832196 SRR5832197 SRR5832198 SRR5832199) # there are 18

line=${ARR[ $SLURM_ARRAY_TASK_ID ]}
echo -e $line
echo "fasterq-dump -e $SLURM_NTASKS --split-files $line -O ../raw_data/"
time fasterq-dump -e $SLURM_NTASKS --split-files $line -O ../raw_data/
echo "vdb-validate $line"
vdb-validate $line # this failed once for me even though the file had downloaded successfully

```

Now I activate the conda environment


```
source /curc/sw/anaconda3/latest
conda activate 2022dsci512
```

Finally I submitted the batch file 

```

sbatch SRA_fetch_async.sbatch

```

Now we can take a look at the folder structure by: 

```
export SINGULARITY_TMPDIR=/scratch/summit/$USER/.singularity/tmp
export SINGULARITY_CACHEDIR=/scratch/summit/$USER/.singularity/cache
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR
export NXF_SINGULARITY_TMPDIR=$SINGULARITY_TMPDIR
export SINGULARITYENV_TINI_SUBREAPER=1

module load singularity

singularity run -B $PWD:$PWD docker://chancsu/dockerutil:shellutil tree -L 2 $PWD

```

And here is the results 

```
singularity run -B $PWD:$PWD docker://chancsu/dockerutil:shellutil tree -L 2 $pwd
INFO:    Using cached SIF image
.
|-- notebook.txt
|-- raw_data
|   |-- SRR3567551_1.fastq
|   |-- SRR3567551_2.fastq
|   |-- SRR3567552_1.fastq
|   |-- SRR3567552_2.fastq
|   |-- SRR3567554_1.fastq
|   |-- SRR3567554_2.fastq
|   |-- SRR3567555_1.fastq
|   |-- SRR3567555_2.fastq
|   |-- SRR3567637_1.fastq
|   |-- SRR3567637_2.fastq
|   |-- SRR3567638_1.fastq
|   |-- SRR3567638_2.fastq
|   |-- SRR3567639_1.fastq
|   |-- SRR3567639_2.fastq
|   |-- SRR3567657_1.fastq
|   |-- SRR3567657_2.fastq
|   |-- SRR3567674_1.fastq
|   |-- SRR3567674_2.fastq
|   |-- SRR3567676_1.fastq
|   |-- SRR3567676_2.fastq
|   |-- SRR3567677_1.fastq
|   |-- SRR3567677_2.fastq
|   |-- SRR3567679_1.fastq
|   `-- SRR3567679_2.fastq
`-- scripts
    |-- Accessions.txt
    |-- SRA_fetch_async.sbatch
    |-- SRA_fetch_sync.sbatch
    `-- logfiles

```

***The raw data was downloaded and they were all intact!***


### Building the metadata 

I made a simple python code for creating the metadata. I created a repository on github for this project 
and the jupyter notebook for making the metadata file from SRA's metadata file is provided in the following link:

[Jupyter notebook for generating the metadata](https://github.com/ParsaGhadermazi/RNA_Seq/blob/main/Metadata.ipynb)

And here is the link to the metadata text file:

[Metadata](https://raw.githubusercontent.com/ParsaGhadermazi/RNA_Seq/main/metadata.txt)

# STEP2 (Download and index the genome)
------------------------------------------
## Downloading the genome | 12/4/2022

The second step for this project is to download and index the genome of Saccharomyces cerevisiae.

To do this I make another folder in my scratch space on summit:

```
cd /scratch/summit/$USER
mkdir SC_Genome
cd SC_Genome 

```

Before doing anything I create my path info text file:

```
touch path.txt

```

Then I download the genomes from the following links via rsync or wget:

```
 rsync -avzP https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/chromosomes/ .

```

I checked the files by: 

```
md5sum *.fa.gz  > 120422_sums.txt
diff md5sum.txt >  120422_sums.txt

```

There was no difference between the two files.

Then I unzip the files by: 

```
gunzip *.fa.gz

```
So far, the structure of the genome directory looks like the following: 

```
/scratch/summit/parsa96@colostate.edu/SC_Genome
|-- 120422_sums.txt
|-- BuildSCIndices.sbatch
|-- chrI.fa
|-- chrII.fa
|-- chrIII.fa
|-- chrIV.fa
|-- chrIX.fa
|-- chrM.fa
|-- chrV.fa
|-- chrVI.fa
|-- chrVII.fa
|-- chrVIII.fa
|-- chrX.fa
|-- chrXI.fa
|-- chrXII.fa
|-- chrXIII.fa
|-- chrXIV.fa
|-- chrXV.fa
|-- chrXVI.fa
|-- md5sum.txt
`-- path.txt

```

Next, I created the slurm script named BuildSCIndices.sbatch.
I also make the whole genome by:
```
ls -1 *.fa | sed -z 's/\n/,/g' > chr_list.txt ### To get the list of the genomes
touch BuildSCIndices.sbatch

cat *.fa > sc3_wholegenome.fa
```
I downloaded the annotation file from :

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz

gzip -d 
```




```

#!/usr/bin/env bash
 
#SBATCH --job-name=execute_hisat2-build
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=shas
#SBATCH --qos=normal    
#SBATCH --time=4:00:00   
#SBATCH --output=log_hisat2-build_%J.txt
 
# Build hisat2 indexes for C. elegans
hisat2-build -p ${SLURM_NTASKS} chrI.fa,chrII.fa,chrIII.fa,chrIV.fa,chrIX.fa,chrM.fa,chrV.fa,chrVI.fa,chrVII.fa,chrVIII.fa,chrX.fa,chrXI.fa,chrXII.fa,chrXIII.fa,chrXIV.fa,chrXV.fa,chrXVI.fa  sc3
 
# Check the build
echo -e "\n\nINDEX-BUILD: inspecting indexes:"
hisat2-inspect -s se3
 
# Capture version number
echo -e "\n\nINDEX-BUILD: version:"
hisat2-build --version

```

Now, I download a .gtf file from 


-----------------

## STEP 3: Running the pipeline 2022-12-08

Before running the pipeline I show how my project directory looks like

```

/scratch/summit/parsa96@colostate.edu/RNA_project
|-- notebook.txt
|-- outputs
|-- raw_data
|   |-- SRR3567551_1.fastq
|   |-- SRR3567551_2.fastq
|   |-- SRR3567552_1.fastq
|   |-- SRR3567552_2.fastq
|   |-- SRR3567554_1.fastq
|   |-- SRR3567554_2.fastq
|   |-- SRR3567555_1.fastq
|   |-- SRR3567555_2.fastq
|   |-- SRR3567637_1.fastq
|   |-- SRR3567637_2.fastq
|   |-- SRR3567638_1.fastq
|   |-- SRR3567638_2.fastq
|   |-- SRR3567639_1.fastq
|   |-- SRR3567639_2.fastq
|   |-- SRR3567657_1.fastq
|   |-- SRR3567657_2.fastq
|   |-- SRR3567674_1.fastq
|   |-- SRR3567674_2.fastq
|   |-- SRR3567676_1.fastq
|   |-- SRR3567676_2.fastq
|   |-- SRR3567677_1.fastq
|   |-- SRR3567677_2.fastq
|   |-- SRR3567679_1.fastq
|   |-- SRR3567679_2.fastq
|   `-- metadata.txt
`-- scripts
    |-- Accessions.txt
    |-- RNAseq_analyzer_221126.sh
    |-- RNAseq_cleanup_221126.sh
    |-- SRA_fetch_async.sbatch
    |-- SRA_fetch_sync.sbatch
    |-- execute_RNAseq_pipeline.sbatch
    `-- logfiles

```

Now I go to scripts and run execute RNAseq_pipeline.sbatch

```
sbatch RNAseq_pipeline.sbatch

```
The entire job took 03:33:13. This was extracted by :


```
sacct --starttime 2022-12-01 --format JobID,Elapsed|grep 12150281

```

Now my directories look like :

```/scratch/summit/parsa96@colostate.edu/RNA_project
|-- notebook.txt
|-- outputs
|   `-- 2022-12-08_output
|       |-- 01_fastp
|       |   |-- SC01
|       |   |-- SC010
|       |   |-- SC011
|       |   |-- SC012
|       |   |-- SC02
|       |   |-- SC03
|       |   |-- SC04
|       |   |-- SC05
|       |   |-- SC06
|       |   |-- SC07
|       |   |-- SC08
|       |   `-- SC09
|       |-- 02_hisat2
|       |   |-- SC010_summary.txt
|       |   |-- SC011_summary.txt
|       |   |-- SC012_summary.txt
|       |   |-- SC01_summary.txt
|       |   |-- SC02_summary.txt
|       |   |-- SC03_summary.txt
|       |   |-- SC04_summary.txt
|       |   |-- SC05_summary.txt
|       |   |-- SC06_summary.txt
|       |   |-- SC07_summary.txt
|       |   |-- SC08_summary.txt
|       |   `-- SC09_summary.txt
|       |--  
|       |   |-- counts.txt
|       |   `-- counts.txt.summary
|       `-- 04_samtools
|           |-- SC010_sort.bam
|           |-- SC010_sort.bam.bai
|           |-- SC010_sort.bw
|           |-- SC011_sort.bam
|           |-- SC011_sort.bam.bai
|           |-- SC011_sort.bw
|           |-- SC012_sort.bam
|           |-- SC012_sort.bam.bai
|           |-- SC012_sort.bw
|           |-- SC01_sort.bam
|           |-- SC01_sort.bam.bai
|           |-- SC01_sort.bw
|           |-- SC02_sort.bam
|           |-- SC02_sort.bam.bai
|           |-- SC02_sort.bw
|           |-- SC03_sort.bam
|           |-- SC03_sort.bam.bai
|           |-- SC03_sort.bw
|           |-- SC04_sort.bam
|           |-- SC04_sort.bam.bai
|           |-- SC04_sort.bw
|           |-- SC05_sort.bam
|           |-- SC05_sort.bam.bai
|           |-- SC05_sort.bw
|           |-- SC06_sort.bam
|           |-- SC06_sort.bam.bai
|           |-- SC06_sort.bw
|           |-- SC07_sort.bam
|           |-- SC07_sort.bam.bai
|           |-- SC07_sort.bw
|           |-- SC08_sort.bam
|           |-- SC08_sort.bam.bai
|           |-- SC08_sort.bw
|           |-- SC09_sort.bam
|           |-- SC09_sort.bam.bai
|           `-- SC09_sort.bw
|-- raw_data
|   |-- SRR3567551_1.fastq
|   |-- SRR3567551_2.fastq
|   |-- SRR3567552_1.fastq
|   |-- SRR3567552_2.fastq
|   |-- SRR3567554_1.fastq
|   |-- SRR3567554_2.fastq
|   |-- SRR3567555_1.fastq
|   |-- SRR3567555_2.fastq
|   |-- SRR3567637_1.fastq
|   |-- SRR3567637_2.fastq
|   |-- SRR3567638_1.fastq
|   |-- SRR3567638_2.fastq
|   |-- SRR3567639_1.fastq
|   |-- SRR3567639_2.fastq
|   |-- SRR3567657_1.fastq
|   |-- SRR3567657_2.fastq
|   |-- SRR3567674_1.fastq
|   |-- SRR3567674_2.fastq
|   |-- SRR3567676_1.fastq
|   |-- SRR3567676_2.fastq
|   |-- SRR3567677_1.fastq
|   |-- SRR3567677_2.fastq
|   |-- SRR3567679_1.fastq
|   |-- SRR3567679_2.fastq
|   `-- metadata.txt
`-- scripts
    |-- Accessions.txt
    |-- RNAseq_analyzer_221126.sh
    |-- RNAseq_cleanup_221126.sh
    |-- SRA_fetch_async.sbatch
    |-- SRA_fetch_sync.sbatch
    |-- execute_RNAseq_pipeline.sbatch
    |-- log_RNAseq_pipe_12150267.txt
    |-- log_RNAseq_pipe_12150281.txt
    |-- log_RNAseq_pipe_12151805.txt
    `-- logfiles
        |-- fasterq.tmp.shas0101.rc.int.colorado.edu.185939
        |-- log-download-12124629.out
        |-- log-download-12124679.out
        |-- log-download-12124699.out
        `-- log-download-12124743.out

22 directories, 89 files

```

I copied 01_fastp files to my computer for making QC chartes and counts.txt in 03_feature for making the DESeq plots.

----------------------------------
----------------------------------
Text file including important paths:

```

hisat2path="/scratch/summit/parsa96@colostate.edu/SC_Genome/sc3"

genomefa="/scratch/summit/parsa96@colostate.edu/SC_Genome/sc3_wholegenome.fa"

gtffile="/scratch/summit/parsa96@colostate.edu/SC_Genome/sacCer3.ensGene.gtf"

```