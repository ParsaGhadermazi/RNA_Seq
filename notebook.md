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
