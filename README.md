# SCSit
  A high-efficiency preprocessing tool for single-cell sequencing data from SPLiT-seq
## 1.	Introduction

# What is SCSit

  SPLiT-seq provides a low-cost platform to generate single-cell data by labeling the cellular origin of RNA through four rounds of combinatorial barcoding. Before analysing this sequence to draw biological conclusions you should always automatically and rapidly perform preprocess and quality control from SPLiT-seq, which directly identified and labeled combinatorial barcoding reads.

  SCSit can directly identify combinatorial barcodes and UMI of cell types and obtain more labeled reads, and remarkably enhance the retained data from SCS due to the exact alignment of insertion and deletion. The consistency of identified reads from SCSit increases to 97%, and mapped reads are twice than the original alignment method (e.g. BLAST and BWA). Furthermore, the runtime of SCSit is less than 10% of the original. It can accurately and rapidly analyze SPLiT-seq raw data and obtain labeled reads, as well as effectively improve the single-cell data from SPLiT-seq platform. The data and source of SCSit are available on the GitHub website https://github.com/shang-qian/SCSit.


## 2.	Basic Operations

### 2.1	Installation
  SCSit is part of the Anaconda distribution and can be installed with Anaconda or Miniconda
  
| conda install -c bioconda scsit |  
| ------------------------------- |

  Or, SCSit is part of the Anaconda distribution and can be installed with PyPi
  
| pip install scsit-tools|  
| ---------------------- |

  Or, SCSit is optimized for x86-64 CPUs. You can acquire precompiled C program from the release page with:

| git clone https://github.com/shang-qian/SCSit.git   &&   cd SCSit  && ./make.sh|
| -------------------------------------------------------------------- |

 
### 2.2	Getting started

  SCSit takes a pair of pair-end sequencing read files from SPLiT-seq in FASTQ format, barcodeList and primerList in TXT format, and outputs a pair of files in FASTQ format. SCSit only supports files in the following formats FastQ (all quality encoding variants).

Newly opened files will immediately appear in finished percentage of running files. Because of the size of these files it can take a couple of minutes to open them. SCSit operates a pair of files multithreading at a time. 

Example:

| scsit -r1 input_r1.fastq -r2 input_r2.fastq -p primer.list -b barcode.list -t 4 -o output|  
| ---------------------- |


## Required arguments:

| Parameter  | String |Description |
| ------------- |:------:|-------------|
| -r1 | R1File | Paths to files that contain input read1 of SPLiT-seq pair-end filesright foo     |
| -r2 | R2File  |  Paths to files that contain input read2 of SPLiT-seq pair-end files    |
| -p | primerList | Primer list of all oligonucleotide sequences used  |
| -b | barcodeList | The 96 well plate oligonucleotides used for each round of barcodes  |


## Optional arguments:

| Parameter  | String |Description |
| ------------- |:------:|-------------|
| -o |  output prefix | Paths to files that contain output file prefix [default: ./output]  |
| -t | int  |  Maximum number of threads to use [default: 2]   |
| -h | | Show this help message and exit  |


### 2.3	SCSit output
The raw reads FASTQ format of SPLiT-seq data as input executed SCSit program and output labeled reads with combinatorial barcodes and UMI. The output Read1 FASTQ format file was composed of combinatorial barcodes and UMI, and Read2 file was corresponding sequencing data. 

## 3.	Citation and Contact

MW, Luan. et al. SCSit: A high-efficiency preprocessing tool for single-cell sequencing data from SPLiT-seq. 


Shang-Qian Xie, Email: sqianxie@foxmail.com  
Mei-Wei Luan, Email: meiweiluan@163.com
