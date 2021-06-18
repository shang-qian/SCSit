Introduction

    SCSit:A high-efficiency cell types identification tool for single-cell sequencing data from SPLiT-seq

    Install from executable binaries

    git clone https://github.com/shang-qian/SCSit.git
    cd SCSit/
    ./make.sh
    export PATH=$PATH:$(pwd)

    After installation, all the executable files can be found in SCSit/. The command line

    export PATH=$PATH:$(pwd)

    above is used for adding SCSit/ to the system PATH.
    
    or, SCSit is part of the Anaconda distribution and can be installed with Anaconda or Miniconda
    
    conda install ¨Cc bioconda scsit


Quick Start

    SCSit -r1 input_r1.fastq -r2 input_r2.fastq -p primer.list -b barcode.list -t 4 -o output

    In the above example, 4 CPU threads will be used.
Contact

    Shang-Qian Xie, Email: sqianxie@foxmail.com
    Mei-Wei Luan, Email: meiweiluan@163.com