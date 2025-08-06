# HapFIRE

HapFIRE is a python-based pipeline to calculate allele and genotype frequencies from pool-seq data given founder SNP VCF file.

Please cite the following paper for HapFIRE

* Wu, Xing, et al. "Rapid adaptation and extinction across climates in synchronized outdoor evolution experiments of Arabidopsis thaliana." bioRxiv (2025): 2025-05.

* Kessner, Darren, Thomas L. Turner, and John Novembre. "Maximum likelihood estimation of frequencies of known haplotypes from pooled sequence data." Molecular biology and evolution 30.5 (2013): 1145-1158.


## Installation
HapFIRE relies on HARP to estimate the haplotype likelihood. To make everyone's life easier, we have included a pre-compiled version of harp (harp_linux_140925_103521). If you prefer to install another version, please replace the executable in the src folder. 

Install via mamba
```
mamba env create -f environment.yml
conda activate hapfire
```

## Required input
* Phased SNP vcf file as the reference panel
  ** bealge 
* sorted and index pool seq alignment bam file
* reference genome fasta file


## Usage
The most common usage of hapfire is to estimate allele and genotype frequency from pool sequencing data

```
python3 hapFIRE.py -v example/example.recode.vcf -b example/example.bam -f example/example.fa -o test
```

If you want to estimate haplotype frequency, HapFIRE will first perform genome-wide block partition and then cluster unique haplotypes from each block, and then estimate haplotype cluster frequencies

```
python3 hapFIRE.py -v example/example.recode.vcf -b example/example.bam -f example/example.fa -s True -p bigld -o test
```
