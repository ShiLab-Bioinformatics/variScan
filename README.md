# variScan
__variScan__ is a bioinformatics tool developed in the Shi Lab at the Monash Biomedicine Discovery Institute for efficient mapping of sequencing reads to a reference sequence library.

Unlike conventional read aligners, variScan is specifically designed to map read pairs to the most similar reference among thousands of highly similar, short reference sequences (typically hundreds of bases in length).

## Download and installation
Binary releases of variScan are available from the GitHub repository: <https://github.com/ShiLab-Bioinformatics/variScan>. The program can be installed by simply decompressing the downloaded package. 

variScan currently supports x86-64 Linux systems and requires `/bin/bash`. If `bash` is located elsewhere on your system, update the first line in `run_variScan.sh` to point to the correct path.

## Running variScan
The only entry for running variScan is `run_variScan.sh`. 

`Usage: ./run_variScan.sh <IN:file1.fastq.gz> <IN:file2.fastq.gz> <IN:library.csv> <OUT:output.xlsx>`

The input files `file1.fastq.gz` and `file2.fastq.gz` must correspond to the first and second reads of each paired-end read, respectively. Reads must be generated using a stranded protocol (forward–reverse orientation). 

Referece sequences are provided in input file `library.csv`. This file must contain exactly two columns: the reference sequence (first column) and the sequence name (second column). Do not include a header row or column titles. All values should be plain text without quotation marks. All the reference sequences must be longer than the reads in the two fastq.gz input.

The output spreadsheet is written in `output.xlsx`. 

Below is an example of `library.csv`.
```
ATCGGTCAATCGTAGCTAATCGGTCAATCGTAGCTAATCGGTCAATCGTAGCTA,seq001
TACGGTCAATCGTAGCTAATCGGTCAATCGTAGCTAATCGGTCAATCGTAGCTT,seq002
CCCGGTCAATCGTAGCTAATCGGTCAATCGTAGCTAATCGGTCAATCGTAGCTG,seq003
......
```

The following parameters can be adjusted by editing `run_variScan.sh`:
1. `THREADS`: Number of CPU cores to use. Default: `8`.
1. `MAX_READ_LENGTH`: Maximum allowed read length in the input data. All reads must be ≤ this length. Default: `151`.
1. `MAX_MISMATCH`: Maximum number of mismatched bases permitted. At least one read in each pair must have a number of mismatches ≤ this threshold. Default: `3`.

## Read alignment rules
Each read in a pair is aligned against every reference sequence. Insertions and deletions (**indels**) are not allowed. Alignments are evaluated at all possible positions, including partial overlaps.

For each read–reference combination, the **optimal alignment position** is the position with the highest number of matched bases.

A read pair is considered **mappable** if, at its optimal alignment position, **at least one read end has three or fewer mismatches**.

For each read pair, the final alignment target is the reference sequence with the highest combined number of matched bases across both read ends, using their respective optimal alignment positions.

If multiple reference sequences are equally optimal, the read pair is reported as **unmappable**.
