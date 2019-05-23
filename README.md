# adaptFinder
Nanopore adaptaters inference using approximate kmer count.

## Requirement
- C++ (11+)
- SeqAn (2.4.0+)

## Compiling
In order to compile this file, I recommend using the following command
~~~
g++ -std=c++14 -fopenmp  -O4 -DNDEBUG -march=native  -mtune=native  adaptFinder.cpp -lrt -o adaptFinder
~~~

## Usage
REQUIRED ARGUMENTS

    input_filename STRING

OPTIONS

    -h, --help
          Display the help message.
    -lc, --low_complexity DOUBLE
          low complexity filter threshold (for k=16), default 1.5
    -sn, --sample_n INTEGER
          sample n sequences from dataset, default 10k sequences
    -sl, --sample_length INTEGER
          size of the sampled portion, default 100 bases
    -nt, --nb_thread INTEGER
          Number of thread to work with, default is 4
    -k, --kmer_size INTEGER
          Size of the kmers, default is 16
    -lim, --limit INTEGER
          limit the number of kmer used after initial counting, default is 500
    -v, --verbosity INTEGER
          Level of details printed out (fixed for the moment)
    -e, --exact_file STRING
          path to export the exact k-mer count, if needed. Default: no export
    -o, --out_file STRING
          path to the output file, default is ./out.txt

## Example
adaptFinder file.fasta -k 16 --sample_n 20000 --sample_length 90 -nt 4 lim 1000 -e exact_out.txt -o approx_out.txt

