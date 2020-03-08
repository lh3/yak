**IMPORTANT: Since 3ace4ff, the format of binary k-mer dump is incompatible
with the previous versions. You have to rerun `yak count` to generate k-mer
dump in the new format.**

## Getting Started

```sh
# Download and compile
git clone https://github.com/lh3/yak
cd yak && make

# build k-mer hash table for assembly; count singletons
./yak count -K1.5g -t32 -o asm.yak asm.fa.gz
# build k-mer hash tables for high-coverage reads; discard singletons
./yak count -b37 -t32 -o ccs.yak ccs-reads.fq.gz
# for paired end: to provide two identical streams
./yak count -b37 -t32 -o sr.yak <(zcat sr*.fq.gz) <(zcat sr*.fq.gz)

# compute assembly or reads QV
./yak qv -t32 -p -K3.2g -l100k sr.yak asm.fa.gz > asm-sr.qv.txt
./yak qv -t32 -p sr.yak ccs-reads.fq.gz > ccs-sr.qv.txt
# compute k-mer QV for reads
./yak inspect ccs.yak sr.yak > ccs-sr.kqv.txt
# evaluate the completeness of assembly
./yak inspect sr.yak asm.yak > sr-asm.kqv.txt

# print k-mer histogram
./yak inspect sr.yak > sr.hist
# print k-mers (warning: large output)
./yak inspect -p sr.yak > sr.kmers
```

## Introduction

Yak is initially developed for two specific use cases: 1) to robustly estimate
the base accuracy of CCS reads and assembly contigs, and 2) to investigate the
systematic error rate of CCS reads. It achieves the goals by comparing
sequences to the k-mer spectrum of short reads or by comparing spectra. No
reference genome or truth data is needed.

It is worth noting that estimating base accuracy is tricky. When the accuracy
approaches Q50, both unsampled and erroneous k-mers in short reads may
interfere with a naive estimator. Yak introduces an empirical model to address
this issue. Its estimate is less affected by the coverage and the quality of
short reads.
