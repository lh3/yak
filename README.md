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
