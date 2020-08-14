#!/usr/bin/env bash

BASEDIR="$(pwd)"
WORKDIR="./work"

compile_blast() {
    # Fetch the sources.
    wget "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.1+-src.tar.gz"
    tar xf ncbi-blast-2.10.1+-src.tar.gz

    # Compile BLAST.
    pushd ncbi-blast-2.10.1+-src.tar.gz/c++/
    patch -p1 < "$BASEDIR/timing_info.patch"
    ./configure && make
    popd
}

compile_mermer() {
    # Fetch the sources.
    git clone https://github.com/markstein-lab/genome-enhancer

    # Compile MERMER.
    pushd genome-enhancer
    make -f make.readfasta
    cd enhancer
    make -f make.scanner

    popd
}

fetch_dm6() {
    wget "http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
    gunzip "dm6.fa.gz"
}

compile_blastdb() {
    "./ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin/makeblastdb" -in "dm6.fa" -dbtype "nucl" -out "dm6.db"
}

compile_fasta() {
    echo "dm6.fa" > fastaname.txt
    "./genome-enhancer/readfasta"
}

mkdir -p "$WORKDIR" && pushd "$WORKDIR"

fetch_dm6
compile_blast
compile_mermer
compile_blastdb
compile_fasta

popd
