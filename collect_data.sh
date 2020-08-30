set -e

HARNESS="../target/release/harness"

export PATH="genome-enhancer/enhancer/:ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin/:$PATH"

# Vary sequence count at 30 bp
$HARNESS blast 30 1 > blast_30_1.txt
$HARNESS blast 30 3 > blast_30_3.txt
$HARNESS blast 30 5 > blast_30_5.txt
$HARNESS blast 30 10 > blast_30_10.txt

$HARNESS mermer 30 1 > mermer_30_1.txt
$HARNESS mermer 30 3 > mermer_30_3.txt
$HARNESS mermer 30 5 > mermer_30_5.txt
$HARNESS mermer 30 10 > mermer_30_10.txt

# Vary length at 1 sequence
$HARNESS blast 10 1 > blast_10_1.txt
$HARNESS blast 15 1 > blast_15_1.txt
$HARNESS blast 20 1 > blast_20_1.txt
$HARNESS blast 25 1 > blast_25_1.txt
$HARNESS blast 30 1 > blast_30_1.txt

$HARNESS mermer 10 1 > mermer_10_1.txt
$HARNESS mermer 15 1 > mermer_15_1.txt
$HARNESS mermer 20 1 > mermer_20_1.txt
$HARNESS mermer 25 1 > mermer_25_1.txt
$HARNESS mermer 30 1 > mermer_30_1.txt
