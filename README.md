# Scripts related to Argonaute-bound RNA manuscript

## Prerequisites

* Julia (tested on v1.3.2);
* Julia packages:
  * ArgParse
  * BioAlignments
  * BioSequences
  * CodecZlib
  * FASTX
  * Printf
  * Statistics
* Weblogo (if pictures are needed) (tested on v3.7.4)

## Usage

Generating nucleotide frequency table in CSV format:
```bash
./scripts/fragmentation-bias.jl --rna --left-bp 0 --right-bp 21 --ref reference.fasta --bam sample.bam > sample.csv
```

Generating nucleotide frequency file (.transfac) for Weblogo:
```bash
./scripts/csv2transfac.pl --skip-header sample.csv > sample.transfac
```

Weblogo can be generated by this command:
```bash
weblogo -c classic -P "" -x "nucleotide position" -U probability -A rna -F svg -D transfac -f sample.transfac -o sample.svg
```

Showing all the available options for the script:
```bash
./scripts/fragmentation-bias.jl --help
```
